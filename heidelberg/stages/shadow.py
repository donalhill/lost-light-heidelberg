"""The shadow stage (notes/pipeline_spec.md; notes/requirements.md T1-T7, C1): do the selected sources avoid Earth's
shadow? Per pair, each final's share of the twins' common exposure in shadow at GEO radius seen from Koenigstuhl (a
final at the midpoint of its two detections), and the same for the pair's controls: its injections found on both
plates and paired as candidates are that pass every position veto in the decision (area controls).
Two tests:
- the background-constrained likelihood (lib/likelihood.py; the primary): kappa > 0 with the glint density and every
  background profiled, each background constrained by the pair's chance estimate, the glint term weighted by the pair's
  efficiency x exposure; the p-value is the largest over glint shares (null_shares) of toys drawn on this run's geometry,
  with backgrounds about each chance estimate shrunk to the sample's mean chance estimate and the map drawn from its
  bootstrap uncertainty;
- the Monte Carlo conditioned on each pair's count (the check): expected E = sum_i n_i s_i, observed O, deficit D = E - O,
  its null from n_i map-weighted control positions per pair.
The background in both follows a log-linear map exp(g . r) fitted to the chance passes (shifts 1-7) of all pairs; glints
stay uniform. A claim needs both tests to reject at 0.05. Also reported: eta hat, the chance estimates' common factor,
the Neyman upper limit on the glint share (kappa = 1), and Bruehl et al.'s pooled binomial test as
published on random positions (all injections introduced on both plates), which carry no claim.
POP selects the finals: candidate (the test), chance or injection (the nulls, each with a pass rule: each tail of the
check >= 0.005 and the likelihood's p >= 0.01), or deficit (the positive control: the injection null's finals, each kept
with probability one minus its shadow share, a known deficit; reported with the nulls, deciding nothing). For the chance
null shift 0 plays the finals and shifts 1-7 the constraint and the map; for the injection null and the positive control
at most five passing injections per pair play the finals in both tests.
For POP candidate everything goes to OUT/sealed_*.json, which the record run time-locks.
python -E -s -m heidelberg.stages.shadow PRODUCTS_DIR SELECT_DIR LSW_PLATES.csv OUT_DIR [POP] [CONTROLS] [FRAME]
  POP candidate | chance | injection | deficit; CONTROLS area (the only one kept); FRAME topo | geo (the anti-solar
  comparison, descriptive).
  The parameters (draws, toys, grids, seeds) are in config/analysis.json, section statistics (HD_CONFIG names another
  file, for a test's reduced toys: python -m heidelberg.config override); HD_JOBS: processes (8)."""
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import binom

from heidelberg.config import A, C, DERIVED
from heidelberg.lib import likelihood as LK
from heidelberg.lib import shadow_geom as SG
from heidelberg.lib.efficiency import pair_efficiency
from heidelberg.lib.selection import DECISION, REMOTE
from heidelberg.lib.shadow_geom import POSITION, midpoint

assert set(REMOTE) <= set(POSITION)
T_ = A.statistics
P_PASS, P_TAIL = T_.null_p, T_.null_tail                   # a null fails if the likelihood's p, or either tail of the check, is below these
INJ_FINALS = T_.inj_finals                                 # injection null: passing injections per pair playing finals
NULL_F = [float(x) for x in T_.null_shares]                # glint shares over which the composite null's p-value is the largest
NULL_F_DEFAULT = ",".join(f"{x:g}" for x in NULL_F)        # as a string, for the code before the rebuild (the equivalence)
BELT_F = np.array(T_.belt_grid, float)                     # the upper limit's grid of shares
B_TUNING = DERIVED["tuning_background"]             # the tuning set's background per pair: the calibration's forecast of the sample's (the stage uses the sample's own)
BRUEHL_N = T_.bruehl_n                                     # random positions per pair in Bruehl et al.'s test



def belt_limit(grid, ok):
    """The upper limit from the belt: the largest share not excluded (ok), whether no share was excluded (the limit lies
    above the grid; the grid's largest share is a lower bound on it), and whether every share was (it lies below the grid;
    the grid's smallest share is an upper bound on it). Never 0."""
    if not ok.any():
        return float(grid[0]), False, True
    return float(grid[ok].max()), bool(ok[-1]), False


def thin_by_shadow(s, rng):
    """The positive control's thinning: each final kept with probability one minus its shadow share."""
    return rng.random(len(s)) >= np.asarray(s)


def toy_q(args):
    """One likelihood toy's statistic (run in a pool of processes)."""
    c, sbar, fin, E, ns, w = args
    return LK.q_stat(c, sbar, fin, E=E, ns=ns, w=w)[0]


def toy_ul(args):
    """One belt toy: the limit's one-sided and two-sided statistics at its own true share (index k0 of the grid)."""
    c, sbar, fin, F0, E, ns, w, k0, Bref = args
    return LK.ul_qt(c, sbar, fin, [F0[k0]], E=E, ns=ns, w=w, Bref=Bref)[:, 0]


def draw(rg, pc, k):
    """Control indices with probabilities pc."""
    return np.searchsorted(np.cumsum(pc), rg.random(k) * (1 - 1e-12))


def check(sf_all, sc_all, pc_all, rng, nmc):
    """The check (T4): per pair n_i finals with shadow shares sf, controls with shares sc drawn with probabilities pc.
    Returns D = E - O, its null draws, the per-pair expected shares, E, O, N."""
    E = O = N = 0.0; nulls = np.zeros(nmc); si_all = []
    for sf, sc, pc in zip(sf_all, sc_all, pc_all):
        n = len(sf); si = float(pc @ sc) if len(sc) else 0.0
        N += n; O += float(sf.sum()); E += n * si; si_all.append(si)
        if n and len(sc):
            nulls += n * si - sc[draw(rng, pc, (nmc, n))].sum(1)
    return E - O, nulls, si_all, E, O, N


def p_values(D, nulls):
    """One-sided p of a deficit and of an excess, from the null draws."""
    n = len(nulls)
    return float((np.sum(nulls >= D) + 1) / (n + 1)), float((np.sum(nulls <= D) + 1) / (n + 1))


def null_pass(p_deficit, p_excess, p_lik):
    """C1: the check two-sided (each tail >= P_TAIL) and the likelihood one-sided (>= P_PASS)."""
    return bool(min(p_deficit, p_excess) >= P_TAIL and p_lik >= P_PASS)


def bruehl(fin_labels, rand_labels):
    """Bruehl et al. (2026, Methods) as published: finals in shadow against the control fraction pooled over pairs, a
    one-sided binomial test of a deficit. fin_labels, rand_labels: per pair, boolean arrays (in shadow at mid-exposure)."""
    kf = int(sum(int(x.sum()) for x in fin_labels)); nf = int(sum(len(x) for x in fin_labels))
    p_pool = float(np.concatenate(rand_labels).mean()) if rand_labels and sum(len(x) for x in rand_labels) else 0.0
    return dict(in_shadow=kf, n=nf, control_fraction=p_pool, expected=nf * p_pool, p=float(binom.cdf(kf, nf, p_pool)) if nf else 1.0)


def load_pairs(PD, SD, L, POP, CTRL, rng0):
    """Per selected pair: the finals, the likelihood's finals, the controls, the chance passes for the map, the
    introduced injections' positions (Bruehl's random positions), E_i, and the common exposure."""
    pairs = []
    xy = lambda d: midpoint(d.ra.to_numpy(float), d.dec.to_numpy(float), d.ra_bo.to_numpy(float), d.dec_bo.to_numpy(float))
    for f in sorted(glob.glob(f"{SD}/*_select.json")):
        r = json.load(open(f)); p = r["pair"]
        if "refused" in r:
            continue
        T = pd.read_parquet(f"{SD}/{p}_flags.parquet")
        if POP == "candidate" and os.path.exists(f"{SD}/sealed/{p}_candidates.parquet"):      # blind mode: candidates only in the sealed store
            T = pd.concat([T, pd.read_parquet(f"{SD}/sealed/{p}_candidates.parquet")], ignore_index=True)
        fin = T[(T.population == ("injection" if POP == "deficit" else POP)) & T.passed]
        c_lik, ns = float(r["chance_final"]), LK.NS
        gch = T[(T.population == "chance") & T.passed & (T["shift"] > 0)]       # the background's map: chance passes, shifts 1-7
        if POP == "chance":                                     # the null: shift 0 plays the finals, shifts 1-7 the constraint
            fin = fin[fin["shift"] == 0]; ns = LK.NS - 1; c_lik = float(len(gch)) / ns
        fin_lik = fin
        if POP in ("injection", "deficit"):
            fin_lik = fin.iloc[rng0.permutation(len(fin))[:INJ_FINALS]]
        inj = T[T.population == "injection"]
        pos_ok = ~inj[[f"fail_{c}" for c in POSITION if c in DECISION] + ["fail_pairing"]].any(axis=1)   # paired as candidates are, then every position veto
        ctrl = inj[pos_ok]
        a, b = (L.loc[os.path.splitext(x)[0]] for x in json.load(open(f"{PD}/{p}.json"))["plates"])
        lo = max(a.dateObs - a.exposure / 172800, b.dateObs - b.exposure / 172800); hi = min(a.dateObs + a.exposure / 172800, b.dateObs + b.exposure / 172800)
        J = pd.read_parquet(f"{PD}/{p}_inj.parquet")
        intro = J[J.placed_a.astype(bool) & J.placed_b.astype(bool)]
        cx = xy(inj[pos_ok]); ra0, de0 = 0.5 * (cx[0].min() + cx[0].max()), 0.5 * (cx[1].min() + cx[1].max())
        hx, hy = 0.5 * (cx[0].max() - cx[0].min()) * np.cos(np.radians(de0)), 0.5 * (cx[1].max() - cx[1].min())
        box = lambda q: np.c_[((q[0] - ra0 + 180) % 360 - 180) * np.cos(np.radians(de0)) / hx, (q[1] - de0) / hy]   # east, north over half-widths
        pairs.append(dict(pair=p, mjd=0.5 * (lo + hi), exp=(hi - lo) * 86400, fin=xy(fin), fin_lik=xy(fin_lik), ctrl=xy(ctrl), eff=pair_efficiency(J, T), eff_bright=pair_efficiency(J, T, g_max=DERIVED["g_half_median"]),
                          intro=(intro.ra.to_numpy(float), intro.dec.to_numpy(float)), chance=c_lik, ns=ns,
                          r_lik=box(xy(fin_lik)), r_ctrl=box(xy(ctrl)), r_ch=box(xy(gch))))
    return pairs


def main():
    PD, SD, LSW, OD = sys.argv[1:5]
    POP = sys.argv[5] if len(sys.argv) > 5 else "candidate"
    CTRL = sys.argv[6] if len(sys.argv) > 6 else "area"
    assert CTRL == "area", "the efficiency controls are cut (notes/cut_list.md)"
    FRAME = sys.argv[7] if len(sys.argv) > 7 else "topo"
    os.makedirs(OD, exist_ok=True)
    L = pd.read_csv(LSW).set_index("plateID")
    pairs = load_pairs(PD, SD, L, POP, CTRL, np.random.default_rng(T_.seed_inj_finals))
    if not pairs:
        sys.exit("no selected pairs")
    NMC = T_.check_draws; jobs = int(os.environ.get("HD_JOBS", C.processes.shadow_jobs))
    sh = lambda xy, q: SG.dark(*xy, q["mjd"], q["exp"], FRAME)
    FIN = "fin_lik" if POP in ("injection", "deficit") else "fin"   # at most INJ_FINALS per pair play the finals throughout
    if POP == "deficit":                                        # the positive control: a known deficit in shadow
        rd = np.random.default_rng(T_.seed_deficit)
        for q in pairs:
            keep = thin_by_shadow(sh(q["fin_lik"], q), rd); q["fin_lik"] = tuple(np.asarray(x)[keep] for x in q["fin_lik"]); q["r_lik"] = q["r_lik"][keep]
    sc_all = [sh(q["ctrl"], q) for q in pairs]; sf_all = [sh(q[FIN], q) for q in pairs]

    # the background's map and its uncertainty (bootstrap over pairs)
    GRAD = LK.fit_gradient([q["r_ch"] for q in pairs], [q["r_ctrl"] for q in pairs])
    rb = np.random.default_rng(T_.seed_grad_boot); ix = [rb.integers(0, len(pairs), len(pairs)) for _ in range(T_.grad_boot)]
    GRAD_SD = np.std([LK.fit_gradient([pairs[i]["r_ch"] for i in b], [pairs[i]["r_ctrl"] for i in b]) for b in ix], 0) if GRAD.any() else np.zeros(2)
    bw = lambda q, r: np.exp(r @ GRAD) / np.mean(np.exp(q["r_ctrl"] @ GRAD)) if len(q["r_ctrl"]) else np.ones(len(r))
    BGW = POP not in ("injection", "deficit")                   # finals shaped like the background (candidates under the null, chance)
    pc_all = [(lambda w: w / w.sum())(bw(q, q["r_ctrl"]) if BGW else np.ones(len(q["r_ctrl"]))) if len(q["r_ctrl"]) else np.zeros(0) for q in pairs]

    # the check
    D, nulls, si_all, E, O, N = check(sf_all, sc_all, pc_all, np.random.default_rng(T_.seed_check), NMC)
    pval, pval_excess = p_values(D, nulls)

    # Bruehl et al. as published: random positions are the injections introduced on both plates, found or not
    lab = lambda x, q: SG.dark(x[0], x[1], q["mjd"], 1.0, FRAME) > 0 if len(x[0]) else np.zeros(0, bool)
    rb2 = np.random.default_rng(T_.seed_bruehl); rand = []
    for q in pairs:
        a_ = q["intro"]
        if len(a_[0]):
            ix_ = rb2.choice(len(a_[0]), min(BRUEHL_N, len(a_[0])), replace=False); rand.append(lab((a_[0][ix_], a_[1][ix_]), q))
    BRUEHL = bruehl([lab(q[FIN], q) for q in pairs], rand)

    # the likelihood and its composite null
    NTOY = T_.candidate_toys if POP == "candidate" else T_.null_toys; trng = np.random.default_rng(T_.seed_toys)
    sfl = [sh(q["fin_lik"], q) for q in pairs]
    sbar = np.array([x.mean() if len(x) else 0.0 for x in sc_all]); c_obs = np.array([q["chance"] for q in pairs], float)
    ns = pairs[0]["ns"]
    Ew = np.array([q["eff"] * q["exp"] for q in pairs], float); Ew = Ew / Ew.mean() if Ew.mean() > 0 else np.ones(len(pairs))
    wl = [bw(q, q["r_lik"]) for q in pairs]
    q_obs, k_obs = LK.q_stat(c_obs, sbar, sfl, E=Ew, ns=ns, w=wl)
    x_hat = LK.fit(c_obs, sbar, sfl, E=Ew, ns=ns, w=wl)[1]; eta_hat = float(np.exp(x_hat[0]))
    s_hat = LK.scale_hat(c_obs, sfl)                         # the chance estimates' common factor (#53)
    B_SAMPLE = float(c_obs.mean())                              # the sample's mean chance estimate, from part one (before any test): the toys'
    k_ = T_.shrink_chance_parts                                 # shrink target and the glint share's reference; each pair's background for the
    Bt = np.maximum((k_ * c_obs + B_SAMPLE) / (k_ + 1), T_.toy_b_min)   # toys is its chance, shrunk to it
    shares = NULL_F if POP == "candidate" else [0.0]
    args = []
    for F in shares:
        eta = Bt.mean() * F / (1 - F)
        for _ in range(NTOY):                                     # toys drawn here in order (reproducible), fitted in parallel
            fin_t, w_t = [], []; g_t = GRAD + trng.normal(0, 1, 2) * GRAD_SD
            for q, sc, pc, b, e, x in zip(pairs, sc_all, pc_all, Bt, Ew, sfl):
                if not len(sc):
                    fin_t.append(np.zeros(0)); w_t.append(np.zeros(0)); continue
                if POP == "candidate":                          # background from a map drawn about the fit, glints uniform and unsuppressed
                    pt = np.exp(q["r_ctrl"] @ g_t); pt = pt / pt.sum()
                    k = np.r_[draw(trng, pt, trng.poisson(b)), trng.integers(0, len(sc), trng.poisson(eta * e))].astype(int)
                else:                                           # the nulls: each pair's count kept, positions redrawn
                    k = draw(trng, pc, len(x))
                fin_t.append(sc[k]); w_t.append(bw(q, q["r_ctrl"][k]))
            c_t = np.maximum(trng.poisson(LK.NS * Bt) / LK.NS + trng.normal(0, np.sqrt(LK.EXCESS * Bt / LK.NS)), 0) if POP == "candidate" else c_obs   # the extra variance
            args.append((c_t, sbar, fin_t, Ew, ns, w_t))
    with ProcessPoolExecutor(jobs) as ex:
        qa = np.array(list(ex.map(toy_q, args, chunksize=C.processes.toy_chunk))).reshape(len(shares), NTOY)
    by_f = {str(F): dict(p=float((np.sum(q_ >= q_obs) + 1) / (NTOY + 1)), q95=float(np.percentile(q_, 95)), q99=float(np.percentile(q_, 99))) for F, q_ in zip(shares, qa)}
    p_lik = max(v["p"] for v in by_f.values()); qn = qa[int(np.argmax([v["q95"] for v in by_f.values()]))]
    LIK = dict(q=float(q_obs), kappa_hat=float(k_obs), eta_hat=eta_hat, scale_hat=s_hat, p=p_lik, toys=NTOY, null_q95=float(np.percentile(qn, 95)),
               null_q99=float(np.percentile(qn, 99)), by_glint_share=by_f)
    if POP == "candidate":
        # the upper limit on the glint share (kappa = 1), by a Neyman belt on this geometry; one reference background, the
        # sample's mean chance estimate (a constant, not rescaled by the fit, so that a share is one glint density in every
        # toy and for the observed; #53), sets the share in the toys and for the observed (T5)
        Bref = B_SAMPLE; NUL = T_.belt_toys; uargs = []
        for k0, F in enumerate(BELT_F):
            eta = Bref * F / (1 - F)
            for _ in range(NUL):
                fin_t, w_t = [], []; g_t = GRAD + trng.normal(0, 1, 2) * GRAD_SD
                for q, sc, pc, b, e in zip(pairs, sc_all, pc_all, Bt, Ew):
                    if not len(sc):
                        fin_t.append(np.zeros(0)); w_t.append(np.zeros(0)); continue
                    pt = np.exp(q["r_ctrl"] @ g_t); pt = pt / pt.sum()
                    kg = trng.integers(0, len(sc), trng.poisson(eta * e)); kg = kg[trng.random(len(kg)) >= sc[kg]]    # glints vanish in shadow
                    k = np.r_[draw(trng, pt, trng.poisson(b)), kg].astype(int); fin_t.append(sc[k]); w_t.append(bw(q, q["r_ctrl"][k]))
                c_t = np.maximum(trng.poisson(LK.NS * Bt) / LK.NS + trng.normal(0, np.sqrt(LK.EXCESS * Bt / LK.NS)), 0)
                uargs.append((c_t, sbar, fin_t, BELT_F, Ew, ns, w_t, k0, Bref))
        with ProcessPoolExecutor(jobs) as ex:
            QTt = np.array(list(ex.map(toy_ul, uargs, chunksize=C.processes.toy_chunk))).reshape(len(BELT_F), NUL, 2)
        qcrit, tcrit = (np.percentile(QTt[..., j], 100 * T_.confidence, axis=1) for j in (0, 1))
        (qo, to), s_ul = LK.ul_qt(c_obs, sbar, sfl, BELT_F, E=Ew, ns=ns, w=wl, Bref=Bref, scale=True)
        f95, above, below = belt_limit(BELT_F, qo <= qcrit)
        LIK["upper_limit"] = dict(share_95=f95, above_grid=above, below_grid=below, grid=BELT_F.tolist(), q=qo.tolist(), q_crit=qcrit.tolist(), toys=NUL, reference_background=Bref,
                                  scale_hat=s_ul, rate_per_pair=Bref * f95 / (1 - f95))
    per = [dict(pair=q["pair"], n=len(sf), n_shadow=int((sf > T_.shadow_split).sum()), n_sunlit=int((sf <= T_.shadow_split).sum()), O=float(sf.sum()), s_ctrl=si,
                chance=q["chance"], eff=q["eff"], eff_bright=q["eff_bright"], exp=q["exp"]) for q, sf, si in zip(pairs, sf_all, si_all)]
    out = dict(population=POP, controls=CTRL, frame=FRAME, pairs=len(pairs), gradient=GRAD.tolist(), gradient_sd=GRAD_SD.tolist(), N=N, null_draws=NMC,
               p0=E / N if N else 0.0, E=E, O=O, D=D, p=pval, p_excess=pval_excess,
               null_D=dict(mean=float(nulls.mean()), sd=float(nulls.std()), q05=float(np.percentile(nulls, 5)), q95=float(np.percentile(nulls, 95))),
               likelihood=LIK, bruehl_replication=BRUEHL, per_pair=per,
               pairs_chance_above_bound=int(sum(q["chance"] > LK.B_HI for q in pairs)))   # kept at the bound, and reported
    if POP == "candidate":
        json.dump(out, open(f"{OD}/sealed_{POP}_{CTRL}_{FRAME}.json", "w"), indent=1, default=float); print("sealed", POP, CTRL, FRAME, len(pairs), "pairs")
    else:
        out["null_pass"] = null_pass(pval, pval_excess, p_lik)
        json.dump(out, open(f"{OD}/null_{POP}_{CTRL}_{FRAME}.json", "w"), indent=1, default=float)
        print(json.dumps({k: v for k, v in out.items() if k != "per_pair"}, indent=1, default=float))


if __name__ == "__main__":
    main()
