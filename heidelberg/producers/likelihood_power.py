"""The power of the background-constrained likelihood (lib/likelihood.py, the primary) and of the claim rule (the
likelihood AND the conditional check both rejecting at 0.05) on the record geometry, by toy Monte Carlo. Per pair a
background B_i (gamma about B) from the background's map exp(g . r) (as heidelberg/producers/likelihood_null.py), its chance
estimate (8 shifts, overdispersed), and glints of density eta per pair, uniform and suppressed in shadow by kappa.
Reports, against CRIT (the frozen critical value of q), the likelihood's rejection rate with glints present but not
suppressed (must stay at or below 0.05), and at glint shares f with kappa = 1 the likelihood's power and the claim's
(the conditional check from cal_check_draws draws per toy, its controls weighted by the map). Toy t of case k from its
own generator (cal_seed_power, k, t); HD_JOBS processes (default 2). Mode power KAPPA: the power table at that
suppression in shadow. Mode size CASE: the background drawn from the map case's multiple of the fitted gradient
(cal_map_cases) while both tests keep the fitted one, and the cases are the false-rejection ones (cal_size_shares,
glints unsuppressed), with the check's own rejection rate beside the claim's. Mode refused: the power table at kappa 1
with each pair refused at the tuning pairs' rate (pairs not selected, heidelberg/checks/P1_numbers.json), from a
generator of its own. Mode blur PREFIX: the power table at kappa 1 with the glints suppressed by the shadow at times in
error (producers/blurred_geometry.py, set t mod cal_blur_sets) while both tests keep the catalogue's times. The design
is in config/analysis.json, section calibration.
python -E -s -m heidelberg.producers.likelihood_power RECORD_GEOM_S.npy CRIT power KAPPA | size CASE | refused | blur PREFIX"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from heidelberg.config import A
from heidelberg.lib import likelihood as LK
from heidelberg.stages.shadow import B_TUNING
from heidelberg.producers.likelihood_null import GRAD, _P1, box

K = A.calibration


def toy(args):
    path, B0, OD, f, kappa, k, t, gt, *opt = args; truth, refuse = (*opt, None, 0.0)[:2]     # the true geometry (blur), the share refused
    Sg = np.asarray(np.load(path, mmap_mode="r")); St = Sg if truth is None else np.asarray(np.load(truth, mmap_mode="r"))
    if refuse:                                                    # from a generator of its own, the other draws unchanged
        keep = np.random.default_rng((K.cal_seed_power, k, t, 1)).random(len(Sg)) >= refuse; Sg, St = Sg[keep], St[keep]
    M, G = Sg.shape; sbar = Sg.mean(1)
    rng = np.random.default_rng((K.cal_seed_power, k, t)); w = np.exp(box(G) @ GRAD); w = w / w.mean(); p = w / w.sum()
    pt = p if gt is None else np.exp(box(G) @ np.asarray(gt)) / np.exp(box(G) @ np.asarray(gt)).sum()
    gs = K.cal_bg_gamma_shape; B = rng.gamma(gs, B0 / gs, M); eta = B0 * f / (1 - f) if f > 0 else 0.0
    c = rng.poisson(LK.NS * B) / LK.NS
    c = np.maximum(c + rng.normal(0, np.sqrt((OD - 1) * B / LK.NS + 1e-12)), 0)
    ks = []
    for i, s in enumerate(Sg):
        kg = rng.integers(0, G, rng.poisson(eta)); kg = kg[rng.random(len(kg)) >= kappa * St[i][kg]]
        ks.append(np.r_[rng.choice(G, rng.poisson(B[i]), p=pt), kg].astype(int))
    q = LK.q_stat(c, sbar, [s[k_] for s, k_ in zip(Sg, ks)], od=OD, w=[w[k_] for k_ in ks])[0]
    sp = Sg @ p; D = sum(len(k_) * sp[i] - Sg[i, k_].sum() for i, k_ in enumerate(ks))          # the conditional check
    nd = K.cal_check_draws; nulls = np.zeros(nd)
    for i, k_ in enumerate(ks):
        if len(k_):
            nulls += len(k_) * sp[i] - Sg[i, rng.choice(G, (nd, len(k_)), p=p)].sum(1)
    return q, (np.sum(nulls >= D) + 1) / (nd + 1)


if __name__ == "__main__":
    PATH, CRIT, MODE = sys.argv[1], float(sys.argv[2]), sys.argv[3]; B0, OD, TOYS = B_TUNING, LK.OD, K.power_toys
    GT, K0, opt = None, 0, lambda t: ()
    if MODE == "power":                                           # the suppression in shadow for the power cases
        KP = float(sys.argv[4])
        cases = [(f, 0.0) for f in K.cal_unsuppressed_shares] + [(f, KP) for f in K.cal_power_shares]
    elif MODE == "size":                                          # the size cases: the true map a multiple of the fitted gradient
        GT = tuple(dict(K.cal_map_cases)[sys.argv[4]] * GRAD)
        cases = [(f, 0.0) for f in K.cal_size_shares]
    elif MODE == "refused":
        P1 = json.load(open(_P1)); R_ = 1 - P1["selected"] / P1["pairs"]; K0, opt = 50, lambda t: (None, R_)
        cases = [(f, 1.0) for f in K.cal_power_shares]
    else:                                                         # blur: the shadow at the clock's times, the tests at the catalogue's
        K0, opt = 60, lambda t: (f"{sys.argv[4]}_{t % K.cal_blur_sets}.npy",)
        cases = [(f, 1.0) for f in K.cal_power_shares]
    with ProcessPoolExecutor(int(os.environ.get("HD_JOBS", 2))) as ex:
        for k, (f, kappa) in enumerate(cases):
            r = np.array(list(ex.map(toy, [(PATH, B0, OD, f, kappa, K0 + k, t, GT, *opt(t)) for t in range(TOYS)], chunksize=4))); lik = r[:, 0] > CRIT; chk = r[:, 1] < A.statistics.alpha
            if MODE == "refused":
                print(f"pairs refused at {R_:.3f}: share {f:.2f}, likelihood {lik.mean():.3f}, check {chk.mean():.3f}, claim {(lik & chk).mean():.3f} ({TOYS} toys)", flush=True)
            elif MODE == "blur":
                print(f"timing blur ({K.cal_blur_sets} sets): share {f:.2f}, likelihood {lik.mean():.3f}, check {chk.mean():.3f}, claim {(lik & chk).mean():.3f} ({TOYS} toys)", flush=True)
            elif kappa:
                print(f"glint share {f:.2f}, kappa {kappa:.2f}: likelihood {lik.mean():.3f}, check {chk.mean():.3f}, both (the claim) {(lik & chk).mean():.3f} ({TOYS} toys)", flush=True)
            else:
                print(f"glint share {f:.2f}, no shadow suppression: likelihood rejected {lik.mean():.3f}, check rejected {chk.mean():.3f}, both {(lik & chk).mean():.3f} ({TOYS} toys)", flush=True)
