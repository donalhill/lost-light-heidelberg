"""Stage 4, all pairs (laptop): the cut-flow with three columns (candidates, chance pairs per shift, injected glints), the
chance background and the number of finals per pair and in total, the injection fidelity pooled over pairs (the
pre-registered 2-point criterion), the efficiency by target G, duration, and sharpening, and the refused pairs with
their reasons. Shadow roles are not assigned here (blinding: the roles come from the record set, and the role report
is heidelberg/stages/role_report.py).
python -I -m heidelberg.stages.summarise PRODUCTS_DIR SELECT_DIR"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

from heidelberg.config import A
from heidelberg.lib import fidelity as FI
from heidelberg.lib.selection import CUTS

PD, SD = sys.argv[1:3]
R = [json.load(open(f)) for f in sorted(glob.glob(f"{SD}/*_select.json"))]
ok = [r for r in R if "refused" not in r]; ref = [r for r in R if "refused" in r]
print(f"{len(R)} pairs: {len(ok)} selected, {len(ref)} refused")
for r in ref:
    print(f"  refused {r['pair']}: {str(r['refused'])[:140]}")
if not ok:
    sys.exit()
BLIND = any(r.get("blind") for r in ok)          # the record run: no candidate results here (the shadow stage reads sealed/)
if BLIND:
    print("blind: candidate columns, finals, and the per-pair table are withheld (OUT_DIR/sealed/)")
rows = []
for c in CUTS:
    row = dict(cut=c)
    for pop in ("candidate", "chance", "injection"):
        fl = [r["cutflow"].get(pop, {}).get(c, [0, 0]) for r in ok]; nd = [max(r["chance_dirs"], 1) if pop == "chance" else 1 for r in ok]
        row[f"{pop}_removes"] = sum(a[0] / n for a, n in zip(fl, nd)); row[f"{pop}_alone"] = sum(a[1] / n for a, n in zip(fl, nd))
    rows.append(row)
T = pd.DataFrame(rows).round(2)
tot = {pop: sum(r["cutflow"].get(pop, {}).get("rows", 0) / (max(r["chance_dirs"], 1) if pop == "chance" else 1) for r in ok) for pop in ("candidate", "chance", "injection")}
fin = {pop: sum(r["cutflow"].get(pop, {}).get("passed", 0) / (max(r["chance_dirs"], 1) if pop == "chance" else 1) for r in ok) for pop in ("candidate", "chance", "injection")}
pd.set_option("display.width", 200)
print(("" if BLIND else f"\ncandidates {tot['candidate']:.0f} -> {fin['candidate']:.0f} final; ") + f"chance pairs per shift {tot['chance']:.1f} -> {fin['chance']:.1f}; "
      f"injected glints found on both {tot['injection']:.0f} -> {fin['injection']:.0f}")
if BLIND:
    T = T[[c for c in T.columns if not c.startswith("candidate")]]
print(T.to_string(index=False)); T.to_csv(f"{SD}/cutflow.csv", index=False)
P = pd.DataFrame([dict(pair=r["pair"], epoch=r["epoch"], candidates=r.get("candidates", np.nan), final=r.get("final", np.nan),   # blind: withheld
                       chance=r["chance_final"], injected=r["injected"],
                       injections_passed=r["injections_passed"], G_eff_half=r["G_eff_half"], source=r["source"], vetoes_complete=r["vetoes_complete"]) for r in ok])
P.to_csv(f"{SD}/pairs.csv", index=False)
rng_ = np.random.default_rng(A.statistics.summary_seed); bi = rng_.integers(0, len(P), (A.statistics.summary_boot, len(P)))   # bootstrap over pairs
ch_se = float(np.std(P.chance.to_numpy()[bi].sum(1)))
if BLIND:
    P = P.drop(columns="final"); print(f"\nchance background {P.chance.sum():.1f} +- {ch_se:.1f} (bootstrap over pairs) over {len(P)} pairs")
else:
    dx = (P.final - P.chance).to_numpy()[bi].sum(1)
    print(f"\nfinals minus chance {P.final.sum() - P.chance.sum():+.1f} +- {np.std(dx):.1f} (bootstrap over pairs, both together); chance alone +- {ch_se:.1f}")
    print(f"finals {P.final.sum()} against a chance background of {P.chance.sum():.1f} (Poisson excess {P.final.sum() - P.chance.sum():+.1f}); "
          f"finals vs chance across pairs: Spearman {P[['final', 'chance']].corr(method='spearman').iloc[0, 1]:.2f}")
# injection fidelity pooled over ALL solved pairs (not only those passing the per-pair screen), each
# criterion in equivalence form: the estimate plus 2 standard errors within the limit (the per-pair checks are screens
# for gross failure only)
F = []
for r in R:
    fp = f"{PD}/{r['pair']}_fidelity.json"                    # recomputed with the current code when present
    if os.path.exists(fp):
        F.append(json.load(open(fp)))
    elif "solve not accepted" not in str(r.get("refused", "")):
        F.append(json.load(open(f"{PD}/{r['pair']}.json")).get("injection_fidelity", {}))
F = [f for f in F if f.get("pooled")]
if F:
    P_ = [f["pooled"] for f in F]; ni, ns = sum(f["n_inj"] for f in P_), sum(f["n_star"] for f in P_)
    ai = sum(f["agree_inj"] * f["n_inj"] for f in P_ if f["n_inj"]) / ni; as_ = sum(f["agree_star"] * f["n_star"] for f in P_ if f["n_star"]) / ns
    ssh = np.sqrt(ai * (1 - ai) / ni + as_ * (1 - as_) / ns)
    def ivw(k, ke):
        v = np.array([(f[k], f[ke]) for f in F if np.isfinite(f.get(k, np.nan)) and np.isfinite(f.get(ke, np.nan)) and f[ke] > 0]).reshape(-1, 2)
        w = 1 / v[:, 1] ** 2; return (v[:, 0] * w).sum() / w.sum(), 1 / np.sqrt(w.sum())
    ax, sax = ivw("minor_axis_offset_dex", "minor_axis_se"); nz, snz = ivw("background_noise_ratio", "background_noise_se")
    k_ = A.fidelity.fid_pooled_sigmas                        # equivalence form: the estimate plus k_ standard errors within the limit
    print(f"injection fidelity pooled over {len(F)} solved pairs (equivalence form: estimate + {k_} s.e. within the limit):")
    print(f"  agreement shortfall {as_ - ai:+.4f} +- {ssh:.4f} (glints {ni}, stars {ns}) -> {'PASS' if as_ - ai + k_ * ssh <= FI.PASS_PTS else 'FAIL'} (<= {FI.PASS_PTS})")
    print(f"  minor axis offset {ax:+.4f} +- {sax:.4f} dex -> {'PASS' if abs(ax) + k_ * sax <= FI.AXIS_DEX else 'FAIL'} (|x| <= {FI.AXIS_DEX})")
    print(f"  ring-noise ratio (paired) {nz:.4f} +- {snz:.4f} -> {'PASS' if nz + k_ * snz <= FI.NOISE else 'FAIL'} (<= {FI.NOISE})")
E = pd.concat([pd.read_csv(f).assign(pair=os.path.basename(f)[:-8]) for f in glob.glob(f"{SD}/*_eff.csv")])
if len(E):
    g = E.groupby([c for c in ("dur", "sharpen", "Gbin") if c in E]).apply(lambda d: pd.Series(dict(n=d.n.sum(), eff=(d.eff * d.n).sum() / max(d.n.sum(), 1))), include_groups=False)
    g.round(3).to_csv(f"{SD}/efficiency.csv"); print("\nefficiency, star-like glints (no sharpening, no duration):"); print(g.xs((0.0, 1.0), level=("dur", "sharpen")).round(3).to_string())
