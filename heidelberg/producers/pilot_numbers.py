"""The pilot numbers of the pre-registration (notes/preregistration.md, "Pilot data"), from the tuning pairs' primary
selection (blind: candidates' counts in SELECT_DIR/sealed/), its summary, the different-night controls' selection, and
the background map fitted on it (the shadow stage's chance null, MAP.json).
Writes OUT.json with each number and the phrase the draft states it in; requirement P1 checks that every phrase is in
the draft.
The overdispersion (the pooled ratio of the variance of the shifts' chance counts to their mean, its error a bootstrap
over pairs) is the one heidelberg/stages/derived.py measured on the same selection (heidelberg.config.DERIVED); the lag
test is the controls' chance estimate at each shift radius (config/analysis.json, section background).
python -E -s -m heidelberg.producers.pilot_numbers SELECT_DIR SUMMARY.txt CONTROL_SELECT_DIR MAP.json OUT.json"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

from heidelberg.config import A, DERIVED
from heidelberg.lib.selection import DECISION
from heidelberg.stages.static_bound import upper_limit

SD, SUMMARY, CD, MAP, OUT = sys.argv[1:6]
NS = A.background.shifts
WORDS = ['no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine']
LAGS = sorted((A.background.chance_shift_arcsec, *A.background.lag_radii_arcsec))   # the estimate's radius and the test's
listed = lambda xs: ", ".join(xs[:-1]) + f", and {xs[-1]}"
R = [json.load(open(f)) for f in sorted(glob.glob(f"{SD}/*_select.json"))]
ok = [r for r in R if "refused" not in r]
solve = [r for r in R if "solve not accepted" not in str(r.get("refused", ""))]
S = {r["pair"]: json.load(open(f"{SD}/sealed/{r['pair']}_select.json")) for r in ok}
finals = sum(S[r["pair"]]["final"] for r in ok); chance = sum(r["chance_final"] for r in ok)
finals_pl = finals + sum(S[r["pair"]]["candidate_cutflow"]["point_like"][1] for r in ok)       # the rows failing it alone
chance_pl = chance + sum(r["cutflow"]["chance"]["point_like"][1] / r["chance_dirs"] for r in ok)
G = pd.Series([r["G_eff_half"] for r in ok])
star = lost = 0                                               # star-like injections: passing, and failing only the point-like cut
for r in ok:
    T = pd.read_parquet(f"{SD}/{r['pair']}_flags.parquet")
    I = T[(T.population == "injection") & (T.dur == 0) & (T.sharpen == 1)]
    others = I[[f"fail_{c}" for c in DECISION if c != "point_like"]].any(axis=1)
    star += int(I.passed.sum()); lost += int((I.fail_point_like & ~others).sum())
txt = open(SUMMARY).read()
num = lambda pat: [float(x) for x in re.search(pat, txt).groups()]
sh, she = num(r"agreement shortfall ([+-][\d.]+) \+- ([\d.]+)"); ax, axe = num(r"minor axis offset ([+-][\d.]+) \+- ([\d.]+) dex")
nz, nze = num(r"ring-noise ratio \(paired\) ([\d.]+) \+- ([\d.]+)")
# the deciding cuts a flags table carries: the controls' selection is the dry run's, made before the pairing cut, which
# fails injections only, so its absence changes no candidate or chance verdict
decides = lambda T: [f"fail_{c}" for c in DECISION if f"fail_{c}" in T]
cf = nc = 0; cc = 0.0
for f in sorted(glob.glob(f"{CD}/*_flags.parquet")):
    T = pd.read_parquet(f); p = ~T[decides(T)].any(axis=1)
    h = T.population == "chance"; nd = int(T[h]["shift"].max() + 1)
    cf += int(p[T.population == "candidate"].sum()); cc += float(p[h].sum()) / nd; nc += 1
od, od_se = DERIVED["overdispersion"], json.load(open(f"{os.path.dirname(os.path.abspath(__file__))}/../derived.json"))["overdispersion"]["se"]
lags = {}
for lag in LAGS:
    pop = "chance" if lag == A.background.chance_shift_arcsec else f"chance_{int(lag)}"
    lags[int(lag)] = sum(float((~T[T.population == pop][decides(T)].any(axis=1)).sum()) / NS
                    for T in (pd.read_parquet(f) for f in sorted(glob.glob(f"{CD}/*_flags.parquet"))))
M = json.load(open(MAP)); (gx, gy), (sx, sy) = M["gradient"], M["gradient_sd"]
sign = lambda x, d: ("−" if x < 0 else "") + f"{abs(x):.{d}f}"
out = dict(pairs=len(R), accepted=len(solve), selected=len(ok), finals=finals, chance=round(chance, 2), finals_no_point_like=finals_pl,
           chance_no_point_like=round(chance_pl, 2), G_half_median=float(G.median()), G_half_range=[float(G.min()), float(G.max())],
           point_like_loss=round(lost / (star + lost), 4), fidelity=dict(shortfall=[sh, she], minor_axis=[ax, axe], noise=[nz, nze]),
           controls=dict(pairs=nc, finals=cf, chance=round(cc, 2)), overdispersion=[od, od_se], lags=lags, gradient=[gx, gy], gradient_sd=[sx, sy])
out["phrases"] = [
    f"Of these, {len(solve)} had both solutions accepted and {len(ok)} passed the selection",
    f"the shortfall of their agreement between plates against Gaia stars was {sh:.4f} ± {she:.4f}",
    f"the offset of their minor axes {sign(ax, 4)} ± {axe:.4f} dex", f"the ratio of their image noise to the stars' {nz:.3f} ± {nze:.3f}",
    f"at a median G of {G.median():.2f} ({G.min():.2f} to {G.max():.2f} across pairs)",
    f"which lowers their efficiency by {100 * lost / (star + lost):.1f}%",
    f"{finals} against {chance:.1f} with the point-like cut, and {finals_pl} against {chance_pl:.1f} before it",
    f"different-night pairs gave {cf} surviving sources against a chance estimate of {cc:.1f}",
    f"(a Poisson upper limit of {upper_limit(cf, cc):.1f} on the {WORDS[nc]} pairs, with the chance estimate taken as exact)",
    f"the scatter among the {WORDS[NS]} shifts is {od:.2f} ± {od_se:.2f} times that of Poisson counts",
    f"shifts of {listed([f'{x:g}' for x in LAGS])} arcseconds gave {listed([f'{v:.1f}' for v in lags.values()])} survivors",
    f"g = {sign(gx, 2)} ± {sx:.2f} east and {sign(gy, 2)} ± {sy:.2f} north, bootstrapped over {M['pairs']} pairs",
]
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
json.dump(out, open(OUT, "w"), indent=1, default=float); print(json.dumps({k: v for k, v in out.items() if k != "phrases"}, default=float))
for p in out["phrases"]:
    print(" ", p)
