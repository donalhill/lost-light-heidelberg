"""The upper limit as a rate after the unlock (notes/requirements.md U5), by code frozen before it: glints per square
degree per hour for star-like flashes of star-equivalent G 11 to 17.3, weighted uniformly (the injected range), and for
flashes brighter than the tuning pairs' median half-efficiency (G_BRIGHT), with each pair's efficiency for them. The
likelihood's glint term on pair i is eta E_i with E_i = eff_i T_i / mean(eff T) over the selected pairs (shadow stage),
so a rate R gives R A_i T_i eff_i glints there and R = eta / (A mean(eff T)), A the pairs' mean area (coverage map) and T
the common exposure in hours. The limit's eta is the stage's rate_per_pair (the share limit at its reference
background); a limit above the belt's grid makes the rate a lower bound on the limit, and one below it an upper bound
(the grid's smallest share), never 0. The reciprocity bracket: a brief
flash records 0.5 to 2 magnitudes brighter than a star of the same energy, so the range in flash energy is G 11.5-13
to 17.8-19.3 (star-equivalent magnitude plus the bracket), carried as a systematic.
python -E -s -m heidelberg.stages.rate SEALED_PRIMARY.json COVERAGE_MAP.csv OUT.json"""
import json
import sys

import numpy as np
import pandas as pd

from heidelberg.config import A, DERIVED

G_RANGE = (A.injections.g_lo, A.injections.g_hi)            # the injections' star-equivalent G
BRACKET = (A.rate.bracket_lo, A.rate.bracket_hi)            # magnitudes by which a brief flash records brighter
G_BRIGHT = (A.injections.g_lo, DERIVED["g_half_median"])     # the brighter flashes: to the tuning pairs' median half-efficiency


def rate(eta, eff, exp_s, area):
    """The rate per square degree per hour that gives eta glints on a pair of mean efficiency x exposure."""
    return float(eta / (np.mean(area) * np.mean(np.asarray(eff) * np.asarray(exp_s) / 3600.0)))


if __name__ == "__main__":
    SP, COV, OUT = sys.argv[1:4]
    o = json.load(open(SP)); U = o["likelihood"]["upper_limit"]; per = pd.DataFrame(o["per_pair"])
    C = pd.read_csv(COV).set_index("pair"); area = C.loc[per.pair, "area_deg2"].to_numpy(float)
    R = rate(U["rate_per_pair"], per.eff, per.exp, area)
    R_b = rate(U["rate_per_pair"], per.eff_bright, per.exp, area)   # the same glints, all brighter than G_BRIGHT's end
    out = dict(rate_95=R, bound="lower bound on the limit (above the grid)" if U.get("above_grid") else "upper bound on the limit (below the grid)" if U.get("below_grid") else "upper limit", share_95=U["share_95"],
               eta_95=U["rate_per_pair"], mean_area_deg2=float(np.mean(area)), mean_eff_exposure_h=float(np.mean(per.eff * per.exp / 3600)),
               rate_95_bright=R_b, bright_G=list(G_BRIGHT), mean_eff_bright=float(np.mean(per.eff_bright)),
               pairs=len(per), star_equivalent_G=list(G_RANGE), flash_energy_G=[[G_RANGE[0] + b, G_RANGE[1] + b] for b in BRACKET])
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print(f"95% {out['bound']}: {R:.3g} glints per square degree per hour (star-equivalent G {G_RANGE[0]}-{G_RANGE[1]}; flash energy "
          f"G {G_RANGE[0] + BRACKET[0]}-{G_RANGE[1] + BRACKET[0]} to {G_RANGE[0] + BRACKET[1]}-{G_RANGE[1] + BRACKET[1]})")
