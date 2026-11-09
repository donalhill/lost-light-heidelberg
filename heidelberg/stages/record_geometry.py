"""Shadow geometry of the candidate record set from GAVO metadata alone: the Bruce twin pairs (heidelberg/lib/sample.py);
per pair the shadow share on a grid of shadow_grid_deg over the test area (topo, GEO), kept if any grid point is
touched. Writes OUT.csv and OUT_s.npy (the touched pairs) and OUT_pairs.csv (every pair).
python -E -s -m heidelberg.stages.record_geometry LSW.csv OUT"""
import sys, numpy as np, pandas as pd
from heidelberg.config import A
from heidelberg.lib import shadow_geom as SG
from heidelberg.lib.sample import twin_pairs
S = A.sample
P = twin_pairs(pd.read_csv(sys.argv[1]))[["a", "b", "mjd", "exp", "ra", "dec", "year"]]
print(len(P), f"Bruce twin pairs {S.first_date[:4]} - before {S.end_before}", flush=True); P.to_csv(sys.argv[2] + "_pairs.csv", index=False)   # every pair, touched or not
h, st = S.test_area_deg / 2, S.shadow_grid_deg; g = np.arange(-h + st / 2, h, st); X, Y = np.meshgrid(g, g); X, Y = X.ravel(), Y.ravel()   # cell centres
out = []
for _, p in P.iterrows():
    dec = p.dec + Y; ra = p.ra + X / np.cos(np.radians(np.clip(dec, -89, 89)))
    s = SG.dark(ra, dec, p.mjd, p.exp)
    if s.max() > 0:
        out.append(dict(p, s_mean=float(s.mean()), partial=float(((s > S.edge_lo) & (s < S.edge_hi)).mean()), full=float((s >= S.edge_hi).mean()), s=s))
O = pd.DataFrame(out); print(len(O), "touched"); O.drop(columns="s").to_csv(sys.argv[2] + ".csv", index=False); np.save(sys.argv[2] + "_s.npy", np.vstack(O.s.values))
print("per-pair mean share in shadow: quartiles", np.percentile(O.s_mean, [25, 50, 75]).round(3), "; mean", round(O.s_mean.mean(), 3))
print(f"share of box partly shadowed ({S.edge_lo} < s < {S.edge_hi}): quartiles", np.percentile(O.partial, [25, 50, 75]).round(3), "; mean", round(O.partial.mean(), 3))
print(f"pairs with the shadow edge crossing the box (partial > {S.edge_area}):", int((O.partial > S.edge_area).sum()))
