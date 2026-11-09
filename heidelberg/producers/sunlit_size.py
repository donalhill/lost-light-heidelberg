"""How the number of sunlit pairs changes what the record run can detect (requirement J6): for each number of sunlit pairs
per shadow pair (cal_sunlit_cases), the record geometry's shadow pairs plus that many pairs with no shadow (a sunlit
pair's test area has shadow share zero), with toys drawn as the calibration draws them (producers/likelihood_power.py and
likelihood_limit.py): the critical value of q from toys with unsuppressed glints at the calibration's critical share
(heidelberg/checks/J3_numbers.json), the power of each test and of the claim at three glint shares (kappa = 1), and the
expected upper limit with no glints from a belt on that geometry. The record set's rule (sample.sunlit_per_shadow) is the
smallest case that reaches the best expected limit.
python -E -s -m heidelberg.producers.sunlit_size RECORD_GEOM_SET.npy J3_NUMBERS.json OUT.json"""
import json
import os
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from heidelberg.config import A
from heidelberg.lib import likelihood as LK
from heidelberg.producers import likelihood_limit as LL
from heidelberg.producers import likelihood_power as LP
from heidelberg.stages.shadow import B_TUNING, BELT_F

K = A.calibration
SHARES = tuple(K.cal_power_shares)
key = lambda x: f"{x:g}".replace(".", "p")

if __name__ == "__main__":
    PATH, J3N, OUT = sys.argv[1:4]
    Sg = np.load(PATH); shadow = Sg[Sg.max(1) > 0]; B0, OD = B_TUNING, LK.OD
    crit_share = json.load(open(J3N))["critical_share"]
    res = dict(shadow=len(shadow), critical_share=crit_share, cases={})
    with tempfile.TemporaryDirectory() as tmp, ProcessPoolExecutor(int(os.environ.get("HD_JOBS", 2))) as ex:
        for c in K.cal_sunlit_cases:
            n = round(c * len(shadow)); geo = f"{tmp}/geometry_{key(c)}.npy"
            np.save(geo, np.vstack([shadow, np.zeros((n, Sg.shape[1]))]))
            q0 = np.array([r[0] for r in ex.map(LP.toy, [(geo, B0, OD, crit_share, 0.0, 100, t, None) for t in range(K.cal_sunlit_null_toys)], chunksize=8)])
            crit = float(np.percentile(q0, 100 * (1 - A.statistics.alpha))); row = dict(per_shadow=c, sunlit=n, critical_value=crit)
            for k, f in enumerate(SHARES):
                r = np.array(list(ex.map(LP.toy, [(geo, B0, OD, f, 1.0, 110 + k, t, None) for t in range(K.cal_sunlit_toys)], chunksize=8)))
                lik, chk = r[:, 0] > crit, r[:, 1] < A.statistics.alpha
                row[f"power_{key(f)}"] = dict(likelihood=float(lik.mean()), check=float(chk.mean()), claim=float((lik & chk).mean()))
            pct = []                                              # the belt on this geometry, then the limit with no glints
            for k, f0 in enumerate(BELT_F):
                Q = np.array(list(ex.map(LL.toy, [(geo, B0, OD, f0, 210 + k, j) for j in range(K.cal_sunlit_toys)], chunksize=8)))[:, 0]
                pct.append(np.percentile(Q[:, k], 100 * A.statistics.confidence))
            Q = np.array(list(ex.map(LL.toy, [(geo, B0, OD, 0.0, 299, j) for j in range(K.cal_sunlit_toys)], chunksize=8)))[:, 0]
            ok = Q <= np.array(pct)[None, :]; ul = np.array([BELT_F[r_].max() if r_.any() else 0.0 for r_ in ok])
            row["limit"] = dict(zip(("lo16", "median", "hi84"), map(float, np.percentile(ul, [16, 50, 84]))))
            res["cases"][key(c)] = row
            print(f"sunlit {c:g} per shadow pair ({n}): critical q {crit:.2f}; claim power " + ", ".join(
                f"{f:.2f}: {row[f'power_{key(f)}']['claim']:.3f}" for f in SHARES) + f"; limit median {row['limit']['median']:.2f}", flush=True)
    best = min(r["limit"]["median"] for r in res["cases"].values())
    res["rule"] = min(r["per_shadow"] for r in res["cases"].values() if r["limit"]["median"] <= best)
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"the smallest number of sunlit pairs per shadow pair with the best expected limit ({best:.2f}): {res['rule']:g}")
