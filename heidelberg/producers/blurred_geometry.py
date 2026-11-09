"""The calibration's geometry with the plate times in error as the asteroid clock found them (requirement J3, the timing
blur): cal_blur_sets sets, each moving every record pair's mid-exposure (both plates of a twin pair share the mount's
clock) by an offset drawn with replacement from the clock's per-plate offsets on the tuning plates less their median,
which the record run's own clock tests (N6). The shadow share on the grid as stages/record_geometry.py computes it, in
the record set's order with the examined pairs left out (stages/geometry_set.py); sunlit pairs are recomputed too. Set
k from its own generator (cal_seed_blur, k), so the sets are computed at once in HD_JOBS processes (default 2) with the
same result; with every offset zero a set is the record geometry set.
python -E -s -m heidelberg.producers.blurred_geometry RECORD_SET.csv CLOCK.json OUT_PREFIX"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from heidelberg.config import A
from heidelberg.lib import shadow_geom as SG

K = A.calibration
h, st = A.sample.test_area_deg / 2, A.sample.shadow_grid_deg                            # the record geometry's grid (cell centres)
g = np.arange(-h + st / 2, h, st); X, Y = np.meshgrid(g, g); X, Y = X.ravel(), Y.ravel()


def shares(R, shift_min):
    """Each pair's shadow share on the grid with its time moved by shift_min (minutes, one per pair)."""
    out = []
    for (_, p), dt in zip(R.iterrows(), shift_min):
        dec = p.dec + Y; ra = p.ra + X / np.cos(np.radians(np.clip(dec, -89, 89)))
        out.append(SG.dark(ra, dec, p.mjd, p.exp, shift_days=dt / 1440))
    return np.vstack(out)


def offsets(clock):
    v = np.array(list(json.load(open(clock))["by_plate"].values()), float)
    return v - np.median(v)


def one(args):
    """Set k: its offsets, its geometry saved to OUT_k.npy, and its line for the log."""
    RS, CLOCK, OUT, k = args
    R = pd.read_csv(RS); R = R[~R.examined]; V = offsets(CLOCK)
    dt = np.random.default_rng((K.cal_seed_blur, k)).choice(V, len(R))
    M = shares(R, dt); np.save(f"{OUT}_{k}.npy", M)
    return (f"blurred set {k}: offsets rms {np.sqrt(np.mean(dt ** 2)):.2f} min; pairs touched {int((M.max(1) > 0).sum())} "
            f"(sunlit now touched {int(((M.max(1) > 0) & (R.role == 'sunlit').to_numpy()).sum())})")


if __name__ == "__main__":
    RS, CLOCK, OUT = sys.argv[1:4]
    with ProcessPoolExecutor(int(os.environ.get("HD_JOBS", 2))) as ex:
        for line in ex.map(one, [(RS, CLOCK, OUT, k) for k in range(K.cal_blur_sets)]):
            print(line, flush=True)
