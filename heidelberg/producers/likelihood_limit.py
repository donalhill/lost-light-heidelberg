"""The upper limit on the glint share that the likelihood (lib/likelihood.py) gives on the record geometry, by a Neyman
belt as the shadow stage computes it (requirement T5): at each share of the stage's grid (belt_grid) toys give the
confidence-level percentile of the one-sided statistic, and on toys at the coverage shares (cal_coverage_shares) the
limit is the largest grid share not excluded: its expected value and spread with no glints, and its coverage. The
record sets are drawn as in heidelberg/producers/likelihood_power.py (background B_i gamma about B from the background's
map, its chance estimate overdispersed, glints uniform and suppressed in shadow by kappa = 1). Toy t of case k from its
own generator (cal_seed_limit, k, t); HD_JOBS processes (default 2). Then, against the same belt, the limit when every
chance estimate is off by a factor (cal_chance_scales), and,
with BLUR_PREFIX, when the glints are suppressed by the shadow at times in error (producers/blurred_geometry.py, set t
mod cal_blur_sets) while the statistic keeps the catalogue's times. The design is in config/analysis.json, section
calibration.
python -E -s -m heidelberg.producers.likelihood_limit RECORD_GEOM_S.npy [BLUR_PREFIX]"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from heidelberg.config import A
from heidelberg.lib import likelihood as LK
from heidelberg.stages.shadow import B_TUNING, BELT_F
from heidelberg.producers.likelihood_null import GRAD, box

K = A.calibration


def toy(args):
    path, B0, OD, f, k, t, *opt = args; truth, scale = (*opt, None, 1.0)[:2]                # the true geometry (blur), the chance estimates' factor
    Sg = np.asarray(np.load(path, mmap_mode="r")); St = Sg if truth is None else np.asarray(np.load(truth, mmap_mode="r")); M, G = Sg.shape; sbar = Sg.mean(1)
    rng = np.random.default_rng((K.cal_seed_limit, k, t)); w = np.exp(box(G) @ GRAD); w = w / w.mean(); p = w / w.sum()
    gs = K.cal_bg_gamma_shape; B = rng.gamma(gs, B0 / gs, M); eta = B0 * f / (1 - f) if f > 0 else 0.0
    c = rng.poisson(LK.NS * B) / LK.NS
    c = np.maximum(c + rng.normal(0, np.sqrt((OD - 1) * B / LK.NS + 1e-12)), 0)
    ks = []
    for i, s in enumerate(Sg):
        kg = rng.integers(0, G, rng.poisson(eta)); kg = kg[rng.random(len(kg)) >= St[i][kg]]
        ks.append(np.r_[rng.choice(G, rng.poisson(B[i]), p=p), kg].astype(int))
    # the stage's statistics (one-sided, two-sided), with the tuning background as the share's one reference (T5)
    return LK.ul_qt(c * scale, sbar, [s[k_] for s, k_ in zip(Sg, ks)], BELT_F, od=OD, w=[w[k_] for k_ in ks], Bref=B0)


def limits(Q, t):
    ok = Q <= t[None, :]; return np.array([BELT_F[r_].max() if r_.any() else 0.0 for r_ in ok])     # largest share not excluded (grid)


if __name__ == "__main__":
    PATH = sys.argv[1]; BLUR = sys.argv[2] if len(sys.argv) > 2 else None; B0, OD, TOYS = B_TUNING, LK.OD, K.cal_belt_toys
    with ProcessPoolExecutor(int(os.environ.get("HD_JOBS", 2))) as ex:
        t = []
        for k, f0 in enumerate(BELT_F):                           # the belt: each share's percentile, from toys at that share
            Q = np.array(list(ex.map(toy, [(PATH, B0, OD, f0, 10 + k, j) for j in range(TOYS)], chunksize=4)))
            t.append(np.percentile(Q[:, 0, k], 100 * A.statistics.confidence))
            print(f"true share {f0:.2f}: 95th percentile of q {t[-1]:.2f}", flush=True)
        t = np.array(t)
        for k, f in enumerate(K.cal_coverage_shares):             # the limit on toys at each coverage share
            Q = np.array(list(ex.map(toy, [(PATH, B0, OD, f, 30 + k, j) for j in range(TOYS)], chunksize=4)))
            ul = limits(Q[:, 0], t); cov = f"; covers in {np.mean(ul >= f - 1e-9):.3f}" if f > 0 else ""
            q16, q50, q84 = np.percentile(ul, [16, 50, 84])
            print(f"true share {f:.2f}: belt upper limit median {q50:.2f} (16-84%: {q16:.2f}-{q84:.2f}){cov} ({TOYS} toys)", flush=True)
        cases = [(f"chance scale {s_:.2f}", (None, s_)) for s_ in K.cal_chance_scales] + ([] if BLUR is None else [("timing blur", "blur")])
        for i, (name, o) in enumerate(cases):                     # against the same belt
            for k, f in enumerate(K.cal_coverage_shares):
                arg = lambda j: (f"{BLUR}_{j % K.cal_blur_sets}.npy", 1.0) if o == "blur" else o
                ul = limits(np.array(list(ex.map(toy, [(PATH, B0, OD, f, 40 + 10 * i + k, j, *arg(j)) for j in range(TOYS)], chunksize=4)))[:, 0], t)
                cov = f", covers in {np.mean(ul >= f - 1e-9):.3f}" if f > 0 else ""
                print(f"{name}, share {f:.2f}: limit median {np.median(ul):.2f}{cov} ({TOYS} toys)", flush=True)
