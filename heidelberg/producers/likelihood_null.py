"""The likelihood's null distribution on the record geometry (lib/likelihood.py as specified): the critical values for
the pre-registration. Per pair a background B_i (gamma about B), its chance estimate (8 shifts, overdispersion 1.4), and
background finals drawn from the background's map exp(g . r) over the box (g fitted to the tuning pairs' chance
passes by lib/likelihood.fit_gradient, as stamped in heidelberg/checks/P1_numbers.json; HD_TOY_GRAD to change it), the fit weighting the background by the same map. The null is
composite (kappa = 0 with any glint density), so F sets a glint share present and NOT suppressed in shadow (uniform);
the frozen critical value is the largest 95th percentile over the F values run. Toy t from its own generator
(cal_seed_null, 100 F, t), so the result does not depend on HD_JOBS (processes, default 2). The design (toys, seeds, the
background's spread) is in config/analysis.json, section calibration.
python -E -s -m heidelberg.producers.likelihood_null RECORD_GEOM_S.npy OUT.npy F"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from heidelberg.config import A
from heidelberg.lib import likelihood as LK
from heidelberg.stages.shadow import B_TUNING

K = A.calibration
_P1 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "checks", "P1_numbers.json")   # the tuning pairs' fit (P1)
GRAD = np.array(json.load(open(_P1))["gradient"])


def box(G):
    """The record geometry's grid (heidelberg/stages/record_geometry.py: cell centres over the test area) in half-widths."""
    h, st = A.sample.test_area_deg / 2, A.sample.shadow_grid_deg
    g = np.arange(-h + st / 2, h, st); X, Y = np.meshgrid(g, g); r = np.c_[X.ravel(), Y.ravel()] / h
    assert len(r) == G; return r


def toy(args):
    path, B0, t, F = args; Sg = np.asarray(np.load(path, mmap_mode="r")); M, G = Sg.shape; sbar = Sg.mean(1)
    rng = np.random.default_rng((K.cal_seed_null, int(round(100 * F)), t)); eta = B0 * F / (1 - F)
    w = np.exp(box(G) @ GRAD); w = w / w.mean(); p = w / w.sum()
    gs = K.cal_bg_gamma_shape; B = rng.gamma(gs, B0 / gs, M); c = rng.poisson(LK.NS * B) / LK.NS
    c = np.maximum(c + rng.normal(0, np.sqrt(LK.EXCESS * B / LK.NS + 1e-12)), 0)   # the chance estimate, overdispersed
    k = [np.r_[rng.choice(G, rng.poisson(b), p=p), rng.integers(0, G, rng.poisson(eta))].astype(int) for b in B]   # glints uniform, unsuppressed
    return LK.q_stat(c, sbar, [s[k_] for s, k_ in zip(Sg, k)], w=[w[k_] for k_ in k])[0]


if __name__ == "__main__":
    PATH, OUT, F = sys.argv[1], sys.argv[2], float(sys.argv[3]); B0, TOYS = B_TUNING, K.cal_null_toys
    with ProcessPoolExecutor(int(os.environ.get("HD_JOBS", 2))) as ex:
        q = np.array(list(ex.map(toy, [(PATH, B0, t, F) for t in range(TOYS)], chunksize=4)))
    np.save(OUT, q)
    print(f"glint share {F}: " + ", ".join(f"{pc}th percentile {np.percentile(q, pc):.2f}" for pc in (90, 95, 99)) + f"; kappa hat > 0 in {(q > 0).mean():.3f} ({TOYS} toys)")
