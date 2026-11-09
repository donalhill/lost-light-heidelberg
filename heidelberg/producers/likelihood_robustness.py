"""The likelihood's false-rejection rate when its background model is wrong (lib/likelihood.py fitted as specified:
the configured overdispersion, no static term, background uniform within each pair), on the record geometry, against the critical
value from the null toys of the specified model. Each case draws record sets with no glints but a background that
departs from the model in one way: the chance estimate's true overdispersion; a static term per pair that the chance
estimate cannot see (flat, as the different-night control bounds it); a background heavier at the edge of each box
(density 1 + a (r / r_max)^2), as defect maps can be.
Every case draws its background from a map (default the fitted one, as heidelberg/producers/likelihood_null.py) and the fit
weights the background by the fitted map. Toy t of case k from its own generator (cal_seed_robust, k, t); HD_JOBS processes
(default 2). The cases and the design are in config/analysis.json, section calibration.
python -E -s -m heidelberg.producers.likelihood_robustness RECORD_GEOM_S.npy CRIT"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from heidelberg.config import A
from heidelberg.lib import likelihood as LK
from heidelberg.stages.shadow import B_TUNING
from heidelberg.producers.likelihood_null import GRAD, box


def draw(rng, Sg, B0, od_true=LK.OD, static=0.0, edge=0.0, glint_f=0.0, eff_shadow=1.0, grad=None):
    M, G = Sg.shape; sbar = Sg.mean(1)
    side = int(round(np.sqrt(G))); yy, xx = np.divmod(np.arange(G), side); r2 = ((xx - (side - 1) / 2) ** 2 + (yy - (side - 1) / 2) ** 2); r2 = r2 / r2.max()
    gs = A.calibration.cal_bg_gamma_shape; B = rng.gamma(gs, B0 / gs, M)
    c = rng.poisson(LK.NS * B) / LK.NS
    c = np.maximum(c + rng.normal(0, np.sqrt(max(od_true - 1, 0) * B / LK.NS + 1e-12)), 0)
    gt = GRAD if grad is None else np.asarray(grad)                 # the true map; the fit always uses GRAD
    w = (1 + edge * r2) * np.exp(box(G) @ gt); w = w / w.sum()
    E = np.where(sbar > A.calibration.cal_shadow_heavy, eff_shadow, 1.0)                     # glint efficiency lower on shadow-heavy pairs (e.g. gegenschein)
    eta = B0 * glint_f / (1 - glint_f) if glint_f else 0.0
    fin = []
    for i, s in enumerate(Sg):
        n = rng.poisson(B[i]); ns_ = rng.poisson(static); ng = rng.poisson(eta * E[i])    # glints NOT suppressed in shadow; statics uniform
        pos = np.r_[rng.choice(G, n, p=w) if n else np.zeros(0, int), rng.integers(0, G, ns_ + ng)].astype(int); fin.append(pos)
    return c, fin, E


def toy(args):
    path, B0, k, t, kw = args; Sg = np.asarray(np.load(path, mmap_mode="r")); kw = dict(kw); wtd = kw.pop("weighted", False)
    c, ks, E = draw(np.random.default_rng((A.calibration.cal_seed_robust, k, t)), Sg, B0, **kw); wf = np.exp(box(Sg.shape[1]) @ GRAD); wf = wf / wf.mean()
    return LK.q_stat(c, Sg.mean(1), [s[k_] for s, k_ in zip(Sg, ks)], E=E if wtd else None, w=[wf[k_] for k_ in ks])[0]


def cases():
    """The robustness cases (config/analysis.json, calibration.cal_robust_cases): each a name and the departures of the
    toys' truth from the fit's model; a map case as a multiple of the fitted gradient."""
    out = []
    for name, kw in A.calibration.cal_robust_cases:
        kw = dict(kw)
        if "grad_scale" in kw:
            kw["grad"] = tuple(kw.pop("grad_scale") * GRAD)
        out.append((name, kw))
    return out


if __name__ == "__main__":
    PATH = sys.argv[1]; CRIT = float(sys.argv[2]); B0, TOYS = B_TUNING, A.calibration.power_toys
    with ProcessPoolExecutor(int(os.environ.get("HD_JOBS", 2))) as ex:
        for k, (name, kw) in enumerate(cases()):
            q = np.array(list(ex.map(toy, [(PATH, B0, k, t, tuple(kw.items())) for t in range(TOYS)], chunksize=4)))
            n = int((q > CRIT).sum()); print(f"{name:66s} rejected {n / TOYS:.3f} ({n} of {TOYS}; nominal 0.05)", flush=True)
