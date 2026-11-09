"""The bound on static sources from the different-night pairs (notes/requirements.md C6): no flash can appear on both
plates of such a pair, so the sources that pass every cut beyond the pair's own chance estimate are static. With n
surviving sources against a chance estimate b summed over the pairs (taken as exact), the classical one-sided 95% upper
limit s on the static sources is the s with P(N <= n | s + b) = 0.05; per pair, s over the number of pairs.
python -E -s -m heidelberg.stages.static_bound CONTROL_SELECT_DIR OUT.json"""
import glob
import json
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.stats import poisson

from heidelberg.config import A

CL = A.statistics.confidence


def upper_limit(n, b, cl=CL):
    """The classical one-sided upper limit on a Poisson signal s over a known background b, given n observed."""
    f = lambda s: poisson.cdf(n, s + b) - (1 - cl)
    return 0.0 if f(0.0) < 0 else float(brentq(f, 0.0, 10.0 * (n + 10)))


def bound(cd):
    """(pairs selected, surviving sources, chance estimate, the limit, the limit per pair) from a control selection."""
    R = [json.load(open(f)) for f in sorted(glob.glob(f"{cd}/*_select.json"))]; ok = [r for r in R if "refused" not in r]
    n = int(sum(r.get("final", 0) for r in ok)); b = float(sum(r["chance_final"] for r in ok))
    s = upper_limit(n, b)
    return dict(pairs=len(ok), refused=len(R) - len(ok), surviving=n, chance=round(b, 3), limit_95=s, per_pair=s / max(len(ok), 1))


if __name__ == "__main__":
    out = bound(sys.argv[1]); json.dump(out, open(sys.argv[2], "w"), indent=1, default=float)
    print(f"{out['surviving']} surviving against a chance estimate of {out['chance']:.1f} on {out['pairs']} pairs: static sources below "
          f"{out['limit_95']:.2f} ({out['per_pair']:.2f} per pair) at 95% confidence")
    np.testing.assert_allclose(upper_limit(4, 4.0), 5.15, atol=0.01)       # the five tuning controls' figure
