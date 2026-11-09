"""The values the record run uses that are measured on the tuning set rather than chosen, from the tuning pairs' primary
selection alone: the background per pair towards which the
shadow stage's toys shrink each pair's chance estimate (the tuning pairs' chance estimate over the pairs selected), and
the chance estimate's overdispersion (the pooled ratio of the variance of the shifts' chance counts to their mean), with
a bootstrap over pairs for its error, and the median over pairs of the G at which star-like flashes are found at half the
peak efficiency (the rate's brighter range ends there). Written to OUT.json (heidelberg/derived.json, committed; the Snakefile rebuilds it
from the tuning selection and compares); the code reads it as heidelberg.config.DERIVED, the calibration and the pilot
numbers included, and the documents state it as derived.*.
python -E -s -m heidelberg.stages.derived SELECT_DIR OUT.json"""
import glob
import json
import sys

import numpy as np
import pandas as pd

from heidelberg.config import A

SD, OUT = sys.argv[1:3]
NS, ST = A.background.shifts, A.statistics
files = sorted(glob.glob(f"{SD}/*_select.json"))
ok = [r for r in (json.load(open(f)) for f in files) if "refused" not in r]
chance = sum(r["chance_final"] for r in ok)
K = np.array([np.bincount(pd.read_parquet(f"{SD}/{r['pair']}_flags.parquet", columns=["population", "passed", "shift"])
                .query("population == 'chance' and passed")["shift"].astype(int), minlength=NS)[:NS] for r in ok], float)
phi = lambda k: float(((k - k.mean(1, keepdims=True)) ** 2).sum() / ((NS - 1) * k.mean(1)).sum())
boot = np.random.default_rng(ST.overdispersion_seed).integers(0, len(K), (ST.overdispersion_boot, len(K)))
od = phi(K); od_se = float(np.std([phi(K[i]) for i in boot]))
out = {"_about": "Values measured on the tuning pairs that the record run uses (heidelberg/stages/derived.py); not to be edited by hand",
       "tuning_background": {"value": chance / len(ok), "unit": "per pair",
                             "source": f"the tuning pairs' chance estimate, {chance:.4f}, over the {len(ok)} pairs selected"},
       "overdispersion": {"value": od, "se": od_se,
                          "source": f"the pooled variance of the {NS} shifts' chance counts over their mean, on the {len(ok)} pairs selected; "
                                    f"its error from {ST.overdispersion_boot} bootstrap resamples of the pairs"},
       "g_half_median": {"value": float(pd.Series([r["G_eff_half"] for r in ok]).median()), "unit": "mag",
                         "source": f"the median over the {len(ok)} pairs selected of the G at which star-like flashes are found at half the peak efficiency"}}
json.dump(out, open(OUT, "w"), indent=1)
print(f"tuning background {chance / len(ok):.4f} per pair ({chance:.4f} over {len(ok)} pairs); overdispersion {od:.4f} +- {od_se:.4f}")
