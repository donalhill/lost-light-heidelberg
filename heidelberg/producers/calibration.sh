#!/bin/bash
# J3 (notes/requirements.md): the calibration with the final likelihood (Newton profile; the belt with one reference
# background) on the record geometry: the null over the ten glint shares and its sup critical value, the power table
# (kappa 1 and 0.5; with pairs refused at the tuning pairs' rate; with the plate times in error as the asteroid clock
# found them, N5), the claim rule's and the within-pair test's false rejection with the map correct and wrong, the
# likelihood's robustness cases, and the belt's expected limit and coverage (also with the chance estimates scaled and
# with the times in error); its numbers to heidelberg/checks/J3_numbers.json. Run through
# python -m heidelberg.tests.requirements --stamp J3 bash heidelberg/producers/calibration.sh (about 50 minutes on 8 processes).
# The design (shares, toys, seeds, cases) is in config/analysis.json (sections statistics and calibration).
set -euo pipefail; cd "$(dirname "$0")/../.."
P=${PY:-~/.venvs/poss/bin/python}; G=work/heidelberg_record_geometry_set.npy; export HD_JOBS=${HD_JOBS:-6}; Q=${CAL_DIR:-work/j3}; mkdir -p $Q   # CAL_DIR, CAL_OUT: elsewhere, for a test's reduced run
{
B=$($P -E -s -c "from heidelberg.stages.shadow import B_TUNING; print(B_TUNING)")   # the tuning background (one definition)
echo "tuning background per pair: $B"
echo "processes: HD_JOBS=$HD_JOBS, nproc $(nproc --all), $(date -u +%H:%M:%S)"   # the parallelism each producer gets
$P -E -s -u -m heidelberg.producers.blurred_geometry work/record_set.csv heidelberg/checks/N5_clock.json $Q/blur   # the times in error
cfg() { $P -E -s -c "from heidelberg.config import A; $1"; }                          # the design: config/analysis.json
for F in $(cfg "print(' '.join(f'{x:g}' for x in A.statistics.null_shares))"); do
  $P -E -s -u -m heidelberg.producers.likelihood_null $G $Q/null_q_f$F.npy $F        # calibration.cal_null_toys toys per share
done
C=$($P -c "import numpy as np,glob; print(round(max(float(np.percentile(np.load(f),95)) for f in glob.glob('$Q/null_q_f*.npy')),3))")
echo "sup critical value (largest 95th percentile over glint shares): $C"
$P -E -s -u - $Q <<'PYEOF2'
import glob, sys, numpy as np
from heidelberg.config import A
K = A.calibration
for f in sorted(glob.glob(f"{sys.argv[1]}/null_q_f*.npy")):
    q = np.load(f); b = np.random.default_rng(K.cal_q95_seed).integers(0, len(q), (K.cal_q95_boot, len(q)))
    print(f, "q95", round(float(np.percentile(q, 95)), 3), "bootstrap se", round(float(np.std(np.percentile(q[b], 95, axis=1))), 3))
PYEOF2
echo "null done $(date -u +%H:%M:%S)"
for KP in $(cfg "print(' '.join(f'{x:g}' for x in A.calibration.cal_power_kappas))"); do
  echo "power (kappa $KP)"; $P -E -s -u -m heidelberg.producers.likelihood_power $G $C power $KP
done
GX=$($P -E -s -c "import json; g = json.load(open('heidelberg/checks/P1_numbers.json'))['gradient']; print(f'{g[0]},{g[1]}')")
echo "the map's gradient (P1): $GX"
echo "claim rule, background map correct (size)"; $P -E -s -u -m heidelberg.producers.likelihood_power $G $C size correct
echo "claim rule, background uniform (the fitted map wrong)"; $P -E -s -u -m heidelberg.producers.likelihood_power $G $C size uniform
echo "claim rule, background map half the fitted gradient"; $P -E -s -u -m heidelberg.producers.likelihood_power $G $C size half
echo "power with refused pairs"; $P -E -s -u -m heidelberg.producers.likelihood_power $G $C refused
echo "power with the times in error"; $P -E -s -u -m heidelberg.producers.likelihood_power $G $C blur $Q/blur
echo "likelihood robustness"; $P -E -s -u -m heidelberg.producers.likelihood_robustness $G $C
echo "belt $(date -u +%H:%M:%S)"; $P -E -s -u -m heidelberg.producers.likelihood_limit $G $Q/blur
echo "done $(date -u +%H:%M:%S)"
} | tee $Q/calibration.txt
$P -E -s -u -m heidelberg.producers.calibration_numbers $Q/calibration.txt ${CAL_OUT:-heidelberg/checks/J3_numbers.json}   # the numbers for the draft (tools/numbers.py)
