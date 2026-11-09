"""The budget's numbers, from what the jobs measured: the dry run's Snakemake log (each rule's start and end), the scan
copies' timings and sizes (work/from_s3/scans/*.txt), the size in bytes of the dry run's products (PRODUCT_BYTES.txt,
from a listing of them on S3), the instance prices
(aws/prices.csv), and S3's prices (heidelberg/facts.json, AWS's price list for eu-north-1). Each item is the measured
time or size scaled to the record run's pairs, at on-demand prices (spot, used for measuring, costs less):
- measuring: the dry run's vCPU-seconds per pair (each measuring job's wall time times its three cores), for the record
  and control pairs on one r7i.48xlarge;
- part one's selection and part two's tests: the dry run's wall time from the last pair measured to the workspace's
  lock, and from the workspace's opening to the seal, scaled by the number of pairs, on one r7i.16xlarge;
- the calibration: its run's wall time, the sunlit study included (its log's last line), on one r7i.48xlarge;
- the scan copy: the mean time per scan at eight at a time, on one c7i.2xlarge;
- storage: from the copy to the report, the compressed scans in S3 Glacier Deep Archive, restored once in bulk for part
  one (measuring and the clock) and the restored copies kept for costs_restore_days at the Standard price; the products in
  S3 Standard-Infrequent Access.
python -m heidelberg.producers.costs DRY_SNAKEMAKE.log SCANS_DIR PRODUCT_BYTES.txt CALIBRATION.log OUT.json"""
import glob
import json
import re
import sys
from datetime import datetime

import pandas as pd

from heidelberg import config

LOG, SCANS, PROD, CAL, OUT = sys.argv[1:6]
F = config.facts()
PRICE = dict(pd.read_csv("aws/prices.csv").values)
MAIN = [pd.read_csv(f) for f in ("heidelberg/record_pairs.csv", "heidelberg/record_control_pairs.csv")]
N_PAIRS = sum(len(d) for d in MAIN)
N_SCANS = len(set().union(*(set(d.a) | set(d.b) for d in MAIN)))   # a different-night pair shares plates with the main sample
VCPU = 192                                                          # measuring on one r7i.48xlarge
N_DRY_CONTROLS = len(pd.read_csv("heidelberg/control_pairs.csv"))
N_DRY = len(pd.read_csv("heidelberg/measure_pairs_tuning.csv")) + N_DRY_CONTROLS
MONTHS = F["costs_storage_months"]


def jobs(log):
    """Each job's rule, start, and end, from a Snakemake log."""
    start, rule, out, t = {}, {}, {}, None
    lines = open(log).read().splitlines()
    for i, ln in enumerate(lines):
        m = re.match(r"^\[(\w{3} \w{3} +\d+ [\d:]+ \d{4})\]$", ln)
        if m:
            t = datetime.strptime(m.group(1), "%a %b %d %H:%M:%S %Y"); continue
        m = re.match(r"^(?:local)?rule (\w+):$", ln)
        if m:
            j = next((re.search(r"jobid: (\d+)", x).group(1) for x in lines[i:i + 8] if "jobid:" in x), None)
            if j:
                start[j], rule[j] = t, m.group(1)
        m = re.match(r"^Finished jobid: (\d+)", ln)
        if m:
            out[m.group(1)] = t
    return [(rule[j], start[j], out[j]) for j in out if j in start]


J = jobs(LOG)
meas = [(b - a).total_seconds() for r, a, b in J if r == "measure"]
vcpu_s_pair = 3 * sum(meas) / len(meas)
end_meas = max(b for r, a, b in J if r == "measure")
lock = max(b for r, a, b in J if r == "lock_workspace"); opened = min(a for r, a, b in J if r == "open_workspace")
seal = max(b for r, a, b in J if r == "seal")
part1_h = (lock - end_meas).total_seconds() / 3600 * N_PAIRS / N_DRY
part2_h = (seal - opened).total_seconds() / 3600 * N_PAIRS / N_DRY
measure_h = N_PAIRS * vcpu_s_pair / VCPU / 3600
cal_h = int(re.findall(r"^wall (\d+) s", open(CAL).read(), re.M)[-1]) / 3600
rows = [open(f).read().split() for f in glob.glob(f"{SCANS}/*.txt")]
dl_s = [float(r[r.index("download_s") + 1]) + float(r[r.index("pack_check_s") + 1]) for r in rows]
packed = [float(r[r.index("packed") + 1]) for r in rows]
copy_h = N_SCANS * (sum(dl_s) / len(dl_s)) / 8 / 3600
scans_gb = N_SCANS * (sum(packed) / len(packed)) / 1e9
prod_gb = N_PAIRS * int(open(PROD).read().split()[0]) / N_DRY / 1e9
dry_h = 2 * ((seal - min(a for r, a, b in J)).total_seconds() / 3600)                  # the two dry runs of both parts
cost = dict(measure=measure_h * PRICE["r7i.48xlarge"], part_one=part1_h * PRICE["r7i.16xlarge"],
            part_two=part2_h * PRICE["r7i.16xlarge"], calibration=cal_h * PRICE["r7i.48xlarge"],
            dry_runs=dry_h * PRICE["r7i.16xlarge"], scan_copy=copy_h * PRICE["c7i.2xlarge"],
            storage=scans_gb * (F["s3_price_gda_gb_month"] * MONTHS + F["s3_price_gda_bulk_gb"]
                                + F["s3_price_standard_gb_month"] * F["costs_restore_days"] / 30)
            + prod_gb * F["s3_price_sia_gb_month"] * MONTHS)
out = dict(pairs=N_PAIRS, scans=N_SCANS, measure_vcpu=VCPU, dry_pairs=N_DRY, dry_controls=N_DRY_CONTROLS, vcpu_s_per_pair=vcpu_s_pair, measure_hours=measure_h, part_one_hours=part1_h, part_two_hours=part2_h,
           calibration_hours=cal_h, scan_copy_hours=copy_h, scans_gb=scans_gb, products_gb=prod_gb, storage_months=MONTHS,
           cost=cost, total=sum(cost.values()))
aws = round(out["total"]); out.update(aws_rounded=aws, rerun_allowance=aws * (F["costs_rerun_factor"] - 1), request=aws * F["costs_rerun_factor"])
json.dump(out, open(OUT, "w"), indent=1)
for k, v in cost.items():
    print(f"{k}: ${v:.0f}")
print(f"request ${out['request']} (the AWS estimate ${aws} x {F['costs_rerun_factor']}); total ${out['total']:.0f} for {N_PAIRS} pairs (measuring {measure_h:.1f} h on {VCPU} vCPU; scans {scans_gb:.0f} GB)")
