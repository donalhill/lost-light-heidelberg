# The Heidelberg Lost Light study as one DAG (see README). Every rule runs in the project's image: on the laptop
# through ./run.sh, on AWS through aws/job.sh. /ref is the previous repository, mounted read-only, for the equality checks.
D = "work/dag"
PY = "python -u"
REF = "/ref"
PCAT = "work/gavo_svn_r9476/plattenkatalog.csv"
LSW_RAW = "work/gavo_snapshot_20261028/lsw_bruce_now.csv"
COV_RAW = "work/heidelberg_tuning/bruce_coverage.csv"
LSW, COV = f"{D}/lsw_timed.csv", f"{D}/coverage_timed.csv"
TUNING = "heidelberg/measure_pairs_tuning.csv"
SAME = f"{PY} tools/same.py"
# the tuning pairs' selection, offline from the 29 October dry run's frozen answers, as the previous repository made it
SEL_TUNING = ("--asteroids work/heidelberg_tuning/ast_h15 --scanner-marks work/dry_run/open/scanner_marks.csv "
              "--snapshots work/dry_catalogues --offline --blind")


rule check_sample:
    input: f"{D}/record_set.ok", f"{D}/controls.ok", f"{D}/geometry_set.ok"


rule check_selection:
    input: f"{D}/tuning_select.ok", f"{D}/derived.ok"


rule inputs:
    # every input against the manifest (heidelberg/inputs.sha256), once; a changed manifest checks again
    input: "heidelberg/inputs.sha256"
    output: f"{D}/inputs.ok"
    shell: "sha256sum -c --quiet {input} && touch {output}"


rule plate_times:
    input: PCAT, LSW_RAW, COV_RAW, "heidelberg/stages/plate_times.py"
    output: lsw=LSW, cov=COV
    shell: "{PY} -m heidelberg.stages.plate_times {PCAT} {LSW_RAW} {output.lsw} {COV_RAW} {output.cov} > {D}/plate_times.log"


rule geometry:
    input: LSW, "heidelberg/stages/record_geometry.py", "heidelberg/lib/shadow_geom.py", "heidelberg/lib/shadow_model.py"
    output: csv=f"{D}/geometry.csv", s=f"{D}/geometry_s.npy", ok=f"{D}/geometry.ok"
    shell: "{PY} -m heidelberg.stages.record_geometry {LSW} {D}/geometry > {D}/geometry.log && "
           "{SAME} {output.csv} work/heidelberg_record_geometry.csv && {SAME} {output.s} work/heidelberg_record_geometry_s.npy && touch {output.ok}"


rule record_set:
    input: LSW, rules.geometry.output.csv, TUNING, "heidelberg/stages/record_set.py", "heidelberg/lib/sample.py"
    output: rs=f"{D}/record_set.csv", ok=f"{D}/record_set.ok"
    shell: "{PY} -m heidelberg.stages.record_set {LSW} {D}/geometry.csv {TUNING} {output.rs} > {D}/record_set.log && "
           "{SAME} {output.rs} work/record_set.csv && {SAME} {output.rs} heidelberg/record_pairs.csv --drop-examined && touch {output.ok}"


rule controls:
    input: rules.record_set.output.rs, LSW, TUNING, "heidelberg/stages/control_pairs.py"
    output: csv=f"{D}/controls.csv", ok=f"{D}/controls.ok"
    shell: "{PY} -m heidelberg.stages.control_pairs {input[0]} {LSW} {TUNING} {output.csv} > {D}/controls.log && "
           "{SAME} {output.csv} heidelberg/record_control_pairs.csv && touch {output.ok}"


rule geometry_set:
    input: rules.record_set.output.rs, rules.geometry.output.csv, rules.geometry.output.s, "heidelberg/stages/geometry_set.py"
    output: s=f"{D}/geometry_set.npy", ok=f"{D}/geometry_set.ok"
    shell: "{PY} -m heidelberg.stages.geometry_set {input[0]} {input[1]} {input[2]} {output.s} && "
           "{SAME} {output.s} work/heidelberg_record_geometry_set.npy && touch {output.ok}"


rule tuning_fidelity:
    # each tuning pair's injection-fidelity report recomputed (PAIR_fidelity.json, written beside the products, where the
    # selection reads it); the verdicts must be the 29 October instance's
    input: f"{D}/inputs.ok", "heidelberg/stages/fidelity.py", "heidelberg/lib/fidelity.py"
    output: f"{D}/tuning_fidelity.ok"
    shell: "{PY} -m heidelberg.stages.fidelity work/dry_prod > {D}/tuning_fidelity.log && "
           "diff <(awk '{{print $1, $2}}' {D}/tuning_fidelity.log) <(head -n $(wc -l < {D}/tuning_fidelity.log) work/dry_run/open/fidelity.log | awk '{{print $1, $2}}') && touch {output}"


rule tuning_select:
    # the tuning pairs' selection, offline from the 29 October dry run's frozen answers
    input: f"{D}/tuning_fidelity.ok", "heidelberg/stages/select.py", "heidelberg/lib/selection.py", "heidelberg/lib/catalogues.py",
           "heidelberg/lib/twin_ratio.py", "heidelberg/lib/decision.py"
    output: directory(f"{D}/tuning_select")
    threads: 6
    shell: "rm -rf {D}/xc_tuning && cp -r work/dry_unsealed/xcache {D}/xc_tuning && "
           "{PY} -m heidelberg.stages.select work/dry_prod {output} {SEL_TUNING} --jobs {threads} --xcache {D}/xc_tuning > {D}/tuning_select.log"


rule tuning_select_check:
    # check 2: equal to the previous repository's (its work/dag/select_primary, made after the last change to the
    # selection code) and to the 29 October instance's
    input: rules.tuning_select.output
    output: f"{D}/tuning_select.ok"
    shell: "{PY} tools/same_dir.py {input} {REF}/work/dag/select_primary > {D}/tuning_select_ref.txt || true; "
           "{PY} -m heidelberg.tests.instance_check {input} work/dry_unsealed/select/sealed work/dry_run/open/select > {D}/tuning_select_instance.txt || true; "
           "tail -1 {D}/tuning_select_ref.txt {D}/tuning_select_instance.txt; "
           "grep -q 'all equal' {D}/tuning_select_ref.txt && grep -q IDENTICAL {D}/tuning_select_instance.txt && touch {output}"


rule derived:
    # the values measured on the tuning pairs (tuning background, overdispersion, median half-efficiency), recomputed and
    # compared with heidelberg/derived.json
    input: rules.tuning_select.output, "heidelberg/stages/derived.py"
    output: f"{D}/derived.ok"
    shell: "{PY} -m heidelberg.stages.derived {input[0]} {D}/derived.json > {D}/derived.log && cmp {D}/derived.json heidelberg/derived.json && touch {output}"


# check 4: the shadow stage on the tuning pairs, the previous repository's code (exported at b34e6f0 into
# work/ref_b34e6f0) against this one, on the same selection with the same seeds and reduced toys; they must agree except
# for the cut items (the two-sided interval, removed from the candidate outputs)
SHADOW_RUNS = {"chance": "chance area topo", "injection": "injection area topo", "deficit": "deficit area topo",
               "candidate_topo": "candidate area topo", "candidate_geo": "candidate area geo"}
QUICK = "candidate_toys=40 null_toys=40 belt_toys=40 check_draws=400"


rule check_shadow:
    input: f"{D}/tuning_shadow.ok"


rule tuning_shadow:
    input: sel=rules.tuning_select.output, lsw=LSW, code=["heidelberg/stages/shadow.py", "heidelberg/lib/likelihood.py"]
    output: directory(f"{D}/tuning_shadow/{{tree}}/{{run}}")
    threads: 6
    params: root=lambda w: "work/ref_b34e6f0" if w.tree == "ref" else ".", args=lambda w: SHADOW_RUNS[w.run]
    shell: "mkdir -p {output} && (cd {params.root} && {PY} -m heidelberg.config override $OLDPWD/{output}/config.json {QUICK}) && "
           "HD_CONFIG=$PWD/{output}/config.json HD_JOBS={threads} PYTHONPATH={params.root} "
           "{PY} -m heidelberg.stages.shadow work/dry_prod {input.sel} {input.lsw} {output} {params.args} > {output}/run.log 2>&1"


rule tuning_shadow_compare:
    input: expand(f"{D}/tuning_shadow/{{tree}}/{{run}}", tree=["ref", "new"], run=SHADOW_RUNS)
    output: f"{D}/tuning_shadow.ok"
    shell: "{PY} tools/same_shadow.py {D}/tuning_shadow/new {D}/tuning_shadow/ref > {D}/tuning_shadow.txt && touch {output}"


# ---- measuring (AWS) ----------------------------------------------------------------------------------------------
# A pair's products go to work/ll/products/MKEY/PAIR/, which aws/job.sh mirrors to s3://BUCKET/ll/. MKEY hashes what the
# measurement reads, listed here and nothing else: its code, the settings sections it reads, and the package lock.
# Unchanged MKEY: the products on S3 are reused; changed: every pair is measured again.
import glob, hashlib, os, json as _json
MEASURE_CODE = ["heidelberg/stages/measure_pair.py", "heidelberg/lib/measure.py", "heidelberg/lib/meridian.py",
                "heidelberg/lib/catalogues.py", "heidelberg/config.py", "requirements.lock", "Dockerfile"]
MEASURE_SECTIONS = ["background", "detection", "injections", "pairing", "physics", "plates", "sample", "selection", "solve"]


def key(files, sections):
    h = hashlib.sha256()
    for f in files:
        h.update(f.encode()); h.update(open(f, "rb").read())
    cfg = _json.load(open("config/analysis.json"))
    h.update(_json.dumps({s: cfg[s] for s in sections}, sort_keys=True).encode())
    return h.hexdigest()[:12]


MKEY = key(MEASURE_CODE, MEASURE_SECTIONS)
PROD = f"work/ll/products/{MKEY}"
# the S3 prefix of the compressed, checked scan copies: the record run's, else the final dry run's (the tuning pairs and
# their different-night controls)
SCANS = config.get("scans", "heidelberg_scans" if config.get("run") == "record" else "heidelberg_scans_dry")
# a pair is named by both its plates, A+B: a different-night control shares its first plate with a pair of the set
pair_ids = lambda f: ["+".join(l.split(",")[:2]) for l in open(f).read().splitlines()[1:]]
TUNING_PAIRS = pair_ids(TUNING)
CONTROL_IDS = set(pair_ids("heidelberg/control_pairs.csv")) | set(pair_ids("heidelberg/record_control_pairs.csv"))
GQ = lambda w: "heidelberg/gq_control" if f"{w.a}+{w.b}" in CONTROL_IDS else "heidelberg/gq_tuning"   # the saved Gaia query (a control's box covers both plates)


rule measure:
    # one pair from its compressed copies, restored with funpack; the saved Gaia query reused
    output: directory(f"{PROD}/{{a}}+{{b}}")
    wildcard_constraints: a="B[0-9]+[ab]", b="B[0-9]+[ab]"
    threads: 3
    resources: mem_gb=24
    params: b=lambda w: w.b, gq=GQ, scans=lambda w: SCANS
    shell: "set -e; if aws s3 ls s3://$BUCKET/ll/products/{MKEY}/{wildcards.a}+{wildcards.b}/{wildcards.a}.json > /dev/null; then "
           "aws s3 sync --quiet s3://$BUCKET/ll/products/{MKEY}/{wildcards.a}+{wildcards.b} {output}; exit 0; fi; "      # measured before: reused
           "t=/tmp/m_{wildcards.a}+{wildcards.b}; rm -rf $t; mkdir -p $t {output}; "      # by both plates: a control shares its first plate
           "for p in {wildcards.a} {params.b}; do aws s3 cp --quiet s3://$BUCKET/{params.scans}/$p.fits.fz $t/ && funpack -O $t/$p.fits $t/$p.fits.fz && rm $t/$p.fits.fz || exit 1; done; "
           "cp {params.gq}/{wildcards.a}_gq_*.npz {output}/ 2>/dev/null || true; "
           "HD_THREADS={threads} {PY} -m heidelberg.stages.measure_pair $t/{wildcards.a}.fits $t/{params.b}.fits {output} > {output}/{wildcards.a}.log 2>&1; "
           "rm -rf $t"


rule restore_check:
    # a tuning scan's copy (27 October) restores the same data and cards as the final dry run's verified copy
    output: "work/ll/restore/{p}.txt"
    shell: "set -e; t=/tmp/r_{wildcards.p}; mkdir -p $t; aws s3 cp --quiet s3://$BUCKET/heidelberg_scans_dry/{wildcards.p}.fits.fz $t/d.fz; "
           "funpack -O $t/d.fits $t/d.fz; aws s3 cp --quiet s3://$BUCKET/heidelberg_tuning/fits/{wildcards.p}.fits.fz $t/t.fz; "
           "funpack -S $t/t.fz | {PY} tools/scan_copy.py $t/d.fits && echo 'restores the same data and cards' > {output}; rm -rf $t"


rule check_products:
    # check 3: every tuning scan's copy restores; the tuning pairs' products equal heidelberg_dry3's within rounding; and
    # the selection made from each (offline, the frozen answers) gives the same verdicts
    input: restore=expand("work/ll/restore/{p}.txt", p=[x for ab in TUNING_PAIRS for x in ab.split("+")]),
           prod=expand(f"{PROD}/{{ab}}", ab=TUNING_PAIRS), sel=["work/ll/check3/select_new", "work/ll/check3/select_dry3"]
    output: "work/ll/check3/check3.ok"
    params: rep="work/ll/check3/check_products.txt"          # the report, kept whatever the outcome
    shell: "{PY} tools/same_products.py {PROD} work/ll/check3/prod_dry3 {input.restore} > {params.rep} 2>&1 || true; "
           "{PY} tools/same_dir.py work/ll/check3/select_new work/ll/check3/select_dry3 >> {params.rep} 2>&1 || true; cat {params.rep}; "
           "grep -q 'equal within rounding' {params.rep} && grep -q 'all equal' {params.rep} && touch {output}"


rule check3_products:
    # the products to select from: the new ones gathered into one folder (links), and heidelberg_dry3's
    input: expand(f"{PROD}/{{ab}}", ab=TUNING_PAIRS)
    output: new=directory("work/ll/check3/prod_new"), dry3=directory("work/ll/check3/prod_dry3")
    shell: "mkdir -p {output.new} && for d in {input}; do cp -l $d/* {output.new}/; done && "
           "aws s3 sync --quiet s3://$BUCKET/heidelberg_dry3/prod {output.dry3} --exclude '*_fidelity.json'"


rule check3_select:
    input: "work/ll/check3/prod_{which}"
    output: directory("work/ll/check3/select_{which}")
    threads: 16
    shell: "{PY} -m heidelberg.stages.fidelity {input} > {input}_fidelity.log && "
           "rm -rf {output}_xc && cp -r work/dry_unsealed/xcache {output}_xc && "
           "{PY} -m heidelberg.stages.select {input} {output} {SEL_TUNING} --jobs {threads} --xcache {output}_xc > {output}.log"


# ---- external data, fetched once (AWS) -----------------------------------------------------------------------------
URLS = {}
for f in ["heidelberg/record_pairs.csv", "heidelberg/record_control_pairs.csv", "heidelberg/measure_pairs_tuning.csv"]:
    h, *rows = open(f).read().splitlines(); h = h.split(",")
    for r in rows:
        r = dict(zip(h, r.split(","))); URLS.setdefault(r["a"], r["ua"]); URLS.setdefault(r["b"], r["ub"])


# every scan the main sample and its different-night pairs need, copied once (aws/launch.sh ... c7i.2xlarge ... scan_copies)
SCAN_PLATES = sorted({p for f in ("heidelberg/record_pairs.csv", "heidelberg/record_control_pairs.csv")
                      for r in open(f).read().splitlines()[1:] for p in r.split(",")[:2]})


rule scan_copies:
    input: [f"work/ll/scans/{p}.txt" for p in SCAN_PLATES]


rule copy_scan:
    # one plate's scan from GAVO, compressed (fpack) and kept at s3://BUCKET/heidelberg_scans/ (Glacier Deep Archive; restored
    # in bulk before part one by tools/restore_scans.sh) only if it restores the
    # scan's pixels byte for byte and its cards (tools/scan_copy.py); the original's SHA-256, the timings, and any
    # non-zero padding in the original (restored as zeros) recorded
    output: "work/ll/scans/{p}.txt"
    resources: gavo=1
    params: url=lambda w: URLS[w.p]
    shell: "set -e; t=/tmp/s_{wildcards.p}; mkdir -p $t; t0=$(date +%s); "
           "curl -sS --fail --retry 5 --retry-delay 20 --max-time 3600 -o $t/{wildcards.p}.fits '{params.url}'; t1=$(date +%s); "
           "fpack -Y $t/{wildcards.p}.fits; note=$(funpack -S $t/{wildcards.p}.fits.fz | {PY} tools/scan_copy.py $t/{wildcards.p}.fits); t2=$(date +%s); "
           "aws s3 cp --quiet --storage-class DEEP_ARCHIVE $t/{wildcards.p}.fits.fz s3://$BUCKET/heidelberg_scans/{wildcards.p}.fits.fz; "
           "echo \"{wildcards.p} sha256 $(sha256sum $t/{wildcards.p}.fits | cut -c1-64) bytes $(stat -c %s $t/{wildcards.p}.fits) "
           "packed $(stat -c %s $t/{wildcards.p}.fits.fz) download_s $((t1 - t0)) pack_check_s $((t2 - t1)) $note restores\" | tr -s ' ' > {output}; rm -rf $t"


rule fetch_vsx:
    # the whole VSX catalogue (VizieR B/vsx/vsx: name and J2000 position), dated, by strips of declination
    # (tools/fetch_vsx.sh: each strip from the first VizieR mirror that answers, with retries, halved when no mirror
    # serves it whole; strips already fetched are kept, so a rerun resumes); read by catalogues.vsx_local in place of
    # the per-field cone queries to AAVSO
    output: "work/ll/vsx/B_vsx_vsx.tsv"
    shell: "bash tools/fetch_vsx.sh $(dirname {output}) && {PY} tools/join_strips.py $(dirname {output})/strips {output} "
           "&& date -u +%Y-%m-%dT%H:%M:%SZ > $(dirname {output})/fetched_utc.txt"


rule profile_measure:
    # one pair measured under cProfile (not kept as products): where the measurement's time goes
    output: "work/ll/profile/{a}+{b}.txt"
    threads: 3
    resources: mem_gb=24
    params: b=lambda w: w.b
    shell: "set -e; t=/tmp/p_{wildcards.a}+{wildcards.b}; mkdir -p $t/out; "
           "for p in {wildcards.a} {params.b}; do aws s3 cp --quiet s3://$BUCKET/{SCANS}/$p.fits.fz $t/ && funpack -O $t/$p.fits $t/$p.fits.fz && rm $t/$p.fits.fz; done; "
           "cp heidelberg/gq_tuning/{wildcards.a}_gq_*.npz $t/out/ 2>/dev/null || true; s=$(date +%s); "
           "HD_THREADS={threads} {PY} -m cProfile -o $t/m.prof -m heidelberg.stages.measure_pair $t/{wildcards.a}.fits $t/{params.b}.fits $t/out > $t/run.log 2>&1; "
           "echo \"wall $(( $(date +%s) - s )) s on $(nproc --all) vCPU\" > {output}; "
           "{PY} -c \"import pstats; s = pstats.Stats('$t/m.prof'); s.sort_stats('tottime').print_stats(30); s.sort_stats('cumulative').print_stats(40)\" >> {output}; rm -rf $t"


# ---- calibration (AWS) ---------------------------------------------------------------------------------------------
# The toys on the record geometry (heidelberg/producers/calibration.sh: the null and its critical value, power, size,
# robustness, and the belt's expected limit and coverage, also with scaled chance estimates and the clock's timing
# errors), then the sunlit-pairs study (J6). Keyed like the measurement. mode quick: a few toys each, to test the
# plumbing on the laptop; full: the configured counts.
CAL_CODE = ["heidelberg/lib/likelihood.py", "heidelberg/stages/shadow.py", "heidelberg/producers/calibration.sh",
            "heidelberg/producers/likelihood_null.py", "heidelberg/producers/likelihood_power.py",
            "heidelberg/producers/likelihood_limit.py", "heidelberg/producers/likelihood_robustness.py",
            "heidelberg/producers/blurred_geometry.py", "heidelberg/producers/calibration_numbers.py",
            "heidelberg/producers/sunlit_size.py", "heidelberg/derived.json", "heidelberg/checks/P1_numbers.json",
            "heidelberg/checks/N5_clock.json", "work/record_set.csv", "work/heidelberg_record_geometry_set.npy",
            "heidelberg/config.py", "requirements.lock", "Dockerfile"]
CAL_SECTIONS = ["statistics", "calibration", "background", "sample"]
CKEY = key(CAL_CODE, CAL_SECTIONS) if os.path.exists("work/heidelberg_record_geometry_set.npy") else "missing"
QUICK_CAL = ("cal_null_toys=40 power_toys=20 cal_belt_toys=20 cal_check_draws=100 cal_blur_sets=2 "
             "cal_sunlit_null_toys=40 cal_sunlit_toys=20 cal_q95_boot=50 cal_pooled_toys=20")


rule calibrate:
    output: directory(f"work/ll/calibration/{CKEY}/{{mode}}")
    params: cfg=lambda w: f"{PY} -m heidelberg.config override $PWD/work/ll/calibration/{CKEY}/{w.mode}/config.json {QUICK_CAL} && export HD_CONFIG=$PWD/work/ll/calibration/{CKEY}/{w.mode}/config.json; " if w.mode == "quick" else ""
    log: f"work/ll/calibration/{CKEY}/{{mode}}.log"        # outside the output folder, which Snakemake removes on failure
    shell: "set -e; mkdir -p {output}; {params.cfg} s=$(date +%s); export HD_JOBS=$(nproc --all); "      # every producer, the sunlit study too
           "PY=python CAL_DIR={output} CAL_OUT={output}/J3_numbers.json bash heidelberg/producers/calibration.sh > {log} 2>&1; "
           "{PY} -m heidelberg.producers.sunlit_size work/heidelberg_record_geometry_set.npy {output}/J3_numbers.json {output}/J6_sunlit.json >> {log} 2>&1; "
           "echo \"wall $(( $(date +%s) - s )) s on $(nproc --all) vCPU\" | tee {output}/timing.txt >> {log}"


rule calibrate_quick:
    input: f"work/ll/calibration/{CKEY}/quick"


rule calibrate_full:
    input: f"work/ll/calibration/{CKEY}/full"


# ---- the run: part one and part two (AWS) --------------------------------------------------------------------------
# One set of rules for the dry run (the tuning pairs, config run=dry) and the record run (run=record). Open outputs go
# under work/ll/runs/RUN/ (mirrored to S3 as the job runs); everything that carries a candidate's result stays in
# sealed/RUN/ on the machine's own disk, which nothing uploads, and leaves it only time-locked (tools/seal.py).
RUN = config.get("run", "dry")
RUNS = {"dry": dict(pairs=TUNING, controls="heidelberg/control_pairs.csv", record_set="work/dry_record_set.csv",
                    snapshots="work/dry_catalogues", lock_part1="+300", lock_part2="+300"),
        "record": dict(pairs="heidelberg/record_pairs.csv", controls="heidelberg/record_control_pairs.csv",
                       record_set="work/record_set.csv", snapshots="work/ll/snapshots/record",
                       lock_part1=config.get("unlock_select_utc", "unset"), lock_part2="2027-04-01T00:00:00Z")}
R_ = RUNS[RUN]; P1, P2, SW = f"work/ll/runs/{RUN}/part1", f"work/ll/runs/{RUN}/part2", f"sealed/{RUN}"
REG = f"work/ll/runs/{RUN}/registered"
LISTS = {"rec": pair_ids(R_["pairs"]), "ctrl": pair_ids(R_["controls"])}
SEL_RUN = (f"--asteroids work/heidelberg_tuning/ast_h15 --scanner-marks {P1}/open/scanner_marks.csv "
           f"--snapshots {R_['snapshots']}")


rule snapshots_record:
    # the record run's catalogue snapshots: the pinned ones, with the whole VSX catalogue in place of the per-field answers
    input: vsx="work/ll/vsx/B_vsx_vsx.tsv"
    output: directory("work/ll/snapshots/record")
    shell: "mkdir -p {output} && cp work/heidelberg_catalogues/*.tsv work/heidelberg_catalogues/*.ecsv {output}/ && cp {input.vsx} {output}/"


rule flat_products:
    # one folder of a list's products (links), each pair's fidelity report added beside them
    input: lambda w: expand(f"{PROD}/{{ab}}", ab=LISTS[w.which])
    output: directory(f"{P1}/prod_{{which}}")
    shell: "mkdir -p {output} && for d in {input}; do cp -l $d/* {output}/; done && {PY} -m heidelberg.stages.fidelity {output} > {output}_fidelity.log"


rule scanner_marks:
    input: f"{P1}/prod_rec"
    output: f"{P1}/open/scanner_marks.csv"
    shell: "{PY} -m heidelberg.stages.backgrounds {input} {P1}/open > {P1}/open/backgrounds.log"


rule coverage:
    # the coverage map, footprints, and sample manifest (registered before any analysis)
    input: prod=f"{P1}/prod_{{which}}", lsw=f"{REG}/rec/plate_times.csv"
    output: [f"{REG}/{{which}}/{n}" for n in ("coverage_map.csv", "coverage_footprints.json", "sample_manifest.csv")]
    params: role=lambda w: "different-night" if w.which == "ctrl" else ""
    shell: "{PY} -m heidelberg.stages.coverage_map {input.prod} {input.lsw} {R_[record_set]} {REG}/{wildcards.which} {params.role} > {REG}/{wildcards.which}_coverage.log"


rule clock:
    # the asteroid clock on the main sample's plates: a failure is reported and changes nothing (rule clock_times)
    input: f"{P1}/prod_rec"
    output: f"{REG}/rec/clock_summary.json"
    shell: "c={P1}/open/clock; mkdir -p $c/parts; "
           "( {PY} -m heidelberg.stages.clock_plates {input} $c/plates.csv && head -1 $c/plates.csv > $c/h.csv "
           "&& tail -n +2 $c/plates.csv | split -n l/$(nproc --all) -d -a 3 - $c/parts/p && for f in $c/parts/p???; do cat $c/h.csv $f > $f.csv; "
           "{PY} -m heidelberg.stages.clock_batch work/heidelberg_tuning/ast_h15 $f.csv $BUCKET/{SCANS} {input} $f.out.csv > $f.log 2>&1 & done; wait; "
           "{PY} -m heidelberg.producers.clock_summary {output} $c/parts/p???.out.csv > $c/summary.txt 2>&1 ) "
           "|| echo '{{\"repeat_tests\": false, \"failed\": true, \"by_plate\": {{}}}}' > {output}"


rule clock_times:
    # the plate times part one registers and part two uses: the catalogue's, shifted by the clock's median offset if clear
    input: sm=f"{REG}/rec/clock_summary.json", lsw=LSW, code="heidelberg/stages/clock_correct.py"
    output: lsw=f"{REG}/rec/plate_times.csv", js=f"{REG}/rec/clock_correction.json"
    shell: "{PY} -m heidelberg.stages.clock_correct {input.sm} {input.lsw} {output.lsw} {output.js} > {REG}/rec_clock_times.log"


rule select_blind:
    # the blind selection, online (catalogue answers saved in its cache), into the sealed workspace
    input: prod=f"{P1}/prod_rec", marks=f"{P1}/open/scanner_marks.csv", snap=R_["snapshots"]
    output: sel=directory(f"{SW}/select"), xc=directory(f"{SW}/xcache")
    threads: workflow.cores
    shell: "{PY} tools/check_ephemerides.py work/heidelberg_tuning/ast_h15/plates.csv {R_[pairs]} {R_[controls]} > {P1}/open/ephemerides_check.txt && "
           "cp -r work/heidelberg_xmatch_cache {output.xc} && "
           "{PY} -m heidelberg.stages.select {input.prod} {output.sel} {SEL_RUN} --blind --jobs {threads} --xcache {output.xc} > {P1}/open/select.log 2>&1"


rule select_repeat:
    # the selection again, offline from the frozen answers: every file must be the same
    input: prod=f"{P1}/prod_rec", sel=f"{SW}/select", xc=f"{SW}/xcache"
    output: f"{P1}/open/reproducibility.txt"
    threads: workflow.cores
    shell: "rm -rf {SW}/select_offline && "
           "{PY} -m heidelberg.stages.select {input.prod} {SW}/select_offline {SEL_RUN} --blind --offline --jobs {threads} --xcache {input.xc} > {P1}/open/select_offline.log 2>&1 && "
           "{PY} tools/same_selection.py {input.sel} {SW}/select_offline > {output}.tmp; rm -rf {SW}/select_offline; "
           "tail -1 {output}.tmp; grep -q 'reproducible' {output}.tmp && mv {output}.tmp {output}"


rule lock_workspace:
    # the transient sample's checksum registered, and the workspace (verdicts and frozen answers) locked until part two
    input: sel=f"{SW}/select", xc=f"{SW}/xcache", rep=f"{P1}/open/reproducibility.txt"
    output: sha=f"{REG}/rec/sample_sha256.txt", tle=f"{P1}/workspace.tgz.tle", tsha=f"{REG}/rec/workspace.tgz.tle.sha256"
    shell: "{PY} -m heidelberg.stages.sample_checksum {input.sel} > {output.sha} && tar -C {SW} -czf {SW}/workspace.tgz select xcache && "
           "{PY} tools/seal.py lock {SW}/workspace.tgz {R_[lock_part1]} {output.tle} && cp {output.tle}.sha256 {output.tsha} && "
           "cp {output.tle}.round {REG}/rec/workspace_round.txt && rm {SW}/workspace.tgz"


rule part_one:
    input: [f"{REG}/{w}/{n}" for w in ("rec", "ctrl") for n in ("coverage_map.csv", "coverage_footprints.json", "sample_manifest.csv")],
           f"{REG}/rec/clock_summary.json", f"{REG}/rec/clock_correction.json", f"{REG}/rec/sample_sha256.txt", f"{P1}/workspace.tgz.tle"
    output: ok=f"{P1}/done.ok", reg=f"{REG}/registered.sha256"
    shell: "(cd {REG} && find . -type f ! -name registered.sha256 | sort | xargs sha256sum) > {output.reg} && touch {output.ok}"


# part two, on a new machine, after the registration is posted and the workspace's round is out
rule open_workspace:
    input: tle=f"{P1}/workspace.tgz.tle", tsha=f"{REG}/rec/workspace.tgz.tle.sha256", sha=f"{REG}/rec/sample_sha256.txt",
           man=f"{REG}/rec/sample_manifest.csv", prod=f"{P1}/prod_rec", man_c=f"{REG}/ctrl/sample_manifest.csv", prod_c=f"{P1}/prod_ctrl",
           reg=f"{REG}/registered.sha256"
    output: directory(f"{SW}/opened")
    shell: "set -e; mkdir -p {output} {P2}/open; (cd {REG} && sha256sum --check --quiet registered.sha256) > {P2}/open/registered_check.txt; "   # every registered file
           "{PY} tools/check_manifest.py {input.man} {input.prod} > {P2}/open/sample_check.txt; "
           "{PY} tools/check_manifest.py {input.man_c} {input.prod_c} > {P2}/open/control_check.txt; "
           + ("{PY} tools/registration.py check {REG}/rec > {P2}/open/registration_check.txt; " if RUN == "record" else "") +
           "{PY} tools/seal.py open {input.tle} $(cat {input.tsha}) {SW}/workspace.tgz && tar -C {output} -xzf {SW}/workspace.tgz && rm {SW}/workspace.tgz && "
           "{PY} -m heidelberg.stages.sample_checksum {output}/select \"$(cat {input.sha})\" > {P2}/open/sample_check_registered.txt"


SHADOW_POPS = {"chance": "chance area topo", "injection": "injection area topo", "deficit": "deficit area topo",
               "candidate_topo": "candidate area topo", "candidate_geo": "candidate area geo"}


rule shadow:
    input: op=f"{SW}/opened", prod=f"{P1}/prod_rec", lsw=f"{REG}/rec/plate_times.csv"
    output: directory(f"{SW}/shadow/{{run}}")
    threads: workflow.cores
    params: args=lambda w: SHADOW_POPS[w.run]
    shell: "HD_JOBS={threads} {PY} -m heidelberg.stages.shadow {input.prod} {input.op}/select {input.lsw} {output} {params.args} > {output}.log 2>&1"


rule part_two_open:
    # what part two reports before the unlock: the null tests, the positive control, the comparison of data quality, the
    # summary, and the different-night controls' static bound (nothing of the candidates)
    input: sh=expand(f"{SW}/shadow/{{run}}", run=["chance", "injection", "deficit"]), op=f"{SW}/opened",
           prod=f"{P1}/prod_rec", ctrl=f"{P1}/prod_ctrl", lsw=f"{REG}/rec/plate_times.csv", marks=f"{P1}/open/scanner_marks.csv"
    output: f"{P2}/open/role_report.csv", f"{P2}/open/static_bound.json"
    threads: workflow.cores
    shell: "set -e; o={P2}/open; for d in {input.sh}; do cp $d/null_*.json $o/; done; "
           "{PY} -m heidelberg.stages.summarise {input.prod} {input.op}/select > $o/summary.txt 2>&1; "
           "{PY} -m heidelberg.stages.role_report {R_[record_set]} {input.prod} {input.op}/select - {input.lsw} $o/role_report.csv > $o/role_report.txt 2>&1; "
           "cp -r work/heidelberg_xmatch_cache {SW}/xc_ctrl && {PY} -m heidelberg.stages.select {input.ctrl} $o/select_control {SEL_RUN} --jobs {threads} --xcache {SW}/xc_ctrl > $o/select_control.log 2>&1; "
           "{PY} -m heidelberg.stages.static_bound $o/select_control $o/static_bound.json > $o/static_bound.txt 2>&1"


rule seal:
    # the candidates' results locked to the unlock's round; the encrypted file's checksum is then posted (seal_record)
    input: sh=expand(f"{SW}/shadow/{{run}}", run=["candidate_topo", "candidate_geo"]), op=f"{SW}/opened",
           open_=f"{P2}/open/role_report.csv"
    output: tle=f"{P2}/sealed.tgz.tle"
    shell: "tar -C {SW} -czf {SW}/sealed.tgz opened/select/sealed shadow/candidate_topo shadow/candidate_geo opened/xcache && "
           "{PY} tools/seal.py lock {SW}/sealed.tgz {R_[lock_part2]} {output.tle} && rm {SW}/sealed.tgz"


rule part_two:
    input: f"{P2}/sealed.tgz.tle", f"{P2}/open/static_bound.json"
    output: f"{P2}/done.ok"
    shell: "touch {output}"


# after the unlock (1 April 2027 for the record run): the seal record, the unsealed results, the verdict and the rate by
# frozen code, and the Zenodo deposit. The record run's seal record is written by hand once the encrypted file's checksum
# is posted (tools/registration.py seal REG TLE URL POSTED_UTC); a dry run writes its own.
U = f"work/ll/runs/{RUN}/after_unlock"
if RUN != "record":
    rule seal_record_dry:
        input: f"{P2}/sealed.tgz.tle"
        output: f"{REG}/seal_record.json"
        shell: "{PY} tools/registration.py seal {REG} {input} dry-run $(date -u +%Y-%m-%dT%H:%M:%SZ)"

rule unlock:
    input: tle=f"{P2}/sealed.tgz.tle", rec=f"{REG}/seal_record.json"
    output: directory(f"{U}/unsealed")
    params: rnd="-" if RUN != "record" else str(_json.load(open("config/analysis.json"))["blind"]["drand_round"]["value"])   # a dry run's own round stands
    shell: "{PY} tools/unseal.py {input.tle} {input.rec} tle $({PY} -c 'from heidelberg import config; print(config.facts()[\"drand_chain\"])') "
           "{params.rnd} {output}"

rule verdict:
    input: uns=f"{U}/unsealed", role=f"{P2}/open/role_report.csv", cov=f"{REG}/rec/coverage_map.csv"
    output: v=f"{U}/results/verdict.json", r=f"{U}/results/rate.json"
    shell: "set -e; d={U}/verdict_in; rm -rf $d; mkdir -p $d; cp {input.uns}/shadow/*/sealed_*.json {P2}/open/null_*.json $d/; "
           "{PY} -m heidelberg.stages.verdict $d {input.role} {output.v} > {U}/results/verdict.txt; "
           "{PY} -m heidelberg.stages.rate $d/sealed_candidate_area_topo.json {input.cov} {output.r} > {U}/results/rate.txt"

rule release:
    input: uns=f"{U}/unsealed", v=f"{U}/results/verdict.json"
    output: directory(f"{U}/release")
    shell: "{PY} tools/release.py $(git rev-parse HEAD) {REG} {P2}/open {input.uns} {U}/results {output} > {U}/release.log"


rule after_unlock:
    input: f"{U}/release"


# ---- results and documents (laptop) --------------------------------------------------------------------------------
rule results:
    input: "tools/results.py", "config/analysis.json", "heidelberg/facts.json", "heidelberg/derived.json",
           [f"heidelberg/checks/{f}.json" for f in ("S5_sample", "P1_numbers", "N5_clock", "J3_numbers", "J6_sunlit", "costs", "run_numbers", "forecast")]
    output: "work/results.json"
    shell: "{PY} tools/results.py {output}"


rule draft:
    input: "draft/preregistration_v2.md"


rule render:
    wildcard_constraints: doc="[a-z0-9_]+"
    input: tpl="draft/{doc}.src.md", res="work/results.json", code="tools/render.py"
    output: "draft/{doc}.md"
    shell: "{PY} tools/render.py {input.tpl} {input.res} {output}"


# ---- the numbers the documents state, from this repository's own outputs --------------------------------------------
# sample: the rebuilt plate table, record set, and geometry, with the dry run's coverage map of the tuning pairs; pilot:
# the tuning selection, its summary, the dry run's different-night controls, and the tuning background map (check 4's
# chance run); clock: the dry run's asteroid clock on the tuning plates; calibration and sunlit: the calibration's run.
# What the AWS jobs made is read from work/from_s3/, a copy of s3://BUCKET/ll/ that no rule writes, taken by hand:
#   aws s3 sync s3://BUCKET/ll/runs/dry work/from_s3/runs/dry; aws s3 sync s3://BUCKET/ll/calibration work/from_s3/calibration
DRY = "work/from_s3/runs/dry"


rule sample_numbers:
    input: lsw=LSW, rs=f"{D}/record_set.csv", geo=f"{D}/geometry.csv", cov=f"{DRY}/registered/rec/coverage_map.csv",
           pairs=f"{D}/geometry_pairs.csv", code="heidelberg/producers/sample_numbers.py"
    output: "heidelberg/checks/S5_sample.json"
    shell: "{PY} -m heidelberg.producers.sample_numbers {input.lsw} {input.rs} {input.geo} {input.cov} {input.pairs} {TUNING} {LSW_RAW} {output} > {D}/sample_numbers.log"


rule tuning_summary:
    input: f"{D}/tuning_select", "heidelberg/stages/summarise.py"
    output: f"{D}/summary_tuning.txt"
    shell: "{PY} -m heidelberg.stages.summarise work/dry_prod {input[0]} > {output} 2> /dev/null"


rule pilot_numbers:
    input: sel=f"{D}/tuning_select", summ=f"{D}/summary_tuning.txt", ctrl=f"{DRY}/part2/open/select_control",
           map=f"{D}/tuning_shadow/new/chance", code="heidelberg/producers/pilot_numbers.py"
    output: "heidelberg/checks/P1_numbers.json"
    shell: "{PY} -m heidelberg.producers.pilot_numbers {input.sel} {input.summ} {input.ctrl} {input.map}/null_chance_area_topo.json {output} > {D}/pilot_numbers.log"


rule clock_numbers:
    input: f"{DRY}/registered/rec/clock_summary.json"
    output: "heidelberg/checks/N5_clock.json"
    shell: "cp {input} {output}"


rule calibration_numbers:
    input: f"work/from_s3/calibration/{CKEY}/full/J3_numbers.json"
    output: "heidelberg/checks/J3_numbers.json"
    shell: "cp {input} {output}"


rule sunlit_numbers:
    input: f"work/from_s3/calibration/{CKEY}/full/J6_sunlit.json"
    output: "heidelberg/checks/J6_sunlit.json"
    shell: "cp {input} {output}"


rule costs:
    # the budget's numbers from the dry run's timings, the scan copies, the products' size, and the calibration's run
    input: log="work/from_s3/jobs/ll_dry2/snakemake.log", scans="work/from_s3/scans", prod="work/from_s3/runs/dry/product_bytes.txt",
           cal=f"work/from_s3/calibration/{CKEY}/full.log", code="heidelberg/producers/costs.py", prices="aws/prices.csv"
    output: "heidelberg/checks/costs.json"
    shell: "{PY} -m heidelberg.producers.costs {input.log} {input.scans} {input.prod} {input.cal} {output} > {D}/costs.log"


rule forecast:
    # what the main sample should hold if no source is a glint, and how far the shadow moves during an exposure
    input: rs="work/record_set.csv", j3="heidelberg/checks/J3_numbers.json", code="heidelberg/producers/forecast_numbers.py",
           model="heidelberg/lib/shadow_model.py"
    output: "heidelberg/checks/forecast.json"
    shell: "{PY} -m heidelberg.producers.forecast_numbers {input.rs} {input.j3} {output}"


rule run_numbers:
    # the different-night controls' bound on static sources (the dry run's part two) and their number
    input: sb=f"{DRY}/part2/open/static_bound.json", ctrl="heidelberg/record_control_pairs.csv"
    output: "heidelberg/checks/run_numbers.json"
    shell: "{PY} -c \"import json, sys; b = json.load(open('{input.sb}')); n = sum(1 for _ in open('{input.ctrl}')) - 1; "
           "json.dump(dict(pilot=dict(static_upper=b['limit_95'], static_upper_per_pair=b['per_pair'], static_pairs=b['pairs']), "
           "sample=dict(control_pairs=n)), open('{output}', 'w'), indent=1)\""
