"""The public release: a deposit directory for Zenodo with the code at the release commit (git
archive), the input manifest, the registered sample manifest and coverage map, the transient sample (the canonical
table), the open nulls, the unsealed results with their checksums and the seal record (U3), the verdict, the rate, and
PROVENANCE.json naming the commit and the SHA-256 of every file in the deposit. The upload is made by hand, once approved.
With PRODUCTS_DIR the measured products are packed as N_ARCHIVES tar archives of whole pairs (Zenodo: at most 100
files per record, and 50 GB unless a quota increase, up to 200 GB, is granted).
python -E -s tools/release.py COMMIT REGISTERED_DIR OPEN_DIR UNSEALED_DIR RESULTS_DIR OUT_DIR [PRODUCTS_DIR]"""
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from heidelberg.stages.sample_checksum import canonical

COMMIT, REG, OPEN, UNSEALED, RES, OUT = sys.argv[1:7]; PRODUCTS = sys.argv[7] if len(sys.argv) > 7 else None
N_ARCHIVES, MAX_FILES = 50, 100
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
shutil.rmtree(OUT, ignore_errors=True); os.makedirs(f"{OUT}/results")
subprocess.run(["git", "-C", ROOT, "archive", "--format=tar.gz", f"--output={OUT}/code.tar.gz", COMMIT], check=True)
shutil.copy(f"{ROOT}/heidelberg/inputs.sha256", f"{OUT}/inputs.sha256")
os.makedirs(f"{OUT}/registered")                      # every registered file (rec and ctrl), with the records
for f in sorted(glob.glob(f"{REG}/**", recursive=True)):
    if os.path.isfile(f) and not f.endswith(".log"):
        r_ = os.path.relpath(f, REG); os.makedirs(os.path.dirname(f"{OUT}/registered/{r_}"), exist_ok=True); shutil.copy(f, f"{OUT}/registered/{r_}")
shutil.copy(f"{UNSEALED}/MANIFEST.sha256", f"{OUT}/unsealed_MANIFEST.sha256")    # the unsealed files' checksums (U3)
open(f"{OUT}/transient_sample.csv", "w").write(canonical(f"{UNSEALED}/opened/select"))
for f in glob.glob(f"{OPEN}/null_*.json") + glob.glob(f"{UNSEALED}/shadow/*/sealed_*.json") + glob.glob(f"{RES}/*.json"):
    shutil.copy(f, f"{OUT}/results/{os.path.basename(f)}")
if PRODUCTS:                                            # whole pairs per archive, in order, so a pair is never split
    import tarfile
    pairs = sorted({os.path.basename(f).split("_")[0].split(".")[0] for f in glob.glob(f"{PRODUCTS}/B*")})
    os.makedirs(f"{OUT}/products")
    for k, grp in enumerate(pairs[i::N_ARCHIVES] for i in range(N_ARCHIVES)):
        if not grp:
            continue
        with tarfile.open(f"{OUT}/products/products_{k:02d}.tar", "w") as tf:
            for p in grp:
                for f in sorted(glob.glob(f"{PRODUCTS}/{p}[._]*")):
                    tf.add(f, arcname=os.path.basename(f))
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
files = sorted(os.path.relpath(f, OUT) for f in glob.glob(f"{OUT}/**", recursive=True) if os.path.isfile(f))
json.dump(dict(commit=subprocess.run(["git", "-C", ROOT, "rev-parse", COMMIT], capture_output=True, text=True).stdout.strip(),
               files={f: sha(f"{OUT}/{f}") for f in files}), open(f"{OUT}/PROVENANCE.json", "w"), indent=1)
assert len(files) + 1 <= MAX_FILES, f"{len(files) + 1} files: Zenodo takes at most {MAX_FILES} per record"
print(f"{OUT}: {len(files)} files and PROVENANCE.json")
