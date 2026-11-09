"""The unlock: the encrypted results opened only if the file is the one whose checksum was
posted after part two (the seal record, tools/registration.py seal) and is locked to the drand chain and round the
record names, which for the record run must be the registered round, posted before that round was emitted. The file is
then decrypted with tlock, unpacked, and every unsealed file's SHA-256 written to OUT/MANIFEST.sha256, which the
release deposits with the record.
python -E -s tools/unseal.py SEALED.tle SEAL_RECORD.json TLE CHAIN ROUND OUT
  ROUND: the registered round, or - for a dry run (the record's own round then stands)."""
import calendar
import hashlib
import json
import os
import subprocess
import sys
import tarfile
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from heidelberg import config

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()


def header(tle):
    """(round, chain) from the tlock stanza of the encrypted file's age header."""
    with open(tle, "rb") as f:
        for ln in f.read(4096).split(b"\n"):
            if ln.startswith(b"-> tlock "):
                _, _, rnd, chain = ln.decode().split()
                return int(rnd), chain
    return None, None


def problems(tle, rec, chain, rnd):
    r = json.load(open(rec)); bad = [f"the seal record has no {k}" for k in ("url", "posted_utc", "sha256", "round", "chain") if not r.get(k)]
    if bad:
        return bad
    h_round, h_chain = header(tle)
    if sha(tle) != r["sha256"]:
        bad.append("the encrypted file differs from the posted checksum")
    if (h_round, h_chain) != (r["round"], r["chain"]):
        bad.append(f"the file is locked to round {h_round} of chain {h_chain}, not the record's")
    if r["chain"] != chain:
        bad.append("the record's chain is not the pinned chain")
    if rnd != "-":
        F = config.facts(); emitted = F["drand_genesis_s"] + (int(rnd) - 1) * F["drand_period_s"]
        posted = calendar.timegm(time.strptime(r["posted_utc"], "%Y-%m-%dT%H:%M:%SZ"))
        if r["round"] != int(rnd):
            bad.append(f"the record's round {r['round']} is not the registered round {rnd}")
        if posted >= emitted:
            bad.append("the checksum was posted after the round was emitted")
    return bad


if __name__ == "__main__":
    TLE_FILE, REC, TLE, CHAIN, ROUND, OUT = sys.argv[1:7]
    bad = problems(TLE_FILE, REC, CHAIN, ROUND)
    if bad:
        sys.exit("not unsealed: " + "; ".join(bad))
    os.makedirs(OUT, exist_ok=True)
    subprocess.run([TLE, "-d", "-c", CHAIN, "-o", f"{OUT}/sealed.tgz", TLE_FILE], check=True)
    with tarfile.open(f"{OUT}/sealed.tgz") as t:
        t.extractall(OUT, filter="data")
    os.remove(f"{OUT}/sealed.tgz")
    files = sorted(os.path.relpath(os.path.join(d, n), OUT) for d, _, ns in os.walk(OUT) for n in ns if n != "MANIFEST.sha256")
    open(f"{OUT}/MANIFEST.sha256", "w").write("".join(f"{sha(f'{OUT}/{f}')}  {f}\n" for f in files))
    print(f"unsealed {len(files)} files from round {header(TLE_FILE)[0]}; their checksums in {OUT}/MANIFEST.sha256")
