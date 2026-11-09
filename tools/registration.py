"""The registration record: REG_DIR/registration_record.json, committed once the checksums
are posted as the organisers instruct, with the posting's URL, its time (UTC), and the SHA-256 of each registered file
(the coverage map, the footprints, the sample manifest, the transient sample's checksum, the locked workspace's hash,
the clock's summary and correction, the plate times used, and registered.sha256, the list of every registered file).
python -E -s tools/registration.py write REG_DIR URL POSTED_UTC
python -E -s tools/registration.py check REG_DIR   exits 1 unless the record exists and every file matches it
python -E -s tools/registration.py seal REG_DIR SEALED.tle URL POSTED_UTC   the seal record (U3): REG_DIR/seal_record.json,
  committed once part two's encrypted file's checksum is posted, with the round and chain it is locked to
Part two of the record run runs check first and stops on failure; the unlock (tools/unseal.py) reads the seal record."""
import hashlib
import json
import os
import sys

FILES = ("coverage_map.csv", "coverage_footprints.json", "sample_manifest.csv", "sample_sha256.txt", "workspace.tgz.tle.sha256", "clock_summary.json",
         "clock_correction.json", "plate_times.csv", "../registered.sha256")   # the last lists every registered file, rec and ctrl
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()


def check(reg):
    f = f"{reg}/registration_record.json"
    if not os.path.exists(f):
        return ["no registration record"]
    r = json.load(open(f)); bad = [k for k in ("url", "posted_utc", "files") if not r.get(k)]
    for n in FILES:
        if not os.path.exists(f"{reg}/{n}"):
            bad.append(f"missing {n}")
        elif r.get("files", {}).get(n) != sha(f"{reg}/{n}"):
            bad.append(f"{n} differs from the record")
    return bad


if __name__ == "__main__":
    cmd, reg = sys.argv[1:3]
    if cmd == "seal":
        from unseal import header
        tle = sys.argv[3]; rnd, chain = header(tle)
        json.dump(dict(url=sys.argv[4], posted_utc=sys.argv[5], sha256=sha(tle), round=rnd, chain=chain), open(f"{reg}/seal_record.json", "w"), indent=1)
        print(f"{reg}/seal_record.json written (round {rnd})")
    elif cmd == "write":
        json.dump(dict(url=sys.argv[3], posted_utc=sys.argv[4], files={n: sha(f"{reg}/{n}") for n in FILES}), open(f"{reg}/registration_record.json", "w"), indent=1)
        print(f"{reg}/registration_record.json written")
    else:
        bad = check(reg)
        print("registration as recorded" if not bad else "registration: " + "; ".join(bad)); sys.exit(1 if bad else 0)
