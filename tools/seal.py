"""Time locks with drand quicknet (tle, in the image): a file encrypted to the round emitted at a given time, kept only if
tle then refuses to decrypt it as too early and names that round; or a locked file opened once its round is out.
python tools/seal.py lock FILE WHEN OUT.tle      WHEN: an ISO 8601 UTC time, or +SECONDS from now (a dry run)
    writes OUT.tle, OUT.tle.sha256, and OUT.tle.round (chain, round, and the time it is emitted)
python tools/seal.py open FILE.tle EXPECTED_SHA256 OUT   decrypts (retrying while the round is not yet out) to OUT
The chain, its genesis, and its period are the cited facts in heidelberg/facts.json."""
import calendar
import hashlib
import subprocess
import sys
import time

sys.path.insert(0, ".")
from heidelberg import config

F = config.facts(); CHAIN, GENESIS, PERIOD = F["drand_chain"], F["drand_genesis_s"], F["drand_period_s"]
TLE = "tle"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()


def lock(src, when, out):
    t = time.time() + float(when[1:]) if when.startswith("+") else calendar.timegm(time.strptime(when, "%Y-%m-%dT%H:%M:%SZ"))
    rnd = (int(t) - GENESIS) // PERIOD + 1
    subprocess.run([TLE, "-e", "-c", CHAIN, "-r", str(rnd), "-o", out, src], check=True)
    r = subprocess.run([TLE, "-d", "-c", CHAIN, "-o", "/dev/null", out], capture_output=True, text=True)
    msg = r.stdout + r.stderr
    if r.returncode == 0 or "too early" not in msg or str(rnd) not in msg:
        sys.exit(f"{out} did not refuse decryption as too early for round {rnd}: {msg.strip()}")
    open(out + ".sha256", "w").write(sha(out) + "\n")
    emitted = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(GENESIS + (rnd - 1) * PERIOD))
    open(out + ".round", "w").write(f"chain {CHAIN}\nround {rnd}\nemitted {emitted}\n")
    print(f"locked to round {rnd} (emitted {emitted}): {sha(out)}")


def open_(tle, expected, out, tries=40, wait=30):
    if sha(tle) != expected.strip():
        sys.exit(f"{tle} differs from its registered checksum")
    for _ in range(tries):
        r = subprocess.run([TLE, "-d", "-c", CHAIN, "-o", out, tle], capture_output=True, text=True)
        if r.returncode == 0:
            print(f"opened {tle}"); return
        time.sleep(wait)
    sys.exit(f"{tle} did not open: {(r.stdout + r.stderr).strip()}")


if __name__ == "__main__":
    cmd, *a = sys.argv[1:]
    lock(*a) if cmd == "lock" else open_(*a)
