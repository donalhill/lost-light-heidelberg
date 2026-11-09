"""Part one stops unless every pair's plate a has an asteroid ephemeris: the asteroid veto skips a plate without one, so a
missing ephemeris would let known asteroids through on that pair alone.
python -E -s tools/check_ephemerides.py EPHEMERIS_PLATES.csv PAIRS.csv [PAIRS.csv ...]"""
import csv
import sys

have = {r["plate"] for r in csv.DictReader(open(sys.argv[1]))}
need = [r["a"] for f in sys.argv[2:] for r in csv.DictReader(open(f))]
miss = sorted(set(need) - have)
if miss:
    sys.exit(f"{len(miss)} of {len(set(need))} plates a have no asteroid ephemeris, e.g. {', '.join(miss[:5])}")
print(f"every plate a has an asteroid ephemeris ({len(set(need))} plates)")
