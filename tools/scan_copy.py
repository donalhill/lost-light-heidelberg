"""Whether a scan's compressed copy restores it (requirement R10): the restored file, read from standard input (funpack
-S), has the scan's pixel bytes (the data array, NAXIS1 x NAXIS2 x ... x |BITPIX|/8 bytes) byte for byte and the scan's
header cards, none lost and none added, apart from the integrity cards (CHECKSUM, DATASUM), which some funpack versions
write with the time of writing and are left out on both sides. The padding that fills the data unit's last block should
be zeros and is outside the image; a few of GAVO's scans carry non-zero bytes there, which a copy restores as zeros
(Donal, 2 November 2026, #4). The copy then still restores the scan, and the count of those bytes is printed for its
record ("padding_nonzero N of M"); any other difference in the padding fails. Streamed in blocks, so neither file is
held in memory.
funpack -S SCAN.fits.fz | python -E -s tools/scan_copy.py SCAN.fits      exits 1, saying why, unless it does"""
import hashlib
import sys
from math import prod

BLOCK = 2880                                              # a FITS block: 36 cards of 80 bytes
STAMPS = {"CHECKSUM", "DATASUM"}


def read(f):
    """(the header's cards, the SHA-256 of the pixel bytes, the padding after them) of a stream holding one FITS image."""
    cards = []
    while True:
        b = f.read(BLOCK)
        if len(b) < BLOCK:
            raise ValueError("the header ends before its END card")
        new = [b[i:i + 80].decode("ascii") for i in range(0, BLOCK, 80)]
        cards += new
        if any(c.rstrip() == "END" for c in new):
            break
    end = next(i for i, c in enumerate(cards) if c.rstrip() == "END")
    kv = {c[:8].strip(): c[10:].split("/")[0].strip() for c in cards[:end] if c[8:10] == "= "}
    left = abs(int(kv["BITPIX"])) // 8 * prod(int(kv[f"NAXIS{i}"]) for i in range(1, int(kv["NAXIS"]) + 1))
    h = hashlib.sha256()
    while left:
        chunk = f.read(min(left, BLOCK * 1024))
        if not chunk:
            raise ValueError("the data end before the last pixel")
        h.update(chunk); left -= len(chunk)
    keep = [c for c in cards[:end] if c.strip() and c[:8].strip() not in STAMPS]
    return keep, h.hexdigest(), f.read()


def problems(scan, restored):
    """(what fails, a note for the record)"""
    (ca, ha, pa), (cb, hb, pb) = scan, restored
    bad = [] if ha == hb else ["the pixels differ"]
    lost, added = [c.strip() for c in ca if c not in cb], [c.strip() for c in cb if c not in ca]
    bad += ([f"cards lost: {lost[:3]}"] if lost else []) + ([f"cards added: {added[:3]}"] if added else [])
    note = ""
    if pa != pb:
        if len(pa) == len(pb) and not any(pb):
            note = f"padding_nonzero {sum(1 for x in pa if x)} of {len(pa)}"
        else:
            bad.append("the padding differs")
    return bad, note


if __name__ == "__main__":
    with open(sys.argv[1], "rb") as f:
        bad, note = problems(read(f), read(sys.stdin.buffer))
    if bad:
        sys.exit("the copy does not restore the scan: " + "; ".join(bad))
    if note:
        print(note)
