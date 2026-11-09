"""The registration's freeze: what "this registration freezes the design, code, and settings, and records their SHA-256
checksums" means in practice. Writes FREEZE.json with the commit, the container image (aws/image.txt, by digest), the
SHA-256 of every file git tracks except the documents (code, settings, pair lists, and heidelberg/inputs.sha256, which
pins every input file; not draft/ or the README, so a DOI can be added to the text after the freeze),
and one SHA-256 over that list, the number posted with the registration. Refuses a working tree with uncommitted changes.
python3 tools/freeze.py write FREEZE.json     then commit FREEZE.json, tag the commit, and post its "freeze_sha256"
python3 tools/freeze.py check FREEZE.json     exits 1 unless every listed file still has its frozen checksum"""
import hashlib
import json
import subprocess
import sys

git = lambda *a: subprocess.run(["git", *a], capture_output=True, text=True, check=True).stdout
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()


TEXT = ("draft/", "README.md")      # the documents, not the analysis: a DOI or a wording added later changes no frozen file


def files():
    return sorted(f for f in git("ls-files").splitlines() if f and not f.startswith(TEXT))


if __name__ == "__main__":
    cmd, out = sys.argv[1:3]
    if cmd == "write":
        if git("status", "--porcelain", "--untracked-files=no").strip():
            sys.exit("uncommitted changes: commit them before freezing")
        F = {f: sha(f) for f in files() if f != out}
        listing = "".join(f"{h}  {f}\n" for f, h in F.items())
        rec = dict(commit=git("rev-parse", "HEAD").strip(), image=open("aws/image.txt").read().strip(),
                   files=F, freeze_sha256=hashlib.sha256(listing.encode()).hexdigest())
        json.dump(rec, open(out, "w"), indent=1)
        print(f"{out}: {len(F)} files, commit {rec['commit'][:12]}, freeze checksum {rec['freeze_sha256']}")
    elif cmd == "check":
        rec = json.load(open(out)); bad = [f for f, h in rec["files"].items() if sha(f) != h]
        listing = "".join(f"{h}  {f}\n" for f, h in rec["files"].items())
        if hashlib.sha256(listing.encode()).hexdigest() != rec["freeze_sha256"]:
            bad.append("the record's own checksum")
        print("\n".join(f"differs: {b}" for b in bad) or f"{len(rec['files'])} files as frozen")
        sys.exit(1 if bad else 0)
