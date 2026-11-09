"""The analysis's parameters and the compute settings, read once (requirement W2). config/analysis.json holds every
choice that can change a result, grouped by stage, each entry with its value and unit; config/compute.json holds what
cannot (parallelism, memory, timeouts). The code reads a value as A.section.key (analysis) or C.section.key (compute).
HD_CONFIG names another analysis file (a test's reduced toy counts, for one); PATH and SHA say which file was read, for
the stages' manifests. Lists are returned as tuples, so no code can change a parameter in place. DERIVED holds the
values measured on the tuning set (heidelberg/derived.json, written by heidelberg/stages/derived.py), which are not
choices and are typed nowhere."""
import hashlib
import json
import os
from types import SimpleNamespace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read(path):
    """The sections of a configuration file as namespaces of values, and every entry by its key (unique in the file)."""
    flat, secs = {}, {}
    for s, entries in json.load(open(path)).items():
        if s.startswith("_"):
            continue
        vals = {}
        for k, e in entries.items():
            if k in flat:
                raise KeyError(f"{k} is defined twice in {path}")
            v = e["value"]; v = tuple(v) if isinstance(v, list) else v
            flat[k] = vals[k] = v
        secs[s] = SimpleNamespace(**vals)
    return SimpleNamespace(**secs), flat


PATH = os.environ.get("HD_CONFIG", f"{ROOT}/config/analysis.json")
SHA = hashlib.sha256(open(PATH, "rb").read()).hexdigest()
A, FLAT = read(PATH)
C, CFLAT = read(f"{ROOT}/config/compute.json")


class Derived(dict):
    """The values measured on the tuning set; a missing one names the stage that makes it (the stage itself imports this
    module, so the file may be absent before its first run)."""
    def __missing__(self, key):
        raise KeyError(f"{key} is not in heidelberg/derived.json: run python -m heidelberg.stages.derived on the tuning selection (Snakefile rule derived)")


DERIVED_PATH = f"{ROOT}/heidelberg/derived.json"
DERIVED = Derived({k: e["value"] for k, e in json.load(open(DERIVED_PATH)).items() if not k.startswith("_")} if os.path.exists(DERIVED_PATH) else {})


def facts():
    """The cited facts (heidelberg/facts.json) by key: values no code defines, each with its source."""
    return {k: e["value"] for k, e in json.load(open(f"{ROOT}/heidelberg/facts.json")).items() if not k.startswith("_")}


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "get":
        # python -E -s -m heidelberg.config get NAME: one value, for a job script; NAME is analysis.KEY, compute.KEY,
        # derived.KEY, or facts.KEY (a list printed space-separated)
        src, key = sys.argv[2].split(".", 1)
        v = {"analysis": FLAT, "compute": CFLAT, "derived": DERIVED}[src][key] if src != "facts" else facts()[key]
        print(" ".join(map(str, v)) if isinstance(v, (list, tuple)) else v); sys.exit(0)
    # python -E -s -m heidelberg.config override OUT KEY=VALUE ...: a copy of the analysis configuration with the named
    # values replaced (JSON values), for a test's or a quick run's reduced toys; point HD_CONFIG at OUT
    cmd, out, *kv = sys.argv[1:]
    assert cmd == "override", cmd
    raw = json.load(open(PATH)); changes = dict(x.split("=", 1) for x in kv)
    for s, entries in raw.items():
        if not s.startswith("_"):
            for k in entries:
                if k in changes:
                    entries[k]["value"] = json.loads(changes.pop(k))
    assert not changes, f"no such parameters: {sorted(changes)}"
    raw["_about"] = f"{raw['_about']} OVERRIDE of {os.path.relpath(PATH, ROOT)}: {' '.join(kv)}"
    json.dump(raw, open(out, "w"), indent=1)
