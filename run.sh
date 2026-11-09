#!/bin/bash
# Run a command in the project's image from the repository root, e.g. ./run.sh snakemake -j4 draft. The previous
# repository, where present, is mounted read-only at /ref (the reference of the rebuild's equality checks).
REF="$HOME/Documents/Code/lost-light-heidelberg-v1.nosync"
exec docker run --rm --platform linux/amd64 -v "$PWD":/repo $([ -d "$REF" ] && echo -v "$REF":/ref:ro) \
  -w /repo -e HOME=/tmp lost-light:dev "$@"
