#!/bin/bash
# The body of every job (aws/launch.sh puts it in the machine's user data, after the time cap and the exports BUCKET,
# NAME, IMAGE, TARGET, INPUTS, and REGION). Runs as root on Ubuntu 24.04.
# Paths: an input under work/ comes from the inputs mirror (s3://BUCKET/inputs/PATH) and is checked against the
# manifest; an input under work/ll/ comes from s3://BUCKET/ll/ (the outputs of earlier jobs, and of this one if it is
# resumed). Everything the target writes under work/ll/ goes back to s3://BUCKET/ll/ every minute and at the end, so a
# job stopped by its cap or a spot interruption resumes from what it uploaded. Logs go to s3://BUCKET/jobs/NAME/.
set -u
trap 'shutdown -h now' EXIT                                      # whatever ends the job ends the machine
J=s3://$BUCKET/jobs/$NAME; L=/root/progress.txt
apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq docker.io amazon-ecr-credential-helper git > /dev/null 2>&1
mkdir -p /root/.docker && echo "{\"credHelpers\": {\"${IMAGE%%/*}\": \"ecr-login\"}}" > /root/.docker/config.json
docker pull -q "$IMAGE" > /dev/null || { echo "image pull failed"; exit 1; }
R() { docker run --rm --network host -v /repo:/repo -w /repo -e HOME=/tmp -e AWS_DEFAULT_REGION=$REGION -e BUCKET=$BUCKET "$IMAGE" "$@"; }
log() { echo "$(date -u +%H:%M:%S) $*" >> $L; R aws s3 cp --quiet /repo/progress.txt $J/progress.txt 2>/dev/null; }
mkdir -p /repo && cd /repo
R aws s3 cp --quiet $J/in/repo.bundle /repo/repo.bundle && git clone -q /repo/repo.bundle /tmp/clone && cp -a /tmp/clone/. /repo/ && rm -rf /tmp/clone /repo/repo.bundle
L=/repo/progress.txt; log "started: image $IMAGE, commit $(git -C /repo rev-parse --short HEAD), target $TARGET"
for p in $INPUTS; do
  case $p in work/ll/*) src=s3://$BUCKET/ll/${p#work/ll/} ;; *) src=s3://$BUCKET/inputs/$p ;; esac
  if awk -v p="$p" 'substr($0, 67) == p {f = 1} END {exit !f}' heidelberg/inputs.sha256; then cp_=cp; else cp_=sync; fi   # a file or a folder
  R aws s3 $cp_ --quiet "$src" "/repo/$p" || { log "STOPPED: input $p"; exit 1; }
done
pat=$(printf '%s\n' $INPUTS | grep -v '^work/ll/')
if [ -n "$pat" ]; then
  awk -v pats="$pat" 'BEGIN {n = split(pats, P, "\n")} {f = substr($0, 67); for (i = 1; i <= n; i++) if (f == P[i] || index(f, P[i] "/") == 1) {print; break}}' \
    heidelberg/inputs.sha256 > /tmp/m.sha256
  (cd /repo && sha256sum -c --quiet /tmp/m.sha256) > /repo/inputs_check.txt 2>&1 \
    || { R aws s3 cp --quiet /repo/inputs_check.txt $J/inputs_check.txt; log "STOPPED: inputs differ from the manifest (inputs_check.txt)"; exit 1; }
  log "inputs checked: $(wc -l < /tmp/m.sha256) files"
fi
mkdir -p /repo/work/ll
( while sleep 60; do R aws s3 sync --quiet /repo/work/ll s3://$BUCKET/ll; R aws s3 sync --quiet /repo/work/dag s3://$BUCKET/jobs/$NAME/dag --exclude "*" --include "*.log" --include "*.txt"; done ) &
SYNC=$!
R snakemake -j "$(nproc)" --resources mem_gb=$(( $(free -g | awk '/Mem/ {print $2}') - 8 )) --keep-going $TARGET > /repo/snakemake.log 2>&1; rc=$?
kill $SYNC 2>/dev/null
R aws s3 sync --quiet /repo/work/ll s3://$BUCKET/ll
R aws s3 sync --quiet /repo/work/dag s3://$BUCKET/jobs/$NAME/dag --exclude "*" --include "*.log" --include "*.txt" --include "*.json" --include "*.ok"
R aws s3 cp --quiet /repo/snakemake.log $J/snakemake.log
log "$([ $rc = 0 ] && echo done || echo "FAILED: snakemake exit $rc")"
