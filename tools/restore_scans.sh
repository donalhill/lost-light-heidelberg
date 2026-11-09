#!/bin/bash
# Restore the scan copies from S3 Glacier Deep Archive before part one (bulk tier: ready within 48 hours; the restored
# copies stay readable for DAYS days), or report how many are restored.
#   tools/restore_scans.sh request [DAYS]     ask for every object under s3://BUCKET/PREFIX/ (default heidelberg_scans)
#   tools/restore_scans.sh status             count objects restored, in progress, and not requested
# Part one starts only when status reports every object restored. PREFIX can be set in the environment.
set -euo pipefail
D=$(cd "$(dirname "$0")/.." && pwd); source "$D/aws/local.env"
P="--profile $AWS_PROFILE_NAME --region $AWS_REGION_NAME"; PREFIX=${PREFIX:-heidelberg_scans}
keys() { aws s3api list-objects-v2 $P --bucket "$BUCKET" --prefix "$PREFIX/" --query 'Contents[].Key' --output text | tr '\t' '\n' | grep -v '^None$' || true; }
case "${1:-}" in
  request)
    DAYS=${2:-14}; n=0
    for k in $(keys); do
      aws s3api restore-object $P --bucket "$BUCKET" --key "$k" \
        --restore-request "{\"Days\": $DAYS, \"GlacierJobParameters\": {\"Tier\": \"Bulk\"}}" 2>/tmp/restore_err || \
        grep -q RestoreAlreadyInProgress /tmp/restore_err || { cat /tmp/restore_err >&2; exit 1; }
      n=$((n + 1))
    done
    echo "restore requested for $n objects under $PREFIX/ ($DAYS days, bulk)";;
  status)
    done_=0; going=0; none=0
    for k in $(keys); do
      r=$(aws s3api head-object $P --bucket "$BUCKET" --key "$k" --query Restore --output text)
      case "$r" in *'ongoing-request="false"'*) done_=$((done_ + 1));; *'ongoing-request="true"'*) going=$((going + 1));; *) none=$((none + 1));; esac
    done
    echo "$PREFIX/: restored $done_, in progress $going, not requested $none";;
  *) sed -n 2,6p "$0"; exit 1;;
esac
