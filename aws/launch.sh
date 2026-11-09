#!/bin/bash
# Launch one job: a Snakemake target run in the project's image on one EC2 machine, which terminates when the target
# finishes or at its time cap, whichever comes first.
#   aws/launch.sh NAME TYPE CAP_MIN ondemand|spot DISK_GB TARGET [INPUT ...]
# INPUT: paths under work/ (as in heidelberg/inputs.sha256) the job reads; the machine copies them from the S3 mirror
# (s3://BUCKET/inputs/PATH) and checks them against the manifest before anything runs. The repository goes up as a git
# bundle of HEAD, which must be committed and clean. The image is the one named in aws/image.txt (URI@digest). Outputs
# and logs go to s3://BUCKET/jobs/NAME/ as the job runs (aws/job.sh). Refused: a worst case (price times cap) over $15
# unless APPROVED_COST covers it, an unpriced type, or a machine not set to terminate on shutdown.
set -euo pipefail
D=$(cd "$(dirname "$0")/.." && pwd)
NAME=$1 TYPE=$2 CAP=$3 MARKET=$4 DISK=$5 TARGET=$6; shift 6
source "$D/aws/local.env"
P="--profile $AWS_PROFILE_NAME --region $AWS_REGION_NAME"; B=$BUCKET; S3=s3://$B/jobs/$NAME
[ "$MARKET" = ondemand ] || [ "$MARKET" = spot ] || { echo "market must be ondemand or spot" >&2; exit 1; }
PRICE=$(awk -F, -v t="$TYPE" '$1 == t {print $2}' "$D/aws/prices.csv"); [ -n "$PRICE" ] || { echo "no price for $TYPE in aws/prices.csv" >&2; exit 1; }
WORST=$(awk -v p="$PRICE" -v c="$CAP" 'BEGIN {printf "%.2f", p * c / 60}')
awk -v w="$WORST" -v a="${APPROVED_COST:-15}" 'BEGIN {exit !(w <= a)}' || { echo "worst case \$$WORST is over \$${APPROVED_COST:-15}: set APPROVED_COST if Donal approved it" >&2; exit 1; }
[ -z "$(git -C "$D" status --porcelain | cat)" ] || { echo "the repository has uncommitted changes" >&2; exit 1; }
ACC=$(aws sts get-caller-identity $P --query Account --output text); IMAGE=$ACC.dkr.ecr.$AWS_REGION_NAME.amazonaws.com/$(cat "$D/aws/image.txt"); COMMIT=$(git -C "$D" rev-parse HEAD)
git -C "$D" bundle create /tmp/ll_repo.bundle HEAD 2>/dev/null
aws s3 cp $P --quiet /tmp/ll_repo.bundle $S3/in/repo.bundle
printf "commit %s\nimage %s\ntarget %s\ninputs %s\n" "$COMMIT" "$IMAGE" "$TARGET" "$*" | aws s3 cp $P --quiet - $S3/in/PROVENANCE.txt
VPC=$(aws ec2 describe-vpcs $P --filters Name=is-default,Values=true --query "Vpcs[0].VpcId" --output text)
SG=$(aws ec2 describe-security-groups $P --filters Name=group-name,Values=$SECURITY_GROUP Name=vpc-id,Values=$VPC --query "SecurityGroups[0].GroupId" --output text)
AMI=$(aws ssm get-parameter $P --name /aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id --query Parameter.Value --output text)
UD=$(printf '#!/bin/bash\nshutdown -h +%s\nexport BUCKET=%s NAME=%s IMAGE=%q TARGET=%q INPUTS=%q REGION=%s\n' \
  "$CAP" "$B" "$NAME" "$IMAGE" "$TARGET" "$*" "$AWS_REGION_NAME"; sed 1d "$D/aws/job.sh")
MKT=(); [ "$MARKET" = spot ] && MKT=(--instance-market-options 'MarketType=spot,SpotOptions={SpotInstanceType=one-time,InstanceInterruptionBehavior=terminate}')
ID=$(aws ec2 run-instances $P --image-id "$AMI" --instance-type "$TYPE" --key-name "$KEY_NAME" --security-group-ids "$SG" \
  --iam-instance-profile Name=$INSTANCE_ROLE --instance-initiated-shutdown-behavior terminate ${MKT[@]+"${MKT[@]}"} \
  --metadata-options HttpTokens=required,HttpPutResponseHopLimit=2 \
  --block-device-mappings "DeviceName=/dev/sda1,Ebs={VolumeSize=$DISK,VolumeType=gp3,Throughput=500,Iops=6000,DeleteOnTermination=true}" \
  --tag-specifications "ResourceType=instance,Tags=[{Key=project,Value=$PROJECT_TAG},{Key=Name,Value=$PROJECT_TAG-ll-$NAME}]" "ResourceType=volume,Tags=[{Key=project,Value=$PROJECT_TAG}]" \
  --user-data "$UD" --query "Instances[0].InstanceId" --output text)
aws ec2 wait instance-running $P --instance-ids "$ID"
SD=$(aws ec2 describe-instance-attribute $P --instance-id "$ID" --attribute instanceInitiatedShutdownBehavior --query InstanceInitiatedShutdownBehavior.Value --output text)
[ "$SD" = terminate ] || { echo "$ID is not set to terminate on shutdown: terminating it" >&2; aws ec2 terminate-instances $P --instance-ids "$ID" > /dev/null; exit 1; }
echo "$(date -u +%Y-%m-%dT%H:%M:%SZ),$NAME,$TYPE,$MARKET,$CAP,$WORST,$ID,${COMMIT:0:7}" >> "$D/local/ledger.csv"
echo "$ID $TYPE $MARKET cap ${CAP} min, worst case \$$WORST, commit ${COMMIT:0:7}"
