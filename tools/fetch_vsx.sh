#!/bin/bash
# VSX (VizieR B/vsx/vsx) in strips of declination into DIR/strips/LO_HI.tsv: each strip from the first mirror that
# answers (NAOJ over plain HTTP, whose HTTPS fails; CDS; CfA), three tries each; a strip no mirror serves whole is
# halved, down to a quarter of a degree; strips already present are kept, so a rerun resumes.
#   bash tools/fetch_vsx.sh DIR
set -u
D=$1; mkdir -p $D/strips
MIRRORS="http://vizier.nao.ac.jp https://vizier.cds.unistra.fr https://vizier.cfa.harvard.edu"
one() {   # one strip, lo to hi degrees; 0 if fetched
  local lo=$1 hi=$2 f=$D/strips/${1}_${2}.tsv m
  [ -s $f ] && return 0
  for m in $MIRRORS; do
    curl -sS --fail --retry 2 --retry-delay 10 --max-time 900 -o $f.part "$m/viz-bin/asu-tsv?-source=B/vsx/vsx&-out.max=unlimited&-out=OID,Name,_RAJ2000,_DEJ2000&DEJ2000=$lo..$hi" \
      && grep -q "^OID	" $f.part && tail -c 200 $f.part | grep -q "[0-9]" && mv $f.part $f && echo "$lo $hi $m" >> $D/strips.log && return 0
    rm -f $f.part
  done
  return 1
}
strip() {   # a strip, halved when no mirror serves it whole
  local lo=$1 hi=$2
  one $lo $hi && return 0
  local mid=$(python3 -c "print(($lo + $hi) / 2)")
  python3 -c "import sys; sys.exit(0 if $hi - $lo > 0.3 else 1)" || { echo "strip $lo..$hi: no mirror answered"; return 1; }
  strip $lo $mid && strip $mid $hi
}
for lo in $(seq -90 2 88); do strip $lo $((lo + 2)) || exit 1; done
