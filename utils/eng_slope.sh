#!/bin/bash

# if user specify too few arguments, print help
if [[ $# -lt 1 ]]; then
  echo "Get energy slope from .cpa file using MACSIMUS plot fit"
  echo "Returns slope value and its stderr"
  echo "Call by:"
  echo "    $0 CPAFILE [ENERGYCOL(default=1) TIMECOL(default=7)]"
  exit
fi

if [[ $# -lt 3 ]]; then
  engcol="c1"
  timecol="c7"
else
  engcol="c$2"
  timecol="c$3"
fi

tabproc $timecol $engcol < $1 | plot -b -INtQ  - ":a+b*A" | fgrep "b=" | tail -n1 | sed 's|b=||' | tabproc "A" "B"