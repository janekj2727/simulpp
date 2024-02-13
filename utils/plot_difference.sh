#!/bin/bash

# Universal script parsing options...
if [[ $# -lt 1 ]]; then
    echo "Plot differences between two cpa files"
    echo "Call by:"
    echo "    $0 FILE1 FILE2 [MAX_COLUMN_NO=15]"
    exit
fi

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="$1"
    POSITIONAL+=("$1") # save it in an array for later
    shift              # past argument
done

set -- "${POSITIONAL[@]}" # restore positional parameters

FILE1="${POSITIONAL[0]}"
FILE2="${POSITIONAL[1]}"
if [[ "${POSITIONAL[2]}" == "" ]]; then
    MAX_COL=15
else
    MAX_COL="${POSITIONAL[2]}"
fi

for i in $(seq 1 "$MAX_COL"); do
    mergetab "$FILE1":"$i" "$FILE2":"$i" | PLOTNAME="Column $(grep -oh " $i:\w*" $FILE1) difference" plot -:0:"B-A" 2>/dev/null
done
