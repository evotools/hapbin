#!/bin/bash
#Usage: processdata.selscan.sh [transposed source .hap] [# of individuals]

sed -n "$(awk '{printf $0"p; ";}' numbers.$2.txt)" < $1 > chr22.$2.t.hap
