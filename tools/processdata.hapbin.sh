#!/bin/bash
#Usage: processdata.hapbin.sh [source .hap] [# of individuals]

cut -d" " -f $(cat numbers.$2.txt | paste -sd ",") $1 > $2.hap

