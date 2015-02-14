#!/bin/bash

#$(cat numbers.$1.txt | paste -sd ",")
cut -d" " -f $(cat numbers.$2.txt | paste -sd ",") $1 > chr22.$2.hap

