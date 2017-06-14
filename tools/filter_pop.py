#!/usr/bin/env python3

import sys
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sample", help="File detailing population and group of individuals", required=True)
parser.add_argument("--select", help="Select population or group", required=True)
parser.add_argument("--hap", help="Input haplotype file", required=True)
parser.add_argument("--hapout", help="Output haplotype file", required=True)
parser.add_argument("--map", help="Input map file", required=True)
parser.add_argument("--mapout", help="Output map file", required=True)
args = parser.parse_args()

sample_filename = args.sample

sel = args.select

hap_filename = args.hap
hap_filtered_filename = args.hapout
map_filename = args.map
map_filtered_filename = args.mapout

locations = defaultdict(list)

with open(sample_filename) as sample_file:
    next(sample_file)
    pos = 0
    for line in sample_file:
        data = line.split()
        pop = data[1]
        group = data[2]
        locations[pop].append(pos)
        locations[pop].append(pos+2)
        locations[group].append(pos)
        locations[group].append(pos+2)
        pos += 4

print(locations[sel])

chrs = len(locations[sel])
t1s = " ".join(['1'] * chrs)
t0s = " ".join(['0'] * chrs)
with open(hap_filename) as hap, open(hap_filtered_filename, 'w') as hap_filtered, open(map_filename) as map, open(map_filtered_filename, 'w') as map_filtered:
    for line in hap:
        mapline = map.readline()
        filtered = " ".join([line[i] for i in locations[sel]])
        if filtered != t0s and filtered != t1s:
            hap_filtered.write(filtered + '\n')
            map_filtered.write(mapline)
