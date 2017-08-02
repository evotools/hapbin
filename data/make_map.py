#!/usr/bin/env python3

import sys
import argparse

def genmap_parse(line):
    l = line.split()
    # (position, genetic_position)
    return (int(l[0]), l[2])

def legend_parse(line):
    l = line.split()
    # (position, id)
    return (int(l[1]), l[0])

parser = argparse.ArgumentParser()
parser.add_argument("--chromosome", help="Create map for CHROMOSOME", required=True)
args = parser.parse_args()

chromosome = args.chromosome
with open('1000GP_Phase3_chr' + chromosome + '.legend') as legend, open('genetic_map_chr' + chromosome + '_combined_b37.txt') as genmap, open('chr' + chromosome + '.map', 'w') as outmap:
    next(legend)
    next(genmap)
    (genmap_pos, genmap_genpos) = genmap_parse(genmap.readline())
    count = 0
    for legend_line in legend:
        if len(legend_line) == 0:
            continue
        count+=1
        (pos, id) = legend_parse(legend_line)
        if pos >= genmap_pos:
            line = genmap.readline()
            if len(line) > 0:
                (genmap_pos, genmap_genpos) = genmap_parse(line)
        outmap.write("%s %s %s %s\n" % (chromosome, count, genmap_genpos, pos))
