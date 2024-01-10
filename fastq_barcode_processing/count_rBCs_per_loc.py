#!/usr/bin/python

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description = 'Counts rBCs from each LP')
parser.add_argument('map_bcs', help = 'mapping bc files used for correction')
parser.add_argument('corrected_bcs', help = 'path to barcodes that were corrected')
parser.add_argument('output', help = 'output basename')
parser.add_argument('-b', nargs = '*', help = 'triple barcode file output from bc3_counts.py')
args = parser.parse_args()

def get_corrected_barcodes(filename):
    
    corrected_bcs = {}
    with open(filename, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip('\n').split('\t')
            corrected_bcs[line[1]] = line[2]

    return corrected_bcs

all_corrected_bcs = {}

for bc_file in args.b:
    corrected_bcs_file = args.corrected_bcs + bc_file.split('/')[1].split('.')[0][:-6] + 'grouped_counts_matched_corrected.txt' 
    all_corrected_bcs[bc_file] = get_corrected_barcodes(corrected_bcs_file)

map_bcs = {}

with open(args.map_bcs, 'r') as f:
    header = f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        if len(line[3]) == 16:
            map_bcs[line[3]] = line[0] + ':' + line[1] + '_' + line[-1]

all_rBCs = defaultdict(list)            

for bc_file in args.b:
    with open(bc_file, 'r') as f:
        header = f.readline()
        current_corrected_bcs = all_corrected_bcs[bc_file]
        for line in f:
            rBC, iBC, gBC, count = line.strip('\n').split('\t')
            if gBC in map_bcs.keys():
                all_rBCs[(iBC, gBC)].append(rBC)
            elif gBC in current_corrected_bcs.keys():
                corrected_gBC = current_corrected_bcs[gBC]
                all_rBCs[(iBC, corrected_gBC)].append(rBC)

rBC_tally = {k: len(set(v)) for k, v in all_rBCs.items()}

with open(args.output + '_rBC_counts.txt', 'w') as f:
    for loc, rBC_count in rBC_tally.items():
        f.write('\t'.join(loc) + '\t' + str(rBC_count) + '\n')