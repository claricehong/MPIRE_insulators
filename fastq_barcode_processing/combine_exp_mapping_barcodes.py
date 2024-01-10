#!/usr/bin/python3

import argparse
from collections import defaultdict
from collections import Counter

parser = argparse.ArgumentParser(description = 'Combines mapping and expression barcodes')
parser.add_argument('exp_bcs', help = 'File with expression barcodes')
parser.add_argument('map_bcs', help = 'File with mapping barcodes')
parser.add_argument('out_file', help = 'Output file name')
parser.add_argument('-r', action = 'store_true', help = 'Indicates whether its emerald LP or recombined mScarlet counts')
parser.add_argument('-u', action = 'store_true', help = 'Indicates whether to store the uncorrected barcodes')
args = parser.parse_args()

exp_bcs = defaultdict(dict)

if args.r == True:
    with open(args.exp_bcs, 'r') as f:
        header = f.readline()
        for line in f:
            iBC, gBC, count = line.strip('\n').split('\t')
            if len(gBC) == 16:
                exp_bcs[(gBC, iBC)] = int(count)
else:
    with open(args.exp_bcs, 'r') as f:
        header = f.readline()
        for line in f:
            gBC, count = line.strip('\n').split('\t')
            if len(gBC) == 16:
                exp_bcs[gBC] = int(count)

map_bcs = {}

with open(args.map_bcs, 'r') as f:
    header = f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        if len(line[3]) == 16:
            map_bcs[line[3]] = line[0] + ':' + line[1] + '_' + line[-1]


def count_distance(str1, str2):

    hamming_distance = 0
    for base1, base2 in zip(str1, str2):
        if base1 != base2:
            hamming_distance += 1

    return hamming_distance

def process_edit_distances(counter):
    good_counts = 0
    bad_counts = 0

    total_counts = sum(counter.values())

    # Most will have matches in the 5, 6, 7 ranges so the second filter really isn't great 
    for distance, count in counter.items():
        if distance < 5:
            good_counts += count
        #elif distance > 7:
        #    bad_counts += count
    
    #total_counted = bad_counts + good_counts

    if good_counts == 1:
        return min(counter, key = counter.get)
    return 'bad'

def find_corrected_bc(hammings):
    
    return min(hammings, key = hammings.get)

corrected_exp_barcodes = defaultdict(int)
uncorrected_barcodes = defaultdict(int)

for BCs, count in exp_bcs.items():
    if args.r == True:
        bc = BCs[0]
    else:
        bc = BCs
    all_hamming = {i:count_distance(bc, i) for i in list(map_bcs.keys())}
    if 0 in all_hamming.values():
        corrected_exp_barcodes[BCs] += count
    else:
        edit_distance_counts = Counter(all_hamming.values())
        distance = process_edit_distances(edit_distance_counts)
        if distance != 'bad':
            corrected_bc = find_corrected_bc(all_hamming)
            if args.r == True:
                corrected_exp_barcodes[(corrected_bc, BCs[1])] += count
                uncorrected_barcodes[(BCs[1], bc, corrected_bc, str(distance))] += count
            else:
                corrected_exp_barcodes[corrected_bc] += count

if args.r == True:
    with open(args.out_file.split('.')[0] + '_corrected.txt', 'w') as f:
        f.write('gBC\tiBC\tcount\n')
        for bc, count in corrected_exp_barcodes.items():
            f.write('\t'.join(bc) + '\t' + str(count) + '\n')
else:
    with open(args.out_file + '_corrected.txt', 'w') as f:
        f.write('gBC\tcount\n')
        for bc, count in corrected_exp_barcodes.items():
            f.write(bc + '\t' + str(count) + '\n')

if args.u == True:
    with open(args.out_file.split('.')[0] + '_matched_corrected.txt', 'w') as f:
        f.write('iBC\tuncorrected_gBC\tcorrected_gBC\thamming_distance\tcount\n')
        for bcs, count in uncorrected_barcodes.items():
            f.write('\t'.join(bcs) + '\t' + str(count) + '\n')
