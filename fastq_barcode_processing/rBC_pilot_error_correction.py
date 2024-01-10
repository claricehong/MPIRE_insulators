#!/usr/bin/python

import sys
import re
import argparse
from collections import defaultdict
from collections import Counter
import gzip

parser = argparse.ArgumentParser(description = 'Obtain barcodes from fastq files')
parser.add_argument('seq_R1', help = 'sequence before R1 barcodes')
parser.add_argument('seq_R2', help = 'sequence before R2 barcodes')
parser.add_argument('fastq_R1', help = 'input fastq file R1')
parser.add_argument('fastq_R2', help = 'input fastq file R2')
parser.add_argument('-o', '--output', help = 'output file basename')
args = parser.parse_args()

def process_fastq(fastq_file):

    current_record = {}

    for name, seq, crap, quality in zip(*[iter(fastq_file)]*4):
        current_record['name'] = name.strip('\n')
        current_record['seq'] = seq.strip('\n')
        current_record['quality'] = quality.strip('\n')

        yield current_record

def rev_comp(seq):
    
    new_seq = ''
    conversion = {'A': 'T', 'T': 'A', 'G':'C', 'C':'G', 'N': 'N'}
    
#    if 'N' in seq:
#        return 'invalid_barcode'
#    else:
    for i in reversed(seq):
        new_seq += conversion[i]
            
    return new_seq

def get_fastq_id(fastq_name):

    return fastq_name.split(' ')[0]

def calculate_hamming_dist(rBC1, rBC2):

    hamming = 0

    for base1, base2 in zip(rBC1, rBC2):
        if base1 != base2:
            hamming += 1

    return hamming

# random_file = open('trouble.txt', 'w')

# def correct_rBCs(rBC_list):
#     # rBC_list.sort(key = lambda x: x[1], reverse = True)
#     filtered = {k: v for k, v in rBC_list.items() if v > 5}
#     sorted_list = [(k, v) for k, v in sorted(filtered.items(), key = lambda x: x[1], reverse = True)]

#     corrected_rBCs = defaultdict(int)
#     for rBC, count in sorted_list:
        # corrected = False
        # for corrected_bc in all_rBCs:
            # if calculate_hamming_dist(rBC, corrected_bc) < 3:
                # corrected_rBCs[corrected_bc] += count
                # corrected = True
                # random_file.write('\n')
                # for k, v in corrected_rBCs.items():    
                #     random_file.write(k + ': ' + str(v) + '\t')
                # break
        # if corrected == False:
        #     corrected_rBCs[rBC] += count
    #     if rBC in valid_rBCs:
    #         corrected_rBCs[rBC] += count

    # return corrected_rBCs

# set output files basename
# defaults to empty string if no input given
if args.output is not None:
    basename = args.output + '_'
else:
    basename = ''

R1_discarded = open(basename + 'discarded_R1', 'w')
R2_discarded = open(basename + 'discarded_R2', 'w')

#AGCTGTACAAGTAAGCTAGC
R1_search_seq = re.compile(args.seq_R1 + '([ATCGN]{12})')
R2_search_seq = re.compile(args.seq_R2 + '([ATCGN]{16})')
#GAAGGGCCGGCCACAACTCGAG

# open fastq files
if str(args.fastq_R1).endswith('.gz'):
    R1_open_file = gzip.open(args.fastq_R1, 'rt')
else:
    R1_open_file = open(args.fastq_R1, 'r')

if str(args.fastq_R2).endswith('.gz'):
    R2_open_file = gzip.open(args.fastq_R2, 'rt')
else:
    R2_open_file = open(args.fastq_R2, 'r')

# read fastq files
R1_reader = process_fastq(R1_open_file)
R2_reader = process_fastq(R2_open_file)

barcodes = defaultdict(dict)
discard_count = 0
total_count = 0

for R1_record, R2_record in zip(R1_reader, R2_reader):

    total_count += 1

    if R1_record['name'].split(' ')[0] != R2_record['name'].split(' ')[0]:
        print(R1_record['name'])
        print(R2_record['name'])
        break

    R1_match = R1_search_seq.search(R1_record['seq'])
    R2_match = R2_search_seq.search(R2_record['seq'])

    if R1_match != None and R2_match != None:
        random_bc = R1_match.group(1)
        gBC = R2_match.group(1)
        if random_bc in barcodes[gBC].keys():
            barcodes[gBC][random_bc] += 1
        else:
            barcodes[gBC][random_bc] = 1
    else:
        discard_count += 1
        R1_discarded.write(R1_record['seq'] + '\n')
        R2_discarded.write(R2_record['seq'] + '\n')

# all_rBCs = []
# for gBC, rBCs in barcodes.items():
#     for rBC, count in rBCs.items():
#         if count > 1:
#             all_rBCs.append(rBC)

# counted_rBCs = Counter(all_rBCs)
# valid_rBCs = [k for k, v in counted_rBCs.items() if v == 1]

# corrected_barcodes = {}

# for gBC, rBCs in barcodes.items():
#     corrected_barcodes[gBC] = correct_rBCs(rBCs)

with open(args.output + '_counts_filtered5_uniquerBC.txt', 'w') as f:
    f.write('rBC\tgBC\tcount\n')
    for gBC, rBCs in barcodes.items():
        for rBC, count in rBCs.items():
            if count > 5:
                f.write(rBC + '\t' + gBC + '\t' + str(count) + '\n')

print('Proportion of reads discarded: ' + str(discard_count/total_count))
R1_discarded.close()
R2_discarded.close()
