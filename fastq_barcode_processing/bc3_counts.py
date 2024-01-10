#!/usr/bin/python

import sys
import re
import argparse
import os
from collections import defaultdict
import gzip

parser = argparse.ArgumentParser(description = 'Obtain barcodes from fastq files')
parser.add_argument('seq_R1', help = 'sequence before R1 barcodes')
parser.add_argument('seq_R2', help = 'sequence before R2 barcodes')
parser.add_argument('fastq_R1', help = 'input fastq file R1')
parser.add_argument('fastq_R2', help = 'input fastq file R2')
parser.add_argument('ins_barcodes', help = 'insulator barcodes')
parser.add_argument('-o', '--output', help = 'output file basename')
parser.add_argument('-i', '--ignore', action='store_true', help = 'ignore the rBC and just count iBCs and gBCs')
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

# set output files basename
# defaults to empty string if no input given
if args.output is not None:
    basename = args.output + '_'
else:
    basename = ''

d_ins_barcodes = {}

with open(args.ins_barcodes, 'r') as f:
    header = f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        d_ins_barcodes[line[1]] = line[0]

R1_discarded = open(basename + 'discarded_R1', 'w')
R2_discarded = open(basename + 'discarded_R2', 'w')

R1_search_seq = re.compile(args.seq_R1 + r'([ATCGN]+$)')
R2_search_seq = re.compile(args.seq_R2 + r'([ATCGN]+$)')

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

barcodes = defaultdict(int)
grouped_barcodes = defaultdict(int)
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
        random_bc = R1_match.group(1)[:12]
        iBC = R1_match.group(1)[32:42]
        gBC = R2_match.group(1)[:16]
        if iBC in d_ins_barcodes.keys():
            barcodes[(random_bc, d_ins_barcodes[iBC], rev_comp(gBC))] += 1
            grouped_barcodes[(d_ins_barcodes[iBC], rev_comp(gBC))] += 1
        else:
            discard_count += 1
            R1_discarded.write(R1_record['seq'] + '\n')
            R2_discarded.write(R2_record['seq'] + '\n')
    
    else:
        discard_count += 1
        R1_discarded.write(R1_record['seq'] + '\n')
        R2_discarded.write(R2_record['seq'] + '\n')

if args.ignore == True:
    with open(args.output + '_grouped_counts.tsv', 'w') as f:
        f.write('iBC\tgBC\tcount\n')
        for barcodes, count in grouped_barcodes.items():
            f.write('\t'.join(barcodes) + '\t' + str(count) + '\n')
else:
    with open(args.output + '_counts.tsv', 'w') as f:
        f.write('rBC\tiBC\tgBC\tcount\n')
        for barcodes, count in barcodes.items():
            if count > 0:
                f.write('\t'.join(barcodes) + '\t' + str(count) + '\n')

print('Proportion of reads discarded: ' + str(discard_count/total_count))
R1_discarded.close()
R2_discarded.close()
