#!/usr/bin/python

import sys
import re
import argparse
import os
from collections import defaultdict
import gzip

parser = argparse.ArgumentParser(
    description = 'Demultiplexes Illumina sequencing fastq files')
parser.add_argument(
    'seq_R1', help = 'primer sequence before genome')
parser.add_argument(
    'seq_R2', help = 'sequence before barcode')
parser.add_argument(
    'fastq_R1', help = 'input fastq file R1')
parser.add_argument(
    'fastq_R2', help = 'input fastq file R2')
parser.add_argument(
    '-o', '--output', help = 'output file basename')
args = parser.parse_args()

def process_fastq(fastq_file):

    current_record = {}

    for name, seq, crap, quality in zip(*[iter(fastq_file)]*4):
        current_record['name'] = name.strip('\n')
        current_record['seq'] = seq.strip('\n')
        current_record['quality'] = quality.strip('\n')

        yield current_record

def get_fastq_id(fastq_name):

    return fastq_name.split(' ')[0]

# set output files basename
# defaults to empty string if no input given
if args.output is not None:
    basename = args.output + '_'
else:
    basename = ''

R1_discarded = open(basename + 'discarded_R1', 'w')
R2_discarded = open(basename + 'discarded_R2', 'w')

R1_search_seq = re.compile(args.seq_R1 + r'([ATCGN]+$)')
R2_search_seq = re.compile(args.seq_R2 + r'([ATCGN]+$)')

output_file = open(args.output + '.fastq', 'w')

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

d_random_bcs = defaultdict(int)

for R1_record, R2_record in zip(R1_reader, R2_reader):
    
    if R1_record['name'].split(' ')[0] != R2_record['name'].split(' ')[0]:
        print(R1_record['name'])
        print(R2_record['name'])
        break

    R1_match = R1_search_seq.search(R1_record['seq'])
    R2_match = R2_search_seq.search(R2_record['seq'])

    if R1_match != None and R2_match != None:
        random_bc = R2_match.group(1)[:16]
        d_random_bcs[random_bc] += 1

        R1_name = R1_record['name'].split(' ')[0] + ':' + random_bc
        R1_seq = R1_match.group(1)[:40]
        start_idx = R1_match.span()[0] + 6
        R1_quality = R1_record['quality'][start_idx:start_idx+40]

        if len(R1_seq) > 30:
            output_file.write('{name}\n{seq}\n+\n{quality}\n'.format(
                name = R1_name, seq = R1_seq, quality = R1_quality))
    
    else:
        R1_discarded.write(R1_record['seq'] + '\n')
        R2_discarded.write(R2_record['seq'] + '\n')

print(len(d_random_bcs))
output_file.close()
R1_discarded.close()
R2_discarded.close()
