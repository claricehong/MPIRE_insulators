#!/usr/bin/python3

import sys
import argparse
import os
import re
import gzip
from collections import defaultdict

parser = argparse.ArgumentParser(description = 'Gets BC-BC counts straight from demultiplexed Fastq files')
parser.add_argument('fastq', help = 'input fastq file (gzip or not gzip is fine)')
parser.add_argument('seq', help = 'sequence before labelling barcode')
parser.add_argument('output', help = 'output file basename')
parser.add_argument('-rc', '--rc', help = 'whether to reverse complement the barcode (reverse if fastq is R2)', action='store_true')
args = parser.parse_args()

def rev_comp(seq):
    
    new_seq = ''
    conversion = {'A': 'T', 'T': 'A', 'G':'C', 'C':'G', 'N': 'N'}
    
    for i in reversed(seq):
        new_seq += conversion[i]
        
    return new_seq

# set output files basename
basename = os.path.basename(args.output)

discarded = open(basename + '_discarded', 'w')

iBCs = {}

#generate regex to search sequence
full_seq = '(' + args.seq + ')'
regex1 = re.compile(full_seq + r'([ATCGN]+$)')

seq_check = re.compile(r'([ATCGN]+$)')

if str(args.fastq).endswith('.gz'):
    fastq_file = gzip.open(args.fastq, 'rt')
else:
    fastq_file = open(args.fastq, 'r')

d = defaultdict(int)
discard_count = 0
counted_reads = 0
    
for line in fastq_file:
    line = line.strip()
    one = regex1.search(line)
    if one != None:
        seq = one.group(2)
        if args.rc:
            rBC = rev_comp(seq[:16])
        else:
            rBC = seq[:16]

        if 'N' in rBC:
            pass
            discard_count += 1
        else:
            d[rBC] += 1
            counted_reads += 1

    else:
        n = seq_check.fullmatch(line)
        if n != None:
            discarded.write(line + '\n')
            discard_count += 1

discard_pct = discard_count/(discard_count + counted_reads)
print('number of reads discarded: ' + str(discard_pct))

#write a list of random BC then the count of that particular combination
with open(basename + '_counts.txt', 'w') as f:
    f.write('gBC\tcount\n')
    for key, value in d.items():
        if len(key) == 16:
            f.write(key + '\t' + str(value) + '\n')
