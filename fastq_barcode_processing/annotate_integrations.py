#!/usr/bin/python

import argparse
import pysam
from collections import Counter
from collections import defaultdict

parser = argparse.ArgumentParser(
    description = 'Read bam file and output barcode with annotated locations')
parser.add_argument('-i', nargs = '*', help = 'bamfile')
parser.add_argument('-o', help = 'output file name')
parser.add_argument('-u', action = 'store_true', help = 'whether to output unique locations only')
args = parser.parse_args()

def rev_comp(seq):
    
    new_seq = ''
    conversion = {'A': 'T', 'T': 'A', 'G':'C', 'C':'G'}
    
    for i in reversed(seq):
        new_seq += conversion[i]
        
    return new_seq

def check_Ns(seq):

    valid = ['A', 'T', 'C', 'G']

    for i in seq:
        if i not in valid:
            return False
    
    return True

def get_strand(read):
    
    if read.is_reverse:
        return '-'
    
    return '+'

bamfiles = {}

for i in args.i:
    sample = i.strip('.bam')
    bamfiles[sample] = pysam.AlignmentFile(i, "rb")

unmapped_reads = open(args.o + '_unmapped', 'w')
read_depths = {}

def read_bamfile(bamfile):
    current_locations = defaultdict(list)
    read_depth = 0
    for read in bamfile:
        read_depth += 1
        if not read.is_unmapped:
            trip_bc = read.query_name.split(':')[-1]
            chromosome = read.reference_name
            location = str(read.get_reference_positions()[0])
            if check_Ns(trip_bc) == True:
                strand = get_strand(read)
                current_locations[rev_comp(trip_bc)].append((chromosome, location, strand))
        else:
            unmapped_reads.write(read.get_forward_sequence() + '\n')

    return current_locations, read_depth

def collate_barcodes(current_locations):
    barcode_locations = {}
    for bcs, locations in current_locations.items():
        count_items = Counter(locations)
        top = count_items.most_common(2)
        total = sum(count_items.values())
        if top[0][1]/total > 0.8:
            barcode_locations[bcs] = [top[0][0], top[0][1]]
        elif (top[0][1] + top[1][1])/total > 0.8:
            loc1 = top[0][0]
            loc2 = top[1][0]
            if loc1[0] == loc2[0]:
                if abs(int(loc1[1]) - int(loc2[1])) < 1000:
                    barcode_locations[bcs] = [top[0][0], top[0][1] + top[1][1]]

    return barcode_locations

all_barcode_locations = {}
all_read_depths = {}

for sample, bamfile in bamfiles.items():
    current_locations, current_read_depth = read_bamfile(bamfile)
    all_read_depths[sample] = current_read_depth
    current_locations = {bcs: locations for bcs, locations in current_locations.items() if len(locations) > 2}
    collated_locations = collate_barcodes(current_locations)
    all_barcode_locations[sample] = collated_locations

def compare_locs(stored_loc, current_loc, current_sample):

    stored_loc_norm_count = stored_loc[1]/all_read_depths[stored_loc[2]]
    current_loc_norm_count = current_loc[1]/all_read_depths[current_sample]

    if stored_loc_norm_count/current_loc_norm_count > 0.8:
        return stored_loc
    else:
        return current_loc + [sample]

combined_barcode_locations = {}

for sample, bc_locations in all_barcode_locations.items():
    for bcs, location in bc_locations.items():
        try:
            stored_loc = combined_barcode_locations[bcs]
            if stored_loc[0] == location[0]:
                pass
            else:
                combined_barcode_locations[bcs] = compare_locs(stored_loc, location, sample)
        except:
            combined_barcode_locations[bcs] = location + [sample]
            
if args.u == True:
    all_filtered_locations = {}

    for bcs, location in combined_barcode_locations.items():
        current_location = location[0][:3]
        if current_location in all_filtered_locations.keys():
            current_count = location[1]
            stored_count = all_filtered_locations[current_location][1]
            if current_count > stored_count:
                all_filtered_locations[current_location] = (bcs, location[1])
        else:
            all_filtered_locations[current_location] = (bcs, location[1])

    with open(args.o + '_unique.bed', 'w') as f:
        for location, bc_counts in all_filtered_locations.items():
            f.write('chr{chrom}\t{loc}\t{loc}\t{tBC}\t{strand}\n'.format(
                chrom = location[0], loc = location[1], tBC = bc_counts[0], strand = location[2]))

else:
    with open(args.o + '.bed', 'w') as f:
        for bcs, location in combined_barcode_locations.items():
            f.write('chr{chrom}\t{loc}\t{loc}\t{tBC}\t{count}\t{strand}\n'.format(
                chrom = location[0][0], loc = location[0][1], tBC = bcs, count = location[1], strand = location[0][2]))

bamfile.close()
unmapped_reads.close()
