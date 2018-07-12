#!/usr/local/bin/python
# Last updated: 29/1/2018
# Author: Matt Douglas

from __future__ import print_function
import pysam
import sys
from collections import defaultdict

def sort_by_pos(introns):
    """Sort tuples of introns by chromosome, then start position, then end
    position.
    """
    return sorted(introns, key=lambda x: (x[0], int(x[1]), int(x[2])))


def regions_to_search(bam_path):
    """Count the number of lines in a BAM file."""
    regions = []

    try:
        for line in pysam.idxstats(bam_path).split('\n'):
            if len(line) > 0:
                region = line.split('\t')[0]
                if region != '*':
                    regions.append(region)
    except pysam.utils.SamtoolsError as e:
        eprint('[ERROR]', e)
        eprint('Is the BAM file indexed?')
        sys.exit(1)

    return regions


def get_introns(pos, cigar):
    """Get the pos, end, and size of any introns from the CIGAR string.
    if(  cigar_type == 0): #match
    elif(cigar_type == 1): #insertions
    elif(cigar_type == 2): #deletion
    elif(cigar_type == 3): #skip
    elif(cigar_type == 4): #soft clipping
    elif(cigar_type == 5): #hard clipping
    elif(cigar_type == 6): #padding
    """
    introns = list()

    cigar_type = [i[0] for i in cigar]
    cigar_len = [i[1] for i in cigar]

    for i in [i for i, l in enumerate(cigar_type) if l == 3]:
        size = cigar_len[i]
        start = pos
        for j in range(len(cigar_type[:i])):
            if cigar_type[j] in [0, 2, 3]:
                start += cigar_len[j]
        end = start + size - 1

        introns.append((start, end))

    return introns


def get_linked_introns(bamfile, region):
    chrom_ind = {}
    linked = set()

    for read in bamfile.fetch(region):
        cigar = read.cigartuples
        ind = read.reference_id
        chrom = bamfile.get_reference_name(read.reference_id)
        chrom_ind[ind] = chrom
        if [i[0] for i in cigar].count(3) < 2: # skip if the read is only split once
            continue
        chrom = bamfile.get_reference_name(read.reference_id)
        pos = read.pos + 1  # SAM coordinates are 1-based
        introns = get_introns(pos, cigar)
        try:
            strand = read.get_tag('XS')
        except KeyError:
            strand = '.'
        introns = tuple([(ind, i[0], i[1], strand) for i in introns])
        linked.add(introns)

    return linked, chrom_ind


def get_linked_introns_with_mate(bamfile, region):
    chrom_ind = {}
    mate_a = defaultdict(dict)
    mate_b = defaultdict(dict)
    linked = set()

    for read in bamfile.fetch(region):
        cigar = read.cigartuples
        ind = read.reference_id
        chrom = bamfile.get_reference_name(read.reference_id)
        chrom_ind[ind] = chrom
        if 3 not in [i[0] for i in cigar]: # skip if the read is not split
            continue
        name = read.query_name.rsplit('/', 1)[0] # remove the read number: "/1" or "/2"
        pos = read.pos + 1  # SAM coordinates are 1-based
        introns = get_introns(pos, cigar)
        try:
            strand = read.get_tag('XS')
        except KeyError:
            strand = '.'
        if read.is_read1:
            mate_a[(ind, strand)][name] = introns
        elif read.is_read2:
            mate_b[(ind, strand)][name] = introns

    for cs in mate_a.keys() + mate_b.keys():
        chrom, strand = cs
        for read in mate_a[cs].keys() + mate_b[cs].keys():
            if read in mate_a[cs]:
                introns_a = mate_a[cs][read]
            else:
                introns_a = []
            if read in mate_b[cs]:
                introns_b = mate_b[cs][read]
            else:
                introns_b = []
            both = set(introns_a + introns_b)
            if len(both) > 1:
                both = [(chrom, i[0], i[1], strand) for i in both]
                both = tuple(sort_by_pos(both))
                linked.add(both)

    return linked, chrom_ind


def convert_to_exons(linked):
    exons = set()

    for introns in linked:
        chrom = introns[0][0]
        strand = introns[0][3]
        for n in range(len(introns)-1):
            start = introns[n][2] + 1
            end = introns[n+1][1] - 1
            exon = chrom, start, end, strand
            exons.add(exon)
    exons = list(exons)

    return exons


def print_as_gff3(exons, path):
    with open(path, 'w') as f:
        print('##gff-version 3', file=f)
        for i in exons:
            chrom, start, end, strand = i
            line = chrom, '.', 'exon', start, end, '.', strand, '.', '.'
            print(*line, sep='\t', file=f)


def print_as_BAM(linked, header, path):
    with pysam.AlignmentFile(path, 'wb', header=header) as f:
        for n, introns in enumerate(linked):
            introns = sort_by_pos(introns)
            # calulate the postion, and distance to the next intron
            if len(introns) > 1:
                tlen = introns[-1][2] - introns[0][1] + 1
            else:
                tlen = 0
            # print out each intron as a seperate BAM entry
            for m, i in enumerate(introns):
                chrom, start, end, strand = i
                length = end - start + 1
                if m < len(introns)-1:
                    next_ref = introns[m+1][1]
                else:
                    next_ref = introns[0][1]
                    tlen = -tlen
                a = pysam.AlignedSegment()
                a.query_name           = 'linked'+str(n)
                a.query_sequence       = 'N'*length
                a.flag                 = 0
                a.reference_id         = chrom
                a.reference_start      = start
                a.mapping_quality      = 60 # 60 = unqiuely mapped for HISAT2
                a.cigartuples          = [(0,length)]
                a.next_reference_id    = chrom
                a.next_reference_start = next_ref
                a.template_length      = tlen
                a.query_qualities      = pysam.qualitystring_to_array('/'*length)
                a.tags                 = [('XN', next_ref+1), ('XI', len(introns))]
                f.write(a)


#############
# Main loop #
#############
if __name__ == '__main__':
    bam_path = sys.argv[1]
    bamfile = pysam.AlignmentFile(bam_path, 'rb') # BAM file must be sorted by read name!
    header = bamfile.header.copy()

    # create the output file for introns
    with open('linked_introns.gff3', 'w') as f:
        print('##gff-version 3', file=f)

    for region in regions_to_search(bam_path):
        print('Searching region {}...'. format(region))
        linked, chrom_ind = get_linked_introns(bamfile, region)
        print('  Found {:,} intron chains!'.format(len(linked)))
        # append the ouput to file
        with open('linked_introns.gff3', 'a') as f:
            for group in linked:
                for i in group:
                    chrom = chrom_ind[i[0]]
                    start, end, strand = i[1:]
                    line = chrom, '.', 'intron', start, end, '.', strand, '.', '.'
                    print(*line, sep='\t', file=f)
                print('#'*80, file=f)
        # print('Converting to exons...')
        # exons = convert_to_exons(linked)
        # print_as_gff3(exons, '.'.join(('linked2', region, 'gff3')))
        print('Formatting output as a BAM file...')
        print_as_BAM(linked, header, '.'.join(('linked_introns', region, 'bam')))

    # # quick and dirty, get the chromosome indeces from the BAM file
    # chrom_ind = {}
    # for region in regions_to_search(bam_path):
    #     for read in bamfile.fetch(region):
    #         ind = read.reference_id
    #         chrom = bamfile.get_reference_name(read.reference_id)
    #         chrom_ind[chrom] = ind
    #         break
    #
    # # read in from a GFF3 file that already exists
    # linked = []
    # temp = []
    # print('Reading in from GFF3...')
    # with open(sys.argv[2], 'r') as f:
    #     for line in f:
    #         if line[0] == '#':
    #             if len(temp) > 0:
    #                 linked.append(temp)
    #             temp = []
    #         else:
    #             i = line.split('\t')
    #             chrom = chrom_ind[i[0]]
    #             intron = chrom, int(i[3]), int(i[4]), i[6]
    #             temp.append(intron)
    #
    # print('Formatting output as a BAM file...')
    # print_as_BAM(linked, header, '.'.join(('linked_introns', 'bam')))
