#!/usr/local/bin/python3
# Last updated: 2/11/2017

import argparse, re, sys
from collections import defaultdict

attr_parent   = re.compile('(?<=Parent=Transcript:)[^;]+(?=(;|))')
attr_seq_name = re.compile('(?<=sequence_name=)[^;]+(?=(;|))')
attr_biotype  = re.compile('(?<=biotype=)[^;]+(?=(;|))')

#####################
# Utility functions #
#####################
def eprint(*args, **kwargs):
    """Print to stderr."""
    if not quiet:
        print(*args, file=sys.stderr, **kwargs)


def fr(region):
    return '{}:{}-{} ({} strand)'.format(*region)


def print_as_gff3(exons, out_path):
    if not print_exons:
        return False

    with open(out_path, 'w') as f:
        print('##gff-version 3', file=f)
        for n, exon in enumerate(exons):
            chrom, left, right, strand = exon
            line   = chrom, '.', 'exon', left, right, 0, strand, '.', '.'
            print(*line, sep='\t', file=f)


#######################################
# Functions for doing the actual work #
#######################################
def parse_ref_exons(path):
    ref_pos_set   = set()
    strand_dict   = defaultdict(set)
    ref_five_end  = defaultdict(set)
    ref_three_end = defaultdict(set)

    eprint('Parsing refernce file', path)

    with open(path, 'r') as f:
        for line in f:
            if line[0] != '#':
                i = line.strip().split('\t')
                chrom  = i[0]
                kind   = i[2]
                start  = int(i[3])
                end    = int(i[4])
                strand = i[6]
                # add exon to the set
                pos = chrom, start, end, strand
                ref_pos_set.add(pos)
                # keep track of which strand the exon is on (could be both)
                s_pos = chrom, start, end
                strand_dict[s_pos].add(strand)
                # keep track of the 5' ends
                l_pos = chrom, start, strand
                ref_five_end[l_pos].add(pos)
                # keep track of the 3' ends
                r_pos = chrom, end, strand
                ref_three_end[r_pos].add(pos)

    eprint('  {:,} introns'.format(len(ref_pos_set)))
    return ref_pos_set, strand_dict, ref_five_end, ref_three_end


def parse_input_exons(path):
    exon_set = set()

    eprint('Parsing input file:', path)

    with open(path, 'r') as f:
        for line in f:
            if line[0] != '#':
                i = line.strip().split('\t')
                chrom  = i[0]
                start  = int(i[3])
                end    = int(i[4])
                strand = i[6]
                pos = chrom, start, end, strand
                exon_set.add(pos)

    eprint('  {:,} introns'.format(len(exon_set)))
    return exon_set


def check_for_exact_matches(exon_set, strand_dict, ref_match):
    not_matching = set()
    exact        = set()
    op_strand    = set()

    eprint('\r[Stage 1/3]', end='')

    for exon in exon_set:
        pos, strand = exon[:3], exon[3]
        if pos in strand_dict:              #...same position
            if strand in strand_dict[pos]:  #...and same strand
                exact.add(exon)
                ref_match.add(exon)
            else:                           #...or different strand
                op_strand.add(exon)
                ref_match.add(exon)
        else:
            not_matching.add(exon)

    print_as_gff3(exact, out_prefix+'.exact.gff3')
    print_as_gff3(op_strand, out_prefix+'.op_strand.gff3')

    eprint('\r[Stage 1/3] 100%')
    return not_matching, exact, op_strand, ref_match


def check_for_partial_matches(exon_set, ref_five_end, ref_three_end, ref_match):
    diff_five_end = set()
    diff_three_end = set()
    not_matching = set()

    eprint('\r[Stage 2/3]', end='')

    for exon in exon_set:
        chrom, start, end, strand = exon
        l_pos = chrom, start, strand
        r_pos = chrom, end, strand
        # if the left end matches...
        if l_pos in ref_five_end:
            if strand == '+':
                diff_three_end.add(exon)
            elif strand == '-':
                diff_five_end.add(exon)
            # mark that reference exon as "ref_match"
            for ref in ref_five_end[l_pos]:
                ref_match.add(ref)
        # if the right end matches...
        elif r_pos in ref_three_end:
            if strand == '+':
                diff_five_end.add(exon)
            elif strand == '-':
                diff_three_end.add(exon)
            # mark that reference exon as "ref_match"
            for ref in ref_three_end[r_pos]:
                ref_match.add(ref)
        # if neither end matches...
        else:
            not_matching.add(exon)

    print_as_gff3(diff_five_end, out_prefix+'.diff_five_end.gff3')
    print_as_gff3(diff_three_end, out_prefix+'.diff_three_end.gff3')

    eprint('\r[Stage 2/3] 100%')
    return not_matching, diff_five_end, diff_three_end, ref_match


def check_for_overlap(exon_set, ref_match):
    overlap = set()
    no_overlap = set()

    for n, exon in enumerate(exon_set):
        overlaps_any = False
        chrom, start, end, strand = exon
        # check which (if any) reference exons the exon overlaps
        for ref_exon in ref_exons:
            ref_chrom, ref_start, ref_end, ref_strand = ref_exon
            if chrom == ref_chrom and strand == ref_strand:
                if ref_start <= start <= ref_end or ref_start <= end <= ref_end:
                    overlaps_any = True
                    overlap.add(exon)
                    ref_match.add(ref_exon)
        # if the exon does not overlap any refernce exon, print it
        if not overlaps_any:
            no_overlap.add(exon)
        # report progress
        eprint('\r[Stage 3/3] {0:.2f}%'.format(100 * n / len(exon_set)), end='')

    print_as_gff3(overlap, out_prefix+'.overlaps.gff3')
    print_as_gff3(no_overlap, out_prefix+'.no_overlap.gff3')

    eprint('\r[Stage 3/3] 100%   ')
    return no_overlap, overlap, ref_match


#############
# Main loop #
#############
if __name__ == '__main__':
    qry_prefix = sys.argv[1]
    ref_prefix = '/home2/mattdoug/Thesis/reference/categorized_exons_by_type/coding'
    ref_nc_prefix = '/home2/mattdoug/Thesis/reference/categorized_exons_by_type/non_coding'
    out_prefix = 'validate'
    print_exons = True
    quiet = True

    suffix = ['internal', 'five_term', 'three_term'] # +['single']
    categories = 'total', 'exact', 'op_strand', 'diff_five_end', 'diff_three_end', 'overlap', 'no_overlap', 'total ref', 'missed'
    results = {i+n:{} for i in categories for n in ('', '_nc')}

    for s in suffix: # for each type of exon (internal, 5' terminal, etc...)
        for r, n in ((ref_prefix, ''), (ref_nc_prefix, '_nc')): # check coding, then non-coding
            qry = '.'.join((qry_prefix, s, 'gff3'))
            ref = '.'.join((r, s, 'gff3'))
            # parse each input file
            ref_exons, strand_dict, ref_five_end, ref_three_end = parse_ref_exons(ref)
            exon_set = parse_input_exons(qry)
            ref_match = set()
            results['total' + n][s] = len(exon_set)
            results['total ref' + n][s] = len(ref_exons)
            # check for exons that exactly match the reference
            not_matching, exact, op_strand, ref_match = check_for_exact_matches(exon_set, strand_dict, ref_match)
            results['exact' + n][s] = len(exact)
            results['op_strand' + n][s] = len(op_strand)
            print_as_gff3(exact, '.'.join((out_prefix, s, 'exact' + n, 'gff3')))
            print_as_gff3(op_strand, '.'.join((out_prefix, s, 'op_strand' + n, 'gff3')))
            # check for exons differing at one end
            not_matching, diff_five_end, diff_three_end, ref_match = check_for_partial_matches(not_matching, ref_five_end, ref_three_end, ref_match)
            results['diff_five_end' + n][s] = len(diff_five_end)
            results['diff_three_end' + n][s] = len(diff_three_end)
            print_as_gff3(diff_five_end, '.'.join((out_prefix, s, 'diff_five_end' + n, 'gff3')))
            print_as_gff3(diff_three_end, '.'.join((out_prefix, s, 'diff_three_end' + n, 'gff3')))
            # check for exons that differ at both ends, but overlap
            no_overlap, overlap, ref_match = check_for_overlap(not_matching, ref_match)
            results['overlap' + n][s] = len(overlap)
            results['no_overlap' + n][s] = len(no_overlap)
            print_as_gff3(overlap, '.'.join((out_prefix, s, 'overlap' + n, 'gff3')))
            print_as_gff3(no_overlap, '.'.join((out_prefix, s, 'no_overlap' + n, 'gff3')))
            # count the number of missed reference exons
            results['missed' + n][s] = len(ref_exons) - len(ref_match)

    header = [' '] + suffix
    print(*header, sep='\t')

    for b in ('', '_nc'):
        for a in categories:
            category = a + b
            line = [category]
            for s in suffix:
                val = results[category][s]
                line.append(val)
            print(*line, sep='\t')
