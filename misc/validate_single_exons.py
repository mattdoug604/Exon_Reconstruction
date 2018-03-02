#!/usr/local/bin/python
#Author: Matt Douglas
#Last updated: 19/11/2017

from __future__ import print_function, division
import argparse, pysam, re, sys
from collections import defaultdict


class Exon(object):
    def __init__(self, chrom, start, end, strand):
        self.chrom  = chrom
        self.start  = start
        self.end    = end
        self.strand = strand
        self.len    = end - start + 1
        self.cov    = []

    def __str__(self):
        return fr( (self.chrom, self.start, self.end, self.strand) )

    def print_as_gff3(self, n):
        return '\t'.join(map(str, (self.chrom, '.', 'CDS', self.start, self.end, sum(self.cov), self.strand, '.', 'ID=exon{}'.format(n)) ))


#####################
# Utility functions #
#####################
def eprint(*args, **kwargs):
    """Print to stderr if 'quiet' is False."""
    if not quiet:
        print(*args, file=sys.stderr, **kwargs)


def dprint(*args, **kwargs):
    """Print to stderr if 'debug' is True."""
    if debug:
        print(*args, file=sys.stderr, **kwargs)


def fr(region):
    """Return a tuple of (chromosome, start position, end position, strand) in
    a readable format.
    """
    if len(region) < 4:
        strand = ''
    else:
        strand = ' ({} strand)'.format(region[3])

    return '{}:{:,}-{:,}{}'.format(region[0], region[1], region[2], strand)


def prog(count, total):
    progress = count * 100 / total
    progress = '{0:.2f}%'.format(progress)

    return progress


def parse_region(region_coord):
    """Parse a set of genomic coordinates into a tuple of (chromosome, start
    position, end position).
    """
    region_coord = region_coord.strip().replace(',', '').replace('..', '-')
    chrom, left, right = re.split(':|-|_', region_coord)
    return chrom, int(left), int(right)


def sort_exons(exons):
    """Sort tuples of exons by chromosome, then start position, then end
    position. NOTE: C. elegans uses roman numerals for chromosomes names.
    """
    numerals = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'X':10, 'MtDNA':11}
    try:
        return sorted(exons, key=lambda x: (numerals[x.chrom], x.start, x.end))
    except KeyError:
        return sorted(exons, key=lambda x: (x.chrom, x.start, x.end))


def print_as_gff3(exons, out_path, kind='CDS'):
    """Convert the results from a tuple (chromosome, start, end, strand to GFF3
    format and print.
    """
    dprint('Printing results to "{}"\n'.format(out_path))

    with open(out_path, 'w') as f:
        print('##gff-version 3', file=f)
        for n, exon in enumerate(exons):
            chrom, left, right, strand = exon
            count  = 0
            line   = [chrom, '.', kind, left, right, count, strand, '.', 'ID=exon{}'.format(n+1)]
            print(*line, sep='\t', file=f)


################################################
# Parse command line arguments and input files #
################################################
def parse_arguments():
    """Parse arguments given on the command line."""
    parser = argparse.ArgumentParser(description='Reconstruct exons from a set of introns and mapped reads.')
    parser.add_argument('-e', type=str, nargs='?', help='path to a GFF3 file containing single-exons')
    parser.add_argument('-a', type=str, nargs='?', help='BAM file of aligned reads')
    parser.add_argument('-t', type=float, nargs='?', help='minimum depth of coverage')
    parser.add_argument('-p', type=str, nargs='?', help='prefix for the output file')
    parser.add_argument('-l', type=str, nargs='?', help='limit search to the specified genomic coordinates (e.g. I:1000-2000)')
    parser.add_argument('-q', action='store_true', help='do not print any progress information [incompatible with -d]')
    parser.add_argument('-d', action='store_true', help='print debug information')
    args = parser.parse_args()

    if args.e is None or args.a is None:
        parser.print_help()
        sys.exit(1)

    if args.t is None:
        threshold = 1.0
    else:
        threshold = args.t

    if args.p is None:
        out_path = 'exons'
    else:
        out_path = args.p.rstrip('.gff3')

    if args.l is not None:
        region = parse_region(args.l)
    else:
        region = None

    if args.d:
        args.q = True

    return args.e, args.a, threshold, out_path, region, args.q, args.d


def parse_exons(exon_path, region=None):
    exons = []

    eprint('Parsing single-exons...', end='')
    with open(exon_path, 'r') as f:
        for line in f:
            if line[0] != '#':
                i = line.split('\t')
                chrom, left, right, strand = i[0], int(i[3]), int(i[4]), i[6]
                # TEMP: some introns have len 0 (need to fix)
                if right - left + 1 < 1:
                    continue
                # if the feature is outside the specified region, skip it
                if region is not None:
                    if chrom != region[0]:
                        continue
                    if left > region[2] or right < region[1]:
                        continue
                exon = Exon(chrom, left, right, strand)
                exons.append(exon)

    # sort exons by position (chrom, start, end)
    exons = sort_exons(exons)

    eprint('\rParsing single-exons... Done!')
    return exons


def compute_cov(exon, samfile):
    cov_list = []

    for col in samfile.pileup(exon.chrom, exon.start, exon.end):
        if exon.start <= col.pos <= exon.end:
            cov_list.append(col.n)

    return cov_list


#############
# Main loop #
#############
if __name__ == '__main__':
    # parse command line arguements
    exon_path, bam_path, threshold, out_path, region, quiet, debug = parse_arguments()
    samfile = pysam.AlignmentFile(bam_path, 'rb')

    # report arguments from the command line to stderr
    eprint('Exon file     =', exon_path)
    eprint('BAM file      =', bam_path)
    eprint('Minimum depth =', threshold)
    eprint('Output prefix =', out_path)
    eprint('Search region =', [fr(region) if region is not None else 'all'][0])
    eprint(' ')

    # identify all exons and coverage islands
    exons = parse_exons(exon_path, region)

    # compute the coverage for each exon
    total = len(exons)
    for n, e in enumerate(exons):
        eprint('\rComputing exon coverage...', prog(n, total), end='')
        cov = compute_cov(e, samfile)
        if cov > 0:
            print(e.print_as_gff3(n+1))
    eprint('\rComputing exon coverage... Done!')
