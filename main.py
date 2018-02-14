#!/home2/mattdoug/python3/bin/python3
# Last updated: 1/15/2018
# Author: Matt Douglas

from __future__ import print_function
from get_introns_from_gff3 import get_introns
from parse_exon_index import parse_index
from reconstruct_exons_main import reconstruct_exons
import argparse, re, sys
import pysam

print('### WORK IN PROGRESS ###\n')

def pprint(*args, **kwargs):
    level = kwargs.pop('level', {'level':None})
    if level == 'debug':
        if DEBUG:
            print(*args, **kwargs)
    elif level == 'progress':
        if not QUIET and not DEBUG:
            print(*args, **kwargs)
    else:
        if not QUIET:
            print(*args, **kwargs)


def fr(region):
    """Convert a tuple of chromosome, start position, end position in a
    readable format.
    """
    return '{}:{:,}-{:,}'.format(*region[:3])


def optype(path, op='r'):
    """If the file is BAM formatted, read/write as binary."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.bam':
        op += 'b'
    return op


def parse_tab(path):
    """Parse a GFF3 file and return all intron coordinates."""
    return

    features = set()

    with open(path, 'r') as f:
        next(f)
        for line in f:
            i = line.split('\t')[0]
            i = i.split('-')
            feat = i[0], int(i[1]), int(i[2]), '+'
            features.add(feat)
            feat = i[0], int(i[1]), int(i[2]), '-'
            features.add(feat)

    return features


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                        description='Reconstruct exons from RNA-Seq data.')
    parser.add_argument('-x',
                        dest='index',
                        type=str,
                        nargs='?',
                        help='prefix of the genome index files')
    parser.add_argument('-i',
                        dest='introns',
                        type=str,
                        nargs='?',
                        help='intron in GFF3 format')
    parser.add_argument('-e',
                        dest='ret_introns',
                        type=str,
                        nargs='?',
                        help='intron retention events in GFF3 format')
    parser.add_argument('-r',
                        dest='region',
                        type=str,
                        nargs='?',
                        help='limit search to genomic coordinates (e.g. I:1000-2000)')
    parser.add_argument('-p',
                        dest='prefix',
                        type=str,
                        nargs='?',
                        default='exons', help='prefix for the output files')
    parser.add_argument('-d',
                        dest='debug',
                        action='store_true',
                        default=False,
                        help='print debug information')
    parser.add_argument('-q',
                        dest='quiet',
                        action='store_true',
                        default=False,
                        help='do not print any progress information [incompatible with -d]')

    args = parser.parse_args()
    # the required parameters are the index and the aligned reads
    if not args.index or not args.introns:
        parser.print_help()
        sys.exit(1)

    QUIET = args.quiet

    # convert a string (e.g. "I:1000..2000") to a tuple
    if args.region is not None:
        temp = args.region.strip().replace(',', '').replace('..', '-')
        temp = re.split(':|-|_', temp)
        chrom, start, end = temp[0], int(temp[1]), int(temp[2])
        args.region = chrom, start, end

    # report arguments from the command line
    pprint('Index prefix  =', args.index)
    pprint('Intron file   =', args.introns)
    pprint('Output prefix =', args.prefix)
    pprint('Search region =', [fr(args.region) if args.region is not None else 'all'][0])
    pprint('')

    index_f, index_r = parse_index(args) # positions of start & stop codons

    introns, splice_f, splice_r, site_dict = get_introns(args, min_count=2, strand_only=True)

    ret_introns = parse_tab(args.ret_introns)

    reconstruct_exons(index_f,
                      index_r,
                      introns,
                      ret_introns,
                      splice_f,
                      splice_r,
                      site_dict,
                      args)
