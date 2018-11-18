#!/home2/mattdoug/python3/bin/python3
# Last updated: 25/4/2018
# Author: Matt Douglas

# Purpose: Get the coordinates of each intron from a GFF3 file, and index the
# positions of each individual splice site.

import sys
import pysam
import logging
from collections import defaultdict

log_format = '%(message)s'
logging.basicConfig(format=log_format, level=logging.INFO)
log = logging.getLogger(__name__)


def parse_gff3(path, region=None):
    intron_dict = defaultdict(int)
    f_count = 0
    r_count = 0
    discard = 0

    log.info('\nGetting intron coordinates...')
    with open(path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            line = line.split('\t')
            intron = line[0], int(line[3]), int(line[4]), line[6]
            strand = line[6]
            if strand == '+':
                f_count += 1
            elif strand == '-':
                r_count += 1
            else:
                discard += 1
                continue
            # make sure the score meets the minimum
            if line[5] != '.':
                count = int(line[5])
            else:
                count = 0.0
            # if the feature is outside the specified region, skip it
            if region is not None:
                if intron[0] != region[0]:
                    continue
                if region[1] <= intron[1] <= region[2]:
                    intron_dict[intron] = count
                elif region[1] <= intron[2] <= region[2]:
                    intron_dict[intron] = count
            else:
                intron_dict[intron] = count

    log.info('  Found {:,} valid introns:'.format(f_count + r_count))
    log.info('    {:,} on the + strand'.format(f_count))
    log.info('    {:,} on the - strand'.format(r_count))
    if discard > 0: log.info('  Discarded {:,} introns without a defined strand'.format(discard))

    return intron_dict


def parse_CIGAR(pos, cigar):
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

        introns.append([start, end])

    return introns


def get_introns(args):
    global QUIET
    global DEBUG
    QUIET = args.quiet
    DEBUG = args.debug
    PREFIX = args.output
    gff_path = args.introns
    region = args.region
    ref_set = None
    intron_dict = defaultdict(int)
    splice_f = {}
    splice_r = {}
    site_l = defaultdict(list)
    site_r = defaultdict(list)

    intron_dict = parse_gff3(gff_path, region)

    # build dictionaries tracking: A) the positions of each splice site on the
    # forward and reverse strands, and B) which introns share a splice site.
    log.debug('Indexing splice sites...')
    for intron, count in intron_dict.items():
        chrom, left, right, strand = intron
        if strand == '+':
            for n in range(3):
                try:
                    splice_f[chrom][n].append((left, 1))
                    splice_f[chrom][n].append((right+1, 2))
                except KeyError:
                    splice_f[chrom] = [ [], [], [] ]
                    splice_f[chrom][n].append((left, 1))
                    splice_f[chrom][n].append((right+1, 2))
        elif strand == '-':
            for n in range(3):
                try:
                    splice_r[chrom][n].append((left-1, 1))
                    splice_r[chrom][n].append((right, 2))
                except KeyError:
                    splice_r[chrom] = [ [], [], [] ]
                    splice_r[chrom][n].append((left-1, 1))
                    splice_r[chrom][n].append((right, 2))
        # keep track of introns sharing a splice site
        try:
            site_l[(chrom, left, strand)].append(intron)
        except KeyError:
            site_l[(chrom, left, strand)] = [intron]
        try:
            site_r[(chrom, right, strand)].append(intron)
        except KeyError:
            site_r[(chrom, right, strand)] = [intron]
    log.debug('Indexing splice sites... Done!')

    return intron_dict, splice_f, splice_r, site_l, site_r
