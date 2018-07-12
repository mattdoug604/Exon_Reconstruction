#!/home2/mattdoug/python3/bin/python3
# Author: Matt Douglas
# Last updated: 2/6/2018

# PURPOSE: Use the translation blocks encoded in the genome and splice sites
#          identified from RNA-Seq data to reconstruct protein-coding exons.
# USAGE: reconstruct_exons.py [-q] [-d] [-p exons] [-l I:100..200] -x <prefix of index files> -i <path to introns>

# UPDATES TO THIS VERSION:
# Changed the script to take a list of known intron retention events and use
# those to reconstruct exons. Intron retention events can be detected using
# iREAD (Hong-Dong et al. 2017)

# NOTES:
# I:217,489..219,488 - example of resolved 2-exon gene
# X:2,412,849..2,413,848 - example of incorrect 2-exon gene
# III:3,562,617..3,582,616 - PacBio reads support alt intron w/ 2 support
# I:5,564,625..5,565,624 - missing internal exon due to phase incompatibility?
# X:13,397,323..13,397,775 - gene with 2 introns, but only 1 is part of coding exons
# II:1,889,845..1,890,844 - what's going on here?
# I:506,491..508,490 - False intron retention events?
# II:8,042,500..8,043,100 - intron retention in 5' terminal exon

from collections import defaultdict
from itertools import chain, product
import argparse, re, sys
import pysam

#####################
# Class definitions #
#####################
class TranslationBlock(object):
    __slots__ = ['start', 'end', 'l_sites', 'r_sites', 's_sites']
    def __init__(self):
        self.start = -1
        self.end   = -1
        self.l_sites = [] # 5' splice sites
        self.r_sites = [] # 3' splice sites
        self.s_sites = [] # start codons

    def __len__(self):
        return self.end - self.start + 1

    def num_sites(self):
        return sum([len(i) for i in (self.l_sites, self.r_sites, self.s_sites)])


#####################
# Utility functions #
#####################
def pprint(*args, **kwargs):
    level = kwargs.pop('level', {'level':None})
    if level == 'debug':
        if DEBUG:
            print(*args, file=sys.stderr, **kwargs)
    elif level == 'progress':
        if not QUIET and not DEBUG:
            print(*args, file=sys.stderr, **kwargs)
    else:
        if not QUIET:
            print(*args, file=sys.stderr, **kwargs)


def fr(region):
    """Return a tuple of (chromosome, start position, end position, strand) in
    a more readable format.
    """
    return '{}:{}-{}{}'.format(*region)


def percent(count, total):
    if total == 0:
        return '0%'
    per = count * 100 / total
    per = '{0:.1f}%'.format(per)

    return per


def merge_dicts(x, y):
    z = {}

    for chrom, frames in x.items():
        if chrom not in z:
            z[chrom] = [ [], [], [] ]
        for n, frame in enumerate(frames):
            z[chrom][n] += frame
    for chrom, frames in y.items():
        if chrom not in z:
            z[chrom] = [ [], [], [] ]
        for n, frame in enumerate(frames):
            z[chrom][n] += frame

    return z


def shift_frame(adj_frame, intron, direction):
    """Calculate the frame a given exon uses based on the frame of the adjacent
    exon and the lenght of the adjacent intron.
    """
    intron_len = intron[2]-intron[1]+1
    modulo = intron_len % 3

    if direction in ('+<', '->'):
        new_frame = (adj_frame - modulo) % 3
    elif direction in ('+>', '-<'):
        new_frame = (adj_frame + modulo) % 3

    pprint('    Calculating frame shift:', level='debug')
    pprint('      Direction is {}'.format(direction), level='debug')
    pprint('      Adjacent exon is in frame {}'.format(adj_frame), level='debug')
    pprint('      Intron {} is {}bp long, modulo is {}'.format(fr(intron), intron_len, modulo), level='debug')
    pprint('      New frame is {}'.format(new_frame), level='debug')

    return new_frame


def sort_indeces(index_dict):
    """Take a dictionary of tuples, remove duplicates and sort by position."""
    pprint('\rRemoving duplicates and sorting by position...', end='', level='debug')
    for chrom, frames in index_dict.items():
        for n, f in enumerate(frames):
            index_dict[chrom][n] = sorted(set(f))
    pprint('\rRemoving duplicates and sorting by position... Done!', level='debug')

    return index_dict


def sort_by_pos(features):
    """Sort tuples of exons by chromosome, then start position, then end
    position. NOTE: C. elegans uses roman numerals for chromosomes names.
    """
    numerals = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'X':10, 'MtDNA':11}
    try:
        return sorted(features, key=lambda x: (numerals[x[0]], int(x[1]), int(x[2])))
    except KeyError:
        return sorted(features, key=lambda x: (x[0], int(x[1]), int(x[2])))


def print_as_gff3(features, out_path, kind='CDS'):
    """Convert the results from a tuple (chromosome, start, end, strand to GFF3
    format and print.
    """
    pprint('Printing results to "{}"'.format(out_path), level='debug')

    with open(out_path, 'w') as f:
        print('##gff-version 3', file=f)
        for n, feat in enumerate(features):
            chrom, left, right, strand = feat
            score = '.'
            phase = '.'
            frame = ','.join(map(str, FRAME[feat]))
            notes = [';frame='+frame if frame != '' else ''][0]
            line = chrom, '.', kind, left, right, score, strand, phase, 'ID={}{}{}'.format(kind, n+1, notes)
            print(*line, sep='\t', file=f)


##################
# Main functions #
##################
def parse_csv(path):
    ref_dict = {}

    with open(path, 'r') as f:
        for line in f:
            i = line.strip().split('\t')
            left, right, lphase, rphase, lstat, rstat, intron = i
            # format left
            temp, strand = left[:-1], left[-1]
            temp = re.split(':|-|_', temp)
            temp[1], temp[2] = int(temp[1]), int(temp[2])
            left = tuple(temp + [strand])
            # format right
            temp, strand = right[:-1], right[-1]
            temp = re.split(':|-|_', temp)
            temp[1], temp[2] = int(temp[1]), int(temp[2])
            right = tuple(temp + [strand])
            # format intron
            temp, strand = intron[:-1], intron[-1]
            temp = re.split(':|-|_', temp)
            temp[1], temp[2] = int(temp[1]), int(temp[2])
            intron = tuple(temp + [strand])
            # output
            ref_dict[intron] = [left, right]

    return ref_dict


def parse_index(path, region=None):
    index_f = {}
    index_r = {}
    start_c = 0
    end_c = 0

    pprint('Parsing reference index...', end='')
    for i in range(6):
        index = path + '.' + str(i+1)
        with open(index, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if len(line) < 3: #if no start/stop codons were found, the 3rd column may be empty
                    continue
                chrom, kind, pos_list = line
                chrom = str(chrom)
                kind  = int(kind)
                # if the feature is outside the specified region, skip it
                if region is not None and chrom != region[0]:
                    continue
                pos_list = list(map(int, pos_list.split(',')))
                if region is not None:
                    temp = []
                    for n in pos_list:
                        if region[1] <= n <= region[2]:
                            temp.append(n)
                    pos_list = temp
                # add the features to the dictionaries
                d = [index_f if i < 3 else index_r][0] # 1 == "+", -1 == "-"
                frame = [i-3 if i >= 3 else i][0]
                try:
                    d[chrom][frame] += [(x, kind) for x in pos_list]
                except KeyError:
                    d[chrom] = [ [], [], [] ]
                    d[chrom][frame] += [(x, kind) for x in pos_list]
                # update the counters
                if kind == 0:
                    start_c += len(pos_list)
                elif kind == 3:
                    end_c += len(pos_list)
    pprint('\rParsing reference index... Done!')
    pprint('  {:,} putative start codons'.format(start_c))
    pprint('  {:,} stop codons'.format(end_c))

    return index_f, index_r


def parse_gff3(path, region=None):
    intron_set = set()
    splice_f = {}
    splice_r = {}
    site_dict = {}

    pprint('Parsing specified intron file...', end='')
    with open(path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            line = line.split('\t')
            intron = line[0], int(line[3]), int(line[4]), line[6]
            # if the feature is outside the specified region, skip it
            if region is not None:
                if region[0] == intron[0]:
                    if region[1] <= intron[1] <= region[2]:
                        intron_set.add(intron)
                    elif region[1] <= intron[2] <= region[2]:
                        intron_set.add(intron)
            else:
                intron_set.add(intron)
    pprint('\rParsing specified intron file... Done!')

    pprint('Indexing splice sites...', end='', level='debug')
    f_count, r_count = 0, 0
    for intron in intron_set:
        chrom, left, right, strand = intron
        if strand == '+':
            f_count += 1
            for n in range(3):
                try:
                    splice_f[chrom][n].append((left, 1))
                    splice_f[chrom][n].append((right+1, 2))
                except KeyError:
                    splice_f[chrom] = [ [], [], [] ]
                    splice_f[chrom][n].append((left, 1))
                    splice_f[chrom][n].append((right+1, 2))
        elif strand == '-':
            r_count += 1
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
            site_dict[(chrom, left, strand)].append(intron)
        except KeyError:
            site_dict[(chrom, left, strand)] = [intron]
        try:
            site_dict[(chrom, right, strand)].append(intron)
        except KeyError:
            site_dict[(chrom, right, strand)] = [intron]
    pprint('\rIndexing splice sites... Done!', level='debug')

    pprint('Reporting {:,} unique introns:'.format(f_count + r_count))
    pprint('  {:,} on the + strand'.format(f_count))
    pprint('  {:,} on the - strand'.format(r_count))

    return intron_set, splice_f, splice_r, site_dict


def get_translation_blocks(index_f, index_r):
    """Return dictionaries of TranslationBlock objects. A translation block is
    defined as the region from a stop codon to the end of the next stop codon.

    INPUT: Two dictionaries corresponding to the forward and reverse strands of
    the genome. Keys in each dict are the chromosomes. The values are 3 lists
    corresponding to each reading frame:
    index_f = { chromsome: [ [frame0], [frame1], [frame2] ] }
    index_r = { chromsome: [ [frame0], [frame1], [frame2] ] }

    RETURN: Two dictionaries (blocks_f, blocks_r), structured the same as the
    input. Each frame is a list of TranslationBlock objects, which contain the
    positions of start codons and splice sites within them.
    """
    sort_by_type = lambda x: sorted(x, key=lambda y: (y[0], y[1]))
    blocks_f = {}
    blocks_r = {}
    empty_blocks_f = {}
    empty_blocks_r = {}
    num_blks = 0
    count = 0
    total = len(index_f)*6 #calculate the number of positions to check

    pprint('', level='debug')
    pprint('Identifying all translation blocks', end='')

    for chrom, frames in index_f.items():
        blocks_f[chrom] = [ [], [], [] ]
        empty_blocks_f[chrom] = [ [], [], [] ]
        for frame, i in enumerate(frames):
            pprint('\rIdentifying all translation blocks...', percent(count, total), end='', level='progress')
            pprint('\nScanning + strand in the > direction, reading frame {}'.format(frame), level='debug')
            if len(i) < 1: # skip if there are no start codons or splice sites
                continue
            block = TranslationBlock()
            last = i[0][0]
            count += 1
            for pos, kind in sort_by_type(i[1:]):
                if kind == 0: # if start codon
                    block.s_sites.append(pos)
                elif kind == 1: # left splice site
                    block.r_sites.append(pos)
                elif kind == 2: # right splice site
                    block.l_sites.append(pos)
                elif kind == 3: # if stop codons
                    block.start = last + 1
                    block.end = pos
                    if block.num_sites() >= 0:
                        blocks_f[chrom][frame].append(block)
                    else:
                        empty_blocks_f[chrom][frame].append(block)
                    # DEBUG
                    pprint('  {}:{:,}-{:,} (+ strand)'.format(chrom, block.start, block.end), level='debug')
                    pprint("    5' splicesites:", ' '.join(map(str, block.l_sites)), level='debug')
                    pprint("    3' splicesites:", ' '.join(map(str, block.r_sites)), level='debug')
                    pprint("    start sites:", ' '.join(map(str, block.s_sites)), level='debug')
                    # reset
                    num_blks += 1
                    last = pos
                    block = TranslationBlock()
    ############################################################################
    for chrom, frames in index_r.items():
        blocks_r[chrom] = [ [], [], [] ]
        empty_blocks_r[chrom] = [ [], [], [] ]
        for frame, i in enumerate(frames):
            pprint('\rIdentifying all translation blocks...', percent(count, total), end='', level='progress')
            pprint('\nScanning - strand in the > direction, reading frame {}'.format(frame), level='debug')
            if len(i) < 1: # skip if there are no start codons or splice sites
                continue
            block = TranslationBlock()
            last = i[0][0]-1
            count += 1
            for pos, kind in sort_by_type(i[1:]):
                if kind == 0: # start codon
                    block.s_sites.append(pos)
                elif kind == 1: # left splice site
                    block.r_sites.append(pos)
                elif kind == 2: # right splice site
                    block.l_sites.append(pos)
                elif kind == 3: # stop codon
                    block.start = last + 1
                    block.end = pos
                    if block.num_sites() >= 0:
                        blocks_r[chrom][frame].append(block)
                    else:
                        empty_blocks_r[chrom][frame].append(block)
                    # DEBUG
                    pprint('  {}:{:,}-{:,} (- strand)'.format(chrom, block.start, block.end), level='debug')
                    pprint("    5' splicesites:", ' '.join(map(str, block.l_sites)), level='debug')
                    pprint("    3' splicesites:", ' '.join(map(str, block.r_sites)), level='debug')
                    pprint("    start sites:", ' '.join(map(str, block.s_sites)), level='debug')
                    # reset
                    num_blks += 1
                    last = pos
                    block = TranslationBlock()

    pprint('\rIdentifying all translation blocks... Done!    ')
    pprint('  {:,} translation blocks total'.format(num_blks))

    return blocks_f, blocks_r


def get_putative_exons(blocks_f, blocks_r):
    """Search each translation block and return the coordinates of all possible
    internal and terminal exons. Also note which reading frame(s) each uses and
    which intron(s) the exons are adjacent to (for terminal exons only).

    INPUT: Two dictionaries corresponding to the translation block objects
    present in each reading frame on forward and reverse strands of the genome.
    Keys for each dict are the chromosomes. The values are 3 lists corresponding
    to each reading frame:
    blocks_f = { chromsome: [ [frame0], [frame1], [frame2] ] }
    blocks_r = { chromsome: [ [frame0], [frame1], [frame2] ] }

    RETURN:
    1) a set of "internal exons" as tuples in the format: (chromosome, start,
       end, strand).
    2) two dictionaries of terminal exons where the value is the exon position
       (chromosome, start, end, strand) and the key is the adjacent intron,
    """
    exon_l_term_dict = defaultdict(set)
    exon_r_term_dict = defaultdict(set)
    five_term_c = 0
    three_term_c = 0
    count = 0
    total = len(blocks_f)*6  #calculate the number of positions to check

    pprint('', level='debug')
    pprint('Identifying putative exons positions', end='')

    # scan the + strand of the genome
    for chrom, frames in blocks_f.items():
        for frame, i in enumerate(frames):
            pprint('\rIdentifying putative exon positions...', percent(count, total), end='', level='progress')
            pprint('\nScanning + strand, reading frame {}'.format(frame), level='debug')
            count += 1
            for block in i:
                if block.num_sites() < 1:
                    continue
                l_sites = sorted(block.l_sites)
                r_sites = sorted(block.r_sites)
                s_sites = sorted(block.s_sites)
                pprint('  Checking block: {}:{:,}-{:,} (+ strand)'.format(chrom, block.start, block.end), level='debug')
                pprint('  l_sites:', l_sites, level='debug')
                pprint('  r_sites:', r_sites, level='debug')
                pprint('  s_sites:', s_sites, level='debug')
                # join start codons to 5' splice sites
                for s in s_sites:
                    x = [r for r in r_sites if r > s]
                    for i, j in product([s], x):
                        exon = chrom, i, j-1, '+'
                        r_adj = chrom, j-1, '+'
                        exon_l_term_dict[r_adj].add(exon)
                        FRAME[exon].add(frame)
                        five_term_c += 1
                        pprint("    5' terminal exon found at {}".format(fr(exon)), level='debug')
                # join 3' splice sites to stop codons
                for l in [s for s in l_sites if s < block.end]:
                    exon = chrom, l, block.end, '+'
                    l_adj = chrom, l, '+'
                    exon_r_term_dict[l_adj].add(exon)
                    FRAME[exon].add(frame)
                    three_term_c += 1
                    pprint("    3' terminal exon found at {}".format(fr(exon)), level='debug')
    ############################################################################
    # scan the - strand of the genome
    for chrom, frames in blocks_r.items():
        for frame, i in enumerate(frames):
            pprint('\rIdentifying putative exon positions...', percent(count, total), end='', level='progress')
            pprint('\nScanning - strand, reading frame {}'.format(frame), level='debug')
            count += 1
            for block in i:
                if block.num_sites() < 1:
                    continue
                LCOUNT = 0
                l_sites = sorted(block.l_sites)
                r_sites = sorted(block.r_sites)
                s_sites = sorted(block.s_sites)
                pprint('  Checking block: {}:{:,}-{:,} (- strand)'.format(chrom, block.start, block.end), level='debug')
                pprint('  l_sites:', l_sites, level='debug')
                pprint('  r_sites:', r_sites, level='debug')
                pprint('  s_sites:', s_sites, level='debug')
                # join stop codon to 3' splice sites
                for r in [s for s in r_sites if s > block.start]:
                    exon = chrom, block.start, r, '-'
                    r_adj = chrom, r, '-'
                    exon_l_term_dict[r_adj].add(exon)
                    FRAME[exon].add(frame)
                    three_term_c += 1
                    pprint("    5' terminal exon found at {}".format(fr(exon)), level='debug')
                # join 5' splice sites to start codons
                for l in l_sites:
                    x = [s for s in s_sites if s > l]
                    for i, j in product([l], x):
                        exon = chrom, i+1, j, '-'
                        l_adj = chrom, i+1, '-'
                        exon_r_term_dict[l_adj].add(exon)
                        FRAME[exon].add(frame)
                        five_term_c += 1
                        pprint("    3' terminal exon found at {}".format(fr(exon)), level='debug')

    pprint('\rIdentifying putative exon positions... Done!    ')
    pprint("  {:,} putative 5' terminal exons".format(five_term_c))
    pprint("  {:,} putative 3' terminal exons".format(three_term_c))

    return exon_l_term_dict, exon_r_term_dict


def terminal_exons_in_phase(exon_l_term_dict, exon_r_term_dict):
    # quick function to calculate the length of an exon
    length = lambda x: x[2] - x[1] + 1

    all_dict = {}
    final_dict = {}
    five_term_c = 0
    three_term_c = 0
    exon_five_term = []
    exon_three_term = []
    # for keeping track of progress
    count = 0
    total = len(ADJ)
    possible = 0

    pprint('', level='debug')
    pprint('Checking which terminal exons are in phase', end='')

    for intron, adj in ADJ.items():
        pprint('\rChecking which terminal exons are in phase...', percent(count, total), end='', level='progress')
        count += 1
        ########################################################################
        chrom, l, r, strand = intron
        l_adj, r_adj = adj
        all_dict[intron] = [[], []]
        final_dict[intron] = [[], []]
        ########################################################################
        # if no internal exons on both sides of the intron...
        if len(l_adj) < 1 and len(r_adj) < 1:
            pprint('Looking on both sides of intron at {}'.format(fr(intron)), level='debug')
            results = []
            l_compare = {}
            r_compare = {}
            # which exons are in which frames, on both sides
            if len(exon_l_term_dict[(chrom, l-1, strand)]) * len(exon_r_term_dict[(chrom, r+1, strand)]) == 0:
                pprint('  No blocks found on one or both sides...', level='debug')
                continue
            for exon in exon_l_term_dict[(chrom, l-1, strand)]:
                all_dict[intron][0].append(exon)
                for f in FRAME[exon]:
                    if f not in l_compare:
                        l_compare[f] = set()
                    l_compare[f].add(exon)
            for exon in exon_r_term_dict[(chrom, r+1, strand)]:
                all_dict[intron][1].append(exon)
                for f in FRAME[exon]:
                    if f not in r_compare:
                        r_compare[f] = set()
                    r_compare[f].add(exon)
            pprint('  Left side:', level='debug')
            for frame, exons in l_compare.items():
                pprint('    Frame {}) {}'.format(frame, ' '.join([fr(i) for i in exons])), level='debug')
            pprint('  Right side:', level='debug')
            for frame, exons in r_compare.items():
                pprint('    Frame {}) {}'.format(frame, ' '.join([fr(i) for i in exons])), level='debug')
            # check which exons are in frame with each other
            # check if they are a valid ORF (i.e. modulo 3 == 0)
            for f, l_set in l_compare.items():
                new_f = shift_frame(f, intron, strand+'>')
                if new_f in r_compare:
                    r_set = r_compare[new_f]
                else:
                    continue
                for l_exon in l_set:
                    for r_exon in r_set:
                        if ( length(l_exon) + length(r_exon)) % 3 == 0:
                            results.append((l_exon, r_exon))
            # if multiple are valid, select the longest ORF
            possible += len(results)
            if len(results) > 0:
                final = max(results, key=lambda x: (x[0][2]-x[0][1] + x[1][2]-x[1][1]))
                final_dict[intron][0] = final[0]
                final_dict[intron][1] = final[1]
                five_term_c += 1
                three_term_c += 1
                pprint('  Longest ORF is:', level='debug')
                pprint('    {} is left of intron {}'.format(fr(final[0]), fr(intron)), level='debug')
                pprint('    {} is right of intron {}'.format(fr(final[1]), fr(intron)), level='debug')
            else:
                pprint('  No valid ORF found!', level='debug')
        ########################################################################
        # otherwise, continue...
        else:
            pprint('Skipping intron at {}, is flanked by exons'.format(fr(intron)), level='debug')
            continue

    pprint('\rChecking which terminal exons are in phase... Done!  ')
    pprint("  {:,} 5' terminal exons".format(five_term_c))
    pprint("  {:,} 3' terminal exons".format(three_term_c))
    pprint('{:,} possible pairings'.format(possible))

    return all_dict, final_dict


############
# Run loop #
############
if __name__ == '__main__':
    # quick function to calculate the length of an exon
    length = lambda x: x[2] - x[1] + 1

    index_path = '/home2/mattdoug/Thesis/reference/index_exons_reconstruction/exons.index'
    stats_path = sys.argv[1]
    intron_path = sys.argv[2]
    # region = 'III', 11608000, 11612999
    region = None
    QUIET = False
    DEBUG = False
    FRAME = defaultdict(set)
    NOTE = defaultdict(set)

    # parse the csv file
    ref_dict = parse_csv(stats_path)

    # Build the indeces
    introns, splice_f, splice_r, SITE = parse_gff3(intron_path, region)
    ADJ = {i:[set(), set()] for i in introns}

    index_f, index_r = parse_index(index_path, region)
    index_f = merge_dicts(index_f, splice_f)
    index_r = merge_dicts(index_r, splice_r)
    index_f = sort_indeces(index_f)
    index_r = sort_indeces(index_r)

    # Get all translation blocks
    blocks_f, blocks_r = get_translation_blocks(index_f, index_r)

    # Get all putative internal and terminal positions
    exon_l_term_dict, exon_r_term_dict = get_putative_exons(blocks_f, blocks_r)

    # Check which terminal exons are in frame
    all_dict, final_dict = terminal_exons_in_phase(exon_l_term_dict, exon_r_term_dict)

    # compare to the reference
    l_exact = 0
    r_exact = 0
    b_exact = 0
    for i in introns:
        ref_exons = ref_dict[i]
        final_exons = final_dict[i]
        all_exons = all_dict[i]
        gene = ref_dict[i][0][0], ref_dict[i][0][1], ref_dict[i][1][2], ref_dict[i][0][3]
        if i[3] == '+':
            if length(ref_exons[0]) < length(final_exons[0]):
                l = 'S'
            elif length(ref_exons[0]) > length(final_exons[0]):
                l = 'L'
            else:
                l = 'E'
            if length(ref_exons[1]) < length(final_exons[1]):
                r = 'S'
            elif length(ref_exons[1]) > length(final_exons[1]):
                r = 'L'
            else:
                r = 'E'
            line = fr(gene), fr(ref_exons[0]), fr(final_exons[0]), l, fr(ref_exons[1]), fr(final_exons[1]), r, len(all_exons[0]), len(all_exons[1])
            if ref_exons[0] == final_exons[0]:
                l_exact += 1
            if ref_exons[1] == final_exons[1]:
                r_exact += 1
            if ref_exons[0] == final_exons[0] and ref_exons[1] == final_exons[1]:
                b_exact += 1
        else:
            if length(ref_exons[0]) < length(final_exons[0]):
                r = 'S'
            elif length(ref_exons[0]) > length(final_exons[0]):
                r = 'L'
            else:
                r = 'E'
            if length(ref_exons[1]) < length(final_exons[1]):
                l = 'S'
            elif length(ref_exons[1]) > length(final_exons[1]):
                l = 'L'
            else:
                l = 'E'
            line = fr(gene), fr(ref_exons[1]), fr(final_exons[1]), l, fr(ref_exons[0]), fr(final_exons[0]), r, len(all_exons[1]), len(all_exons[0])
            if ref_exons[0] == final_exons[0]:
                r_exact += 1
            if ref_exons[1] == final_exons[1]:
                l_exact += 1
            if ref_exons[0] == final_exons[0] and ref_exons[1] == final_exons[1]:
                b_exact += 1
        # gene, 5' exon (ref), 5' exon, 3' exon (ref), 3' exon, No.5' exons, No.3' exons
        print(*line, sep='\t')

        print_as_gff3(all_exons[0] + all_exons[1], 'test.gff3', kind='CDS')

    pprint('Total two-exon transcripts = {:,}'.format(len(introns)))
    pprint('5 end matches = {:,}'.format(l_exact))
    pprint('3 end matches = {:,}'.format(r_exact))
    pprint('Both ends matches = {:,}'.format(b_exact))
