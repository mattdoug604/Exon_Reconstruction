#!/home2/mattdoug/python3/bin/python3
# Last updated: 12/4/2018
# Author: Matt Douglas
#
# PURPOSE: Reconstruct protein-coding exons from a reference genome and RNA-seq
# data.
#
# UPDATES: In this version, exons are scored based on the score of the adjacent
# intron(s).
#
# MY NOTES:
# I:217,489..219,488 - example of resolved 2-exon gene
# X:2,412,849..2,413,848 - example of incorrect 2-exon gene
# III:3,562,617..3,582,616 - PacBio reads support alt intron w/ 2 support
# I:5,564,625..5,565,624 - missing internal exon due to phase incompatibility?
# X:13,397,323..13,397,775 - gene with 2 introns, but only 1 is part of coding exons
# II:1,889,845..1,890,844 - what's going on here?
# I:506,491..508,490 - False intron retention events?
# II:8,042,500..8,043,100 - intron retention in 5' terminal exon
# II:1,918,871-1,919,297 - 1bp "intron" introduces stop codon
# III:5,547,132-5,547,217 - Missed terminal exon (fixed)
# X:48,051-48,995 - Two genes, same strand, overlap at the ends - not in frame with each other so the ends get (falsely) discarded as UTRs
# II:6,490,579-6,491,167 - complex region w/ two "internal" 5' terminal exons
# I:348,487..353,486 - terminal exon does not have adjacent intron at the end
#
# TO DO:
# V:14895887..14895997 - intron retention in terminal exons, how to deal with?
#
# NOTS:
# Tried 'bedtools coverage' to get read depth for intron retention events, but
# it seems to be just as slow (if not slower)

from collections import defaultdict
from itertools import chain, product
import argparse, re, sys
import pysam
from main import parse_arguments
from parse_exon_index import parse_index
from introns_from_gff3 import get_introns
import intron_retention

#####################
# Class definitions #
#####################
class TranslationBlock(object):
    """A translation block is an object that represents the region from one stop
    codon to the next stop codon in a reading frame of the genome. Each
    TranslationBlock object tracks the start and end positions of the block, as
    well as all splice sites and putative start codons that lie within."""
    __slots__ = ['start', 'end', 'l_sites', 'r_sites', 's_sites']
    def __init__(self):
        self.start = -1
        self.end = -1
        self.l_sites = [] # 5' splice sites
        self.r_sites = [] # 3' splice sites
        self.s_sites = [] # start codons

    def num_sites(self):
        return sum([sum(i) for i in (self.l_sites, self.r_sites, self.s_sites)])


#####################
# Utility functions #
#####################
def pprint(*args, **kwargs):
    """Print runtime information at various 'levels' that can be turned on or
    off."""
    level = kwargs.pop('level', {'level':None})
    if level == 'debug':
        if DEBUG:
            print(*args, **kwargs)
    elif level == 'progress':
        if True not in (QUIET, NOPROG, DEBUG):
            print(*args, **kwargs)
    else:
        if not QUIET:
            print(*args, **kwargs)


def percent(val, total):
    if total == 0:
        return '0%'
    per = val * 100 / total
    per = '{0:.1f}%'.format(per)

    return per


def merge_dicts(x, y):
    """Merge to dictionaries (x, y) with the structure: key:[[], [], []] into a
    new dictionary (z) with the same structure."""
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


def frame_shift(adj_frame, intron, direction):
    """Calculate the frame a given exon uses based on the frame of the adjacent
    exon and the lenght of the intervening intron.
    """
    i_len = intron[2]-intron[1]+1
    mod = i_len % 3

    if direction in ('+<', '->'):
        new_frame = (adj_frame - mod) % 3
    elif direction in ('+>', '-<'):
        new_frame = (adj_frame + mod) % 3

    pprint('    Calculating frame shift:', level='debug')
    pprint('      Direction is {}'.format(direction), level='debug')
    pprint('      Adjacent exon is in frame {}'.format(adj_frame), level='debug')
    pprint('      Intron {} is {}bp long, modulo is {}'.format(intron, i_len, mod), level='debug')
    pprint('      New frame is {}'.format(new_frame), level='debug')

    return new_frame


def sort_indeces(index_dict):
    """Take a dictionary of tuples, remove duplicates and sort by value."""
    pprint('\rRemoving duplicates and sorting by position...', end='', level='debug')

    for chrom, frames in index_dict.items():
        for n, f in enumerate(frames):
            index_dict[chrom][n] = sorted(set(f))

    pprint('\rRemoving duplicates and sorting by position... Done!', level='debug')

    return index_dict


def sort_by_pos(features):
    """Sort tuples of exons by chromosome, then start position, then end
    position. NOTE: Wormbase uses roman numerals for chromosomes names.
    """
    numerals = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'X':10, 'MtDNA':11}
    try:
        return sorted(features, key=lambda x: (numerals[x[0]], int(x[1]), int(x[2])))
    except KeyError:
        return sorted(features, key=lambda x: (x[0], int(x[1]), int(x[2])))


def calc_score(exon):
    """Exon score is based on the scores of the adjacent introns. Sum the scores
    of all the adjacent introns upstream of the exon, then do the same for the
    adjacent introns downstream of the exon. The exon score is the minimum of
    these two values."""
    chrom, left, right, strand = exon
    left_introns = SITE[(chrom, left-1, strand)]
    right_introns = SITE[(chrom, right+1, strand)]

    if len(left_introns) * len(right_introns) > 0:
        l_sum = sum(INTRON[i] for i in left_introns)
        r_sum = sum(INTRON[i] for i in right_introns)
        return min(l_sum, r_sum)
    elif len(left_introns) > 0:
        return sum(INTRON[i] for i in left_introns)
    elif len(right_introns) > 0:
        return sum(INTRON[i] for i in right_introns)
    else:
        # return 0
        raise KeyError('No adjacent intron(s) found for exon: {}:{}-{}{}'.format(*exon))
        sys.exit(1)


def print_as_gff3(features, out_path, kind='CDS', mode='w'):
    """Take a list of features as a tuple (chromosome, start, end, strand),
    and print them as GFF3 formatted entries.
    """
    pprint('Printing results to "{}"'.format(out_path), level='debug')

    with open(out_path, mode) as f:
        if 'w' in mode: # print the header only if we're opening a new file
            print('##gff-version 3', file=f)
        for n, feat in enumerate(features):
            chrom, left, right, strand = feat
            if kind == 'CDS':
                score = calc_score(feat)
            elif kind == 'intron':
                score = INTRON[feat]
            else:
                score = '.'
            phase = '.'
            frame = ','.join(map(str, FRAME[feat]))
            notes = [';frame='+frame if frame != '' else ''][0]
            line = chrom, '.', kind, left, right, score, strand, phase, 'ID={}{}{}'.format(kind, n+1, notes)
            print(*line, sep='\t', file=f)


##################
# Main functions #
##################
def get_translation_blocks(index_f, index_r):
    """Return two dictionaries of TranslationBlock objects. A 'translation block'
    is defined as the region from the first base of the first codon after a stop
    codon to the last base of the next stop codon in the same reading frame.

    INPUT: Two dictionaries corresponding to the indexed positions of all stop
    codons, putative start codons, and splice sites on the forward and reverse
    strands of the genome, respectively. Keys in each dict are the chromosomes.
    The values are 3 lists corresponding to each reading frame:
    index_f = { chromsome: [ [frame0], [frame1], [frame2] ] }
    index_r = { chromsome: [ [frame0], [frame1], [frame2] ] }

    RETURN: Two dictionaries (blocks_f, blocks_r), structured the same as the
    input where each of the three lists (for each key) is a list of
    TranslationBlocks objects in that reading frame.
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
                    if len(block.l_sites) > 0: pprint("    5' splicesites:", ' '.join(map(str, block.l_sites)), level='debug')
                    if len(block.r_sites) > 0: pprint("    3' splicesites:", ' '.join(map(str, block.r_sites)), level='debug')
                    if len(block.s_sites) > 0: pprint("    start sites:", ' '.join(map(str, block.s_sites)), level='debug')
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
                    if len(block.l_sites) > 0: pprint("    5' splicesites:", ' '.join(map(str, block.l_sites)), level='debug')
                    if len(block.r_sites) > 0: pprint("    3' splicesites:", ' '.join(map(str, block.r_sites)), level='debug')
                    if len(block.s_sites) > 0: pprint("    start sites:", ' '.join(map(str, block.s_sites)), level='debug')
                    # reset
                    num_blks += 1
                    last = pos
                    block = TranslationBlock()

    pprint('\rIdentifying all translation blocks... Done!    ')
    pprint('  {:,} translation blocks'.format(num_blks))

    # DEBUG: output all the translation blocks
    temp_f = merge_dicts(blocks_f, empty_blocks_f)
    temp_r = merge_dicts(blocks_r, empty_blocks_r)
    out_path = '{}.translation_blocks.gff3'.format(PREFIX)
    #print_as_gff3([], out_path, mode='w') # create a new file
    for d, strand, suf in ((temp_f, '+', 'plus'), (temp_r, '-', 'minus')):
        temp = [ [], [], [] ]
        for chrom, frames in d.items():
            for frame, blocks in enumerate(frames):
                temp[frame] += [(chrom, i.start, i.end, strand) for i in blocks]
        for frame, blocks in enumerate(temp):
            blocks = sort_by_pos(blocks)
            #print_as_gff3(blocks, out_path, kind=suf+str(frame), mode='a')

    return blocks_f, blocks_r


def get_putative_exons(blocks_f, blocks_r):
    """Search each translation block and return the genomic coordinates of all
    possible internal and terminal exons, which reading frame(s) each use, and
    which intron(s) the exons are adjacent to.

    INPUT: Two dictionaries of TranslationBlock objects on forward and reverse
    strands of the genome, respectively. Keys for each dict are the chromosomes.
    Values are three lists corresponding to each reading frame:
    blocks_f = { chromsome: [ [frame0], [frame1], [frame2] ] }
    blocks_r = { chromsome: [ [frame0], [frame1], [frame2] ] }

    RETURN:
    1) a set of "internal exons" as tuples in the format: (chromosome, start,
       end, strand).
    2) two dictionaries of terminal exons (one for 5' terminal exons, one for 3'
       terminal) where the key is the adjacent intron(s) and the value is a set
       of exon positions as tuples (chromosome, start, end, strand).
    """
    exon_internal = set()
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
                l_sites = sorted(block.l_sites)
                r_sites = sorted(block.r_sites)
                s_sites = sorted(block.s_sites)
                pprint('  Checking block: {}:{:,}-{:,} (+ strand)'.format(chrom, block.start, block.end), level='debug')
                pprint('  l_sites:', l_sites, level='debug')
                pprint('  r_sites:', r_sites, level='debug')
                pprint('  s_sites:', s_sites, level='debug')
                # join splice site to splice site
                for l in l_sites:
                    x = [r for r in r_sites if r > l]
                    for i, j in product([l], x):
                        exon = chrom, i, j-1, '+'
                        exon_internal.add(exon)
                        FRAME[exon].add(frame)
                        pprint('    internal exon found at {}'.format(exon), level='debug')
                # join start codons to 5' splice sites
                for s in s_sites:
                    x = [r for r in r_sites if r > s]
                    for i, j in product([s], x):
                        exon = chrom, i, j-1, '+'
                        r_adj = chrom, j-1, '+'
                        exon_l_term_dict[r_adj].add(exon)
                        FRAME[exon].add(frame)
                        five_term_c += 1
                        pprint("    5' terminal exon found at {}".format(exon), level='debug')
                # join 3' splice sites to stop codons
                for l in [s for s in l_sites if s < block.end]:
                    exon = chrom, l, block.end, '+'
                    l_adj = chrom, l, '+'
                    exon_r_term_dict[l_adj].add(exon)
                    FRAME[exon].add(frame)
                    three_term_c += 1
                    pprint("    3' terminal exon found at {}".format(exon), level='debug')
    ############################################################################
    # scan the - strand of the genome
    for chrom, frames in blocks_r.items():
        for frame, i in enumerate(frames):
            pprint('\rIdentifying putative exon positions...', percent(count, total), end='', level='progress')
            pprint('\nScanning - strand, reading frame {}'.format(frame), level='debug')
            count += 1
            for block in i:
                l_sites = sorted(block.l_sites)
                r_sites = sorted(block.r_sites)
                s_sites = sorted(block.s_sites)
                pprint('  Checking block: {}:{:,}-{:,} (- strand)'.format(chrom, block.start, block.end), level='debug')
                pprint('  l_sites:', l_sites, level='debug')
                pprint('  r_sites:', r_sites, level='debug')
                pprint('  s_sites:', s_sites, level='debug')
                # join splice site to splice site
                for l in l_sites:
                    x = [r for r in r_sites if r > l]
                    for i, j in product([l], x):
                        exon = chrom, i+1, j, '-'
                        exon_internal.add(exon)
                        FRAME[exon].add(frame)
                        pprint('    internal exon found at {}'.format(exon), level='debug')
                # join stop codon to 3' splice sites
                for r in [s for s in r_sites if s > block.start]:
                    exon = chrom, block.start, r, '-'
                    r_adj = chrom, r, '-'
                    exon_l_term_dict[r_adj].add(exon)
                    FRAME[exon].add(frame)
                    three_term_c += 1
                    pprint("    5' terminal exon found at {}".format(exon), level='debug')
                # join 5' splice sites to start codons
                for l in l_sites:
                    x = [s for s in s_sites if s > l]
                    for i, j in product([l], x):
                        exon = chrom, i+1, j, '-'
                        l_adj = chrom, i+1, '-'
                        exon_r_term_dict[l_adj].add(exon)
                        FRAME[exon].add(frame)
                        five_term_c += 1
                        pprint("    3' terminal exon found at {}".format(exon), level='debug')

    pprint('\rIdentifying putative exon positions... Done!    ')
    pprint("  {:,} putative internal exons".format(len(exon_internal)))
    pprint("  {:,} putative 5' terminal exons".format(five_term_c))
    pprint("  {:,} putative 3' terminal exons".format(three_term_c))

    exon_internal = list(exon_internal)

    return exon_internal, exon_l_term_dict, exon_r_term_dict


def check_adjacency(exon_list):
    """Given a list of exons and a dictionary of splice sites and their
    corresponding introns, return a dictionary where each key is an intron, and
    the value is two lists of exons that are adjacent to the intron - the first
    list is exons at the 5' end, the second list is exons at the 3' end."""
    global ADJ

    pprint('Building adjacency dict...', level='debug')

    for exon in exon_list:
        chrom, left, right, strand = exon
        l_adj = chrom, left - 1, strand
        r_adj = chrom, right + 1, strand
        for intron in SITE[l_adj]:
            ADJ[intron][1].add(exon)
            pprint('  exon is right of intron {}'.format(intron), level='debug')
        for intron in SITE[r_adj]:
            ADJ[intron][0].add(exon)
            pprint('  exon is left of intron {}'.format(intron), level='debug')

    pprint('Building adjacency dict... Done!', level='debug')


def resolve_internal_frames(exon_list, max_iter=20):
    """Sometimes a putative internal exon may appear in more than one reading
    frame. Here we attempt to resolve any abiguous frames by iteratively
    checking the frame of adjacent exons, and using those to calulate what the
    frame should be."""
    unresolved = []
    count = 0
    total = len(exon_list)
    total_change = 0

    pprint('', level='debug')
    pprint('Identifying the correct phase for internal exons', end='')

    for n in range(max_iter):
        pprint('\nStarting iteration {}:'.format(n+1), level='debug')
        new_FRAME = defaultdict(set)
        num_change = 0
        for exon in exon_list:
            if len(FRAME[exon]) > 1:
                pprint('Iteration {}) Block {} has more than one frame ({}), skipping'.format(n+1, exon, FRAME[exon]), level='debug')
                continue
            chrom, left, right, strand = exon
            (frame,) = FRAME[exon]
            left_introns = SITE[(chrom, left-1, strand)]
            right_introns = SITE[(chrom, right+1, strand)]
            pprint('Iteration {}) Using exon {} (frame {}) to resolve adjacent reading frames'.format(n+1, exon, frame), level='debug')
            ####################################################################
            for intron in left_introns:
                pprint('  Looking for exons adjacent to upstream intron {}...'.format(intron), level='debug')
                for adj_exon in ADJ[intron][0]:
                    adj_frames = FRAME[adj_exon]
                    if len(adj_frames) != 1:
                        pprint('  Determining frame for exon {}, frames: {}'.format(adj_exon, adj_frames), level='debug')
                        new_frame = frame_shift(frame, intron, strand+'<')
                        new_FRAME[adj_exon].add(new_frame)
            ####################################################################
            for intron in right_introns:
                pprint('  Looking for exons adjacent to downstream intron {}...'.format(intron), level='debug')
                for adj_exon in ADJ[intron][1]:
                    adj_frames = FRAME[adj_exon]
                    if len(adj_frames) != 1:
                        pprint('  Determining frame for exon {}, frames {}:'.format(adj_exon, adj_frames), level='debug')
                        new_frame = frame_shift(frame, intron, strand+'>')
                        new_FRAME[adj_exon].add(new_frame)

        # update reading frames
        pprint('\nSummary of frame changes:', level='debug')
        for exon, new_frames in new_FRAME.items():
            old_frames = FRAME[exon]
            if old_frames != new_frames:
                pprint('  Block {} frames {} => {}'.format(exon, old_frames, new_frames), level='debug')
                FRAME[exon] = new_frames
                num_change += 1
        total_change += num_change
        if num_change < 1:
            pprint('  None. Stopping at iteration {}'.format(n+1), level='debug')
            break

    # check for exons that still have more than one valid frame
    pprint('\nExons with unknown phase:', level='debug')
    for exon, frames in FRAME.items():
        if len(frames) > 1:
            unresolved.append(exon)
            pprint('  {}, frames: {}'.format(exon, frames), level='debug')

    pprint('\rIdentifying the correct phase for internal exons... Done!')
    pprint('  Resolved {:,} case(s) where an internal exon had >1 possible phase'
           .format(total_change))
    pprint('  The phase of {:,} ({}) internal exons could not be determined'
           .format(len(unresolved), percent(len(unresolved), len(exon_list))))

    # DEBUG: Output internal exons with ambiguous frame
    #print_as_gff3(unresolved, PREFIX+'.internal.multi_frame.gff3')

    return unresolved


def terminal_exons_in_phase(exon_internal, exon_l_term_dict, exon_r_term_dict, max_iter=3):
    """Putative terminal exons can be identified at many loci. Only report ones
    that A) do not overlap internal exons and B) are in-frame with an adjacent
    exon.

    If exons could not be indentified on one (or both) sides of a given intron,
    that intron is ignored, we remove any internal exons that were previously
    indentifed on the one side and instead look for terminal exons at those
    positions."""
    exon_five_term = set()
    exon_three_term = set()
    intron_discard = set()
    exon_discard = set()
    global ADJ
    global SITE
    # for keeping track of progress
    count = 0
    total = len(ADJ)
    # quick function to calculate the length of an exon
    length = lambda x: x[2] - x[1] + 1
    # for DEBUG-ing
    single_intron = []

    pprint('', level='debug')
    pprint('Checking which terminal exons are in phase', end='')

    for n in range(max_iter):
        pprint('\nStarting iteration {}:'.format(n+1), level='debug')
        new_adj = {}
        for intron, adj in ADJ.items():
            pprint('\rChecking which terminal exons are in phase...',
                   percent(count/max_iter, total), end='', level='progress')
            chrom, l, r, strand = intron
            l_adj, r_adj = adj
            left_exons = exon_l_term_dict[(chrom, l-1, strand)]
            right_exons = exon_r_term_dict[(chrom, r+1, strand)]
            l_compare = {}
            r_compare = {}
            expected = set()
            results = []
            match = None
            count += 1
            ####################################################################
            if len(l_adj) < 1 and len(r_adj) < 1:
                pprint('Looking on both sides of intron at {}'.format(intron), level='debug')
                new_adj[intron] = [ set(), set() ]
                # which exons are in which frames, on both sides
                for exon in left_exons:
                    for f in FRAME[exon]:
                        if f not in l_compare:
                            l_compare[f] = set()
                        l_compare[f].add(exon)
                for exon in right_exons:
                    for f in FRAME[exon]:
                        if f not in r_compare:
                            r_compare[f] = set()
                        r_compare[f].add(exon)
                pprint('  Left side:', level='debug')
                for frame, exons in l_compare.items():
                    pprint('    Frame {}) {}'.format(frame, ' '.join([str(i) for i in exons])), level='debug')
                pprint('  Right side:', level='debug')
                for frame, exons in r_compare.items():
                    pprint('    Frame {}) {}'.format(frame, ' '.join([str(i) for i in exons])), level='debug')
                # check which exons are in frame with each other
                # check if they are a valid ORF (i.e. modulo 3 == 0)
                for f, l_set in l_compare.items():
                    new_f = frame_shift(f, intron, strand+'>')
                    if new_f in r_compare:
                        r_set = r_compare[new_f]
                    else:
                        continue
                    for l_exon in l_set:
                        for r_exon in r_set:
                            if ( length(l_exon) + length(r_exon)) % 3 == 0:
                                results.append((l_exon, r_exon))
                # if multiple are valid, select the longest ORF
                if len(results) > 0:
                    final = max(results, key=lambda x: (x[0][2]-x[0][1] + x[1][2]-x[1][1]))
                    new_adj[intron][0].add(final[0])
                    new_adj[intron][1].add(final[1])
                    if strand == '+':
                        exon_five_term.add(final[0])
                        exon_three_term.add(final[1])
                    elif strand == '-':
                        exon_five_term.add(final[1])
                        exon_three_term.add(final[0])
                    pprint('  Longest ORF is:', level='debug')
                    pprint('    {} is left of intron {}'.format(final[0], intron), level='debug')
                    pprint('    {} is right of intron {}'.format(final[1], intron), level='debug')
                    # DEBUG: print out any single-intron genes
                    single_intron.append(final[0])
                    single_intron.append(final[1])
                else:
                    pprint('  No valid ORF found! Delete intron', intron, level='debug')
                    del new_adj[intron]
                    intron_discard.add(intron)
            ####################################################################
            elif len(l_adj) < 1:
                pprint('Looking upstream of intron at {}'.format(intron), level='debug')
                new_adj[intron] = [ set(), r_adj ]
                # determine what the expected frame(s) should be based on
                # adjacent internal exons
                for exon in r_adj:
                    for f in FRAME[exon]:
                        new_f = frame_shift(f, intron, strand+'<')
                        pprint('  Downstream exon {}, frame {}'.format(exon, f), level='debug')
                        expected.add(new_f)
                pprint('  Expect adjacent frames to be one of: {}'.format(expected), level='debug')
                # only keep the longest exon(s) that match the expected frame(s)
                for exon in left_exons:
                    for f in FRAME[exon]:
                        if f not in l_compare:
                            l_compare[f] = set()
                        l_compare[f].add(exon)
                for f in expected:
                    if f in l_compare:
                        match = max(l_compare[f], key=length)
                        new_adj[intron][0].add(match)
                        if strand == '+':
                            exon_five_term.add(match)
                        elif strand == '-':
                            exon_three_term.add(match)
                        pprint('  Blocks {} with frame {} match'
                               .format(', '.join([str(i) for i in l_compare[f]]), f), level='debug')
                        pprint('  Longest is {}'.format(match), level='debug')
                if match is None:
                    pprint('  No valid upstream exon found! Delete intron', intron, level='debug')
                    del new_adj[intron]
                    intron_discard.add(intron)
            ####################################################################
            elif len(r_adj) < 1:
                pprint('Looking downstream of intron at {}'.format(intron), level='debug')
                new_adj[intron] = [ l_adj, set() ]
                # determine what the expected frame(s) should be based on
                # adjacent internal exons
                for exon in l_adj:
                    for f in FRAME[exon]:
                        new_f = frame_shift(f, intron, strand+'>')
                        pprint('  Upstream exon {}, frame {}'.format(exon, f), level='debug')
                        expected.add(new_f)
                pprint('  Expect adjacent frames to be one of: {}'.format(expected), level='debug')
                # only keep the longest exon(s) that match the expected frame(s)
                for exon in right_exons:
                    for f in FRAME[exon]:
                        if f not in r_compare:
                            r_compare[f] = set()
                        r_compare[f].add(exon)
                for f in expected:
                    if f in r_compare:
                        match = max(r_compare[f], key=length)
                        new_adj[intron][1].add(match)
                        if strand == '+':
                            exon_three_term.add(match)
                        elif strand == '-':
                            exon_five_term.add(match)
                        pprint('  Blocks {} with frame {} match'
                               .format(', '.join([str(i) for i in r_compare[f]]), f), level='debug')
                        pprint('  Longest is {}'.format(match), level='debug')
                if match is None:
                    pprint('  No valid downstream exon found! Delete intron', intron, level='debug')
                    del new_adj[intron]
                    intron_discard.add(intron)
            ####################################################################
            else:
                pprint('Intron at {} is already flanked by exons. Continue.'.format(intron), level='debug')
                new_adj[intron] = [ l_adj, r_adj ]

        if new_adj == ADJ:
            pprint('No change. Stopping iterations.', level='debug')
            break

        pprint('Updating adjacency dictionary...', level='debug')
        ADJ = new_adj
        SITE = defaultdict(list)
        for intron in new_adj:
            chrom, left, right, strand = intron
            SITE[(chrom, left, strand)].append(intron)
            SITE[(chrom, right, strand)].append(intron)

        # remove exons from the ADJ dictionary
        pprint('Removing internal exons without adjacent introns...', level='debug')
        temp = set()
        for exon in exon_internal:
            chrom, left, right, strand = exon
            left_introns = SITE[(chrom, left-1, strand)]
            right_introns = SITE[(chrom, right+1, strand)]
            if len(left_introns) * len(right_introns) > 0:
                pprint('  Keeping exon', exon, level='debug')
                temp.add(exon)
            else:
                pprint('  Discarding exon', exon, level='debug')
                exon_discard.add(exon)
                # remove the exon from the adj dict
                for intron in left_introns + right_introns:
                    if intron in ADJ:
                        ADJ[intron][0].discard(exon)
                        ADJ[intron][1].discard(exon)
        exon_internal = temp

        pprint("Removing 5' terminal exons exons without adjacent introns...", level='debug')
        temp = set()
        for exon in exon_five_term:
            chrom, left, right, strand = exon
            if strand == '+':
                right_introns = SITE[(chrom, right+1, strand)]
                if len(right_introns) > 0:
                    pprint('  Keeping exon', exon, level='debug')
                    temp.add(exon)
                else:
                    pprint('  Discarding exon', exon, level='debug')
                    exon_discard.add(exon)
                    # remove the exon from the adj dict
                    for intron in right_introns:
                        if intron in ADJ:
                            ADJ[intron][0].discard(exon)
            elif strand == '-':
                left_introns = SITE[(chrom, left-1, strand)]
                if len(left_introns) > 0:
                    pprint('  Keeping exon', exon, level='debug')
                    temp.add(exon)
                else:
                    pprint('  Discarding exon', exon, level='debug')
                    exon_discard.add(exon)
                    # remove the exon from the adj dict
                    for intron in left_introns:
                        if intron in ADJ:
                            ADJ[intron][1].discard(exon)
        exon_five_term = temp

        pprint("Removing 3' terminal exons exons without adjacent introns...", level='debug')
        temp = set()
        for exon in exon_three_term:
            chrom, left, right, strand = exon
            if strand == '+':
                left_introns = SITE[(chrom, left-1, strand)]
                if len(left_introns) > 0:
                    pprint('  Keeping exon', exon, level='debug')
                    temp.add(exon)
                else:
                    pprint('  Discarding exon', exon, level='debug')
                    exon_discard.add(exon)
                    # remove the exon from the adj dict
                    for intron in left_introns:
                        if intron in ADJ:
                            ADJ[intron][1].discard(exon)
            elif strand == '-':
                right_introns = SITE[(chrom, right+1, strand)]
                if len(right_introns) > 0:
                    pprint('  Keeping exon', exon, level='debug')
                    temp.add(exon)
                else:
                    pprint('  Discarding exon', exon, level='debug')
                    exon_discard.add(exon)
                    # remove the exon from the adj dict
                    for intron in right_introns:
                        if intron in ADJ:
                            ADJ[intron][0].discard(exon)
        exon_three_term = temp

    pprint('\rChecking which terminal exons are in phase... Done!  ')
    pprint("  {:,} 5' terminal exons".format(len(exon_five_term)))
    pprint("  {:,} 3' terminal exons".format(len(exon_three_term)))

    # print_as_gff3(single_intron, PREFIX+'.double.gff3') # DEBUG
    # print_as_gff3(intron_discard, PREFIX+'.non_coding_introns.gff3', kind='intron') # DEBUG
    # print_as_gff3(exon_discard, PREFIX+'.terminal_removed.gff3') # DEBUG

    exon_internal = list(exon_internal)
    exon_five_term = list(exon_five_term)
    exon_three_term = list(exon_three_term)

    return exon_internal, exon_five_term, exon_three_term, intron_discard


def define_gene_ends(features):
    """Return gene boundaries (i.e. contigous regions of internal exons,
    terminal exons, and introns.
    """
    flatten = chain.from_iterable
    temp1 = defaultdict(list)
    temp2 = defaultdict(list)
    ends_f = {}
    ends_r = {}
    genes = []

    pprint('', level='debug')
    pprint('Defining putative gene boundaries', end='')

    # merge sets of introns and internal exons
    for i in features:
        chrom, l, r, strand = i
        temp1[(chrom, strand)].append((l, r))

    # find contiguous regions of intron <- internal exon -> intron on each
    # chromosome on each strand
    for region, ranges in temp1.items():
        ranges = sorted(flatten(((l, 1), (r, -1)) for l, r in ranges))
        c, x = 0, 0
        for value, label in ranges:
            if c == 0:
                x = value
            c += label
            if c == 0:
                temp2[region].append((x, value))

    # treat the ends of genes as "pseudo-splice sites" that we'll use to find
    # single-exon genes
    for region, ranges in temp2.items():
        chrom, strand = region
        for i in ranges:
            left, right = i
            d = [ends_f if strand == '+' else ends_r][0]
            genes.append((chrom, left, right, strand))
            for n in range(3):
                try:
                    d[chrom][n].append((left, 2))
                    d[chrom][n].append((right, 3))
                except KeyError:
                    d[chrom] = [ [], [], [] ]
                    d[chrom][n].append((left, 2))
                    d[chrom][n].append((right, 3))

    pprint('\rDefining putative gene boundaries... Done!')
    pprint('  {:,} multi-exon genes found'.format( sum([len(i) for i in temp2.values()])))

    for d, s in ((ends_f, '+'), (ends_r, '-')):
        for chrom, frames in d.items():
            for pos, kind in frames[0]:
                if kind == 2:
                    pprint('left side of gene at: {}:{:,} ({} strand)'.format(chrom, pos, s), level='debug')
                elif kind == 3:
                    pprint('right side of gene at: {}:{:,} ({} strand)'.format(chrom, pos, s), level='debug')

    #print_as_gff3(genes, PREFIX+'.genes.gff3', kind='gene')

    return ends_f, ends_r, genes


def get_single_exons(index_f, index_r):
    """Single exons are defined as regions from putative start codon to stop
    codon within the same translation block."""
    exon_single = set()
    count = 0
    total = len(index_f)*6 #calculate the number of positions to check

    pprint('', level='debug')
    pprint('Checking intergenic regions for single exons', end='')

    for chrom, frames in index_f.items():
        for frame, i in enumerate(frames):
            pprint('\rChecking intergenic regions for single exons...', percent(count, total), end='', level='progress')
            pprint('\nScanning + strand in the > direction, reading frame {}'.format(frame), level='debug')
            count += 1
            track_list = []
            intergenic = True
            for pos, kind in i:
                # When a start codon is reached...
                if kind == 0:
                    if intergenic is True:
                        track_list.append((pos, 0))
                # When a stop codon is reached...
                elif kind == 3:
                    for t, x in track_list:
                        # if we're tracking from a start codon...
                        if x == 0:
                            exon = chrom, t, pos-1, '+'
                            exon_single.add(exon)
                            FRAME[exon].add(frame)
                            pprint('  single-exon exon found at {}'.format(exon), level='debug')
                    track_list = [] # reset
                # When a left splice site is reached...
                elif kind == 1:
                    track_list = []
                    intergenic = False
                # When a right splice site is reached...
                elif kind == 2:
                    track_list = []
                    intergenic = True
    ############################################################################
    for chrom, frames in index_r.items():
        for frame, i in enumerate(frames):
            pprint('\rChecking intergenic regions for single exons...', percent(count, total), end='', level='progress')
            pprint('\nScanning - strand in the > direction, reading frame {}'.format(frame), level='debug')
            count += 1
            track_list = []
            intergenic = True
            for pos, kind in i:
                # When a start codon is reached...
                if kind == 0:
                    for t, x in track_list:
                        # if we're tracking from a stop codon...
                        if x == 0:
                            exon = chrom, t, pos, '-'
                            exon_single.add(exon)
                            FRAME[exon].add(frame)
                            pprint('  single-exon exon found at {}'.format(exon), level='debug')
                # When a stop codon is reached...
                elif kind == 3:
                    if intergenic:
                        track_list = [ (pos, 0) ]
                # When a left splice site is reached...
                elif kind == 1:
                    track_list = []
                    intergenic = False
                # When a right splice site is reached...
                elif kind == 2:
                    track_list = []
                    intergenic = True

    pprint('\rChecking intergenic regions for single exons... Done!    ')
    pprint('  {:,} single exons'.format(len(exon_single)))

    exon_single = list(exon_single)

    return exon_single


############
# Run loop #
############
if __name__ == '__main__':
    args = parse_arguments()

    QUIET = args.quiet
    NOPROG = args.noprog
    DEBUG = args.debug
    PREFIX = args.prefix
    FRAME = defaultdict(set)
    NOTE = defaultdict(set)

    # Parse the index files (positions of start and stop codons) #
    ##############################################################
    index_f, index_r = parse_index(args)

    # Get the positions of all the introns #
    ########################################
    INTRON, splice_f, splice_r, SITE = get_introns(args)
    ADJ = {i:[set(), set()] for i in INTRON}

    # Get all translation blocks #
    ##############################
    index_f = merge_dicts(index_f, splice_f)
    index_r = merge_dicts(index_r, splice_r)
    index_f = sort_indeces(index_f)
    index_r = sort_indeces(index_r)
    blocks_f, blocks_r = get_translation_blocks(index_f, index_r)

    # Get all putative internal and terminal positions #
    ####################################################
    exon_internal, exon_l_term_dict, exon_r_term_dict = get_putative_exons(blocks_f, blocks_r)
    exon_internal, _ = intron_retention.filter(args, INTRON, exon_internal)
    check_adjacency(exon_internal)

    # Try and resolve ambiguous phase for internal exons #
    ######################################################
    num_unresolved = len(resolve_internal_frames(exon_internal))

    # Check which terminal exons are in frame #
    ###########################################
    exon_internal, exon_five_term, exon_three_term, intron_discard = terminal_exons_in_phase(exon_internal, exon_l_term_dict, exon_r_term_dict)
    for i in intron_discard:
        del INTRON[i]
    # do one more round to use the terminal exons to determine the phase of internal exons
    _ = resolve_internal_frames(exon_internal + exon_five_term + exon_three_term)

    # Define gene possible stuctures #
    #################################
    ends_f, ends_r, genes = define_gene_ends(list(INTRON) + exon_internal + exon_five_term + exon_three_term)

    # Get all single-exons #
    ########################
    # index_f = merge_dicts(index_f, ends_f)
    # index_r = merge_dicts(index_r, ends_r)
    # index_f = sort_indeces(index_f)
    # index_r = sort_indeces(index_r)
    # exon_single = get_single_exons(index_f, index_r)

    # Output the results #
    ######################
    pprint('Outputting results...', end='')
    exon_internal = sort_by_pos(exon_internal)
    exon_five_term = sort_by_pos(exon_five_term)
    exon_three_term = sort_by_pos(exon_three_term)
    exon_final = sort_by_pos(exon_internal + exon_five_term + exon_three_term)
    # exon_single = sort_by_pos(exon_single)
    print_as_gff3(exon_internal, PREFIX+'.internal.gff3')
    print_as_gff3(exon_five_term, PREFIX+'.five_term.gff3')
    print_as_gff3(exon_three_term, PREFIX+'.three_term.gff3')
    # print_as_gff3(exon_single, PREFIX+'.single.gff3')
    print_as_gff3(exon_final, PREFIX+'.final.gff3')
    pprint('\rOutputting results... Done!')
    pprint('  {:,} multi-exon genes, {:,} exons found'.format(len(genes), len(exon_final)))
    pprint('    {:,} internal exons'.format(len(exon_internal)))
    pprint("    {:,} 5' terminal exons".format(len(exon_five_term)))
    pprint("    {:,} 3' terminal exons".format(len(exon_three_term)))
    # pprint('  {:,} single exon genes'.format(len(exon_single)))

    # Finish up #
    #############
    pprint('Finsihed!')
