#!/home2/mattdoug/python3/bin/python3
# Last updated: 8/3/2018
# Author: Matt Douglas

# PURPOSE: Reconstruct protein-coding exons from a reference genome and RNA-seq
# data.

# UPDATES TO THIS VERSION:
# Fixed how the script deals with introns that have no upstream and/or
# downstream exons.

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

# ALSO:
# Tried Bedtools genomecov to get read depth for intron retention events, but
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


def frame_shift(adj_frame, intron, direction):
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
    pprint('      Intron {} is {}bp long, modulo is {}'.format(intron, intron_len, modulo), level='debug')
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
    position. NOTE: Wormbase uses roman numerals for chromosomes names.
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

    # DEBUG: output all the translation blocks
    temp_f = merge_dicts(blocks_f, empty_blocks_f)
    temp_r = merge_dicts(blocks_r, empty_blocks_r)
    for d, strand, suf in ((temp_f, '+', 'plus'), (temp_r, '-', 'minus')):
        temp = [ [], [], [] ]
        for chrom, frames in d.items():
            for frame, blocks in enumerate(frames):
                temp[frame] += [(chrom, i.start, i.end, strand) for i in blocks]
        for frame, blocks in enumerate(temp):
            blocks = sort_by_pos(blocks)
            out_path = '{}.t_blocks.{}{}.gff3'.format(PREFIX, suf, frame)
            print_as_gff3(blocks, out_path, kind='block')

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
    pprint('  Resolved {:,} case(s) where an internal exons had >1 possible phase'
           .format(total_change))
    pprint('  The phase of {:,} ({}) internal exons could not be determined'
           .format(len(unresolved), percent(len(unresolved), len(exon_list))))

    # DEBUG: Output internal exons with ambiguous frame
    print_as_gff3(unresolved, PREFIX + '.internal.multi_frame.gff3')

    return unresolved


def terminal_exons_in_phase(exon_internal, exon_l_term_dict, exon_r_term_dict, max_iter=3):
    # quick function to calculate the length of an exon
    length = lambda x: x[2] - x[1] + 1

    global ADJ
    global SITE
    exon_five_term = set()
    exon_three_term = set()
    five_term_c = 0
    three_term_c = 0
    intron_discard = set()
    exon_discard = set()
    # for keeping track of progress
    count = 0
    total = len(ADJ)

    single_intron = [] # DEBUG

    pprint('', level='debug')
    pprint('Checking which terminal exons are in phase', end='')

    for n in range(max_iter):
        pprint('\nStarting iteration {}:'.format(n+1), level='debug')
        new_ADJ = {}
        for intron, adj in ADJ.items():
            pprint('\rChecking which terminal exons are in phase...',
                   percent(count, total), end='', level='progress')
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
            ########################################################################
            if len(l_adj) < 1 and len(r_adj) < 1:
                pprint('Looking on both sides of intron at {}'.format(intron), level='debug')
                new_ADJ[intron] = [ set(), set() ]
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
                    new_ADJ[intron][0].add(final[0])
                    new_ADJ[intron][1].add(final[1])
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
                    pprint('  No valid ORF found!', level='debug')
                    del new_ADJ[intron]
                    intron_discard.add(intron)
            ########################################################################
            elif len(l_adj) < 1:
                pprint('Looking upstream of intron at {}'.format(intron), level='debug')
                new_ADJ[intron] = [ set(), r_adj ]
                # determine what the expected frame(s) should be based on adjacent
                # internal exons
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
                        new_ADJ[intron][0].add(match)
                        if strand == '+':
                            exon_five_term.add(match)
                        elif strand == '-':
                            exon_three_term.add(match)
                        pprint('  Blocks {} with frame {} match'
                               .format(', '.join([str(i) for i in l_compare[f]]), f), level='debug')
                        pprint('  Longest is {}'.format(match), level='debug')
                if match is None:
                    pprint('  No valid upstream exon found!', level='debug')
                    del new_ADJ[intron]
                    intron_discard.add(intron)
            ########################################################################
            elif len(r_adj) < 1:
                pprint('Looking downstream of intron at {}'.format(intron), level='debug')
                new_ADJ[intron] = [ l_adj, set() ]
                # determine what the expected frame(s) should be based on adjacent
                # internal exons
                for exon in l_adj:
                    for f in FRAME[exon]:
                        new_f = frame_shift(f, intron, strand+'>')
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
                        new_ADJ[intron][1].add(match)
                        if strand == '+':
                            exon_three_term.add(match)
                        elif strand == '-':
                            exon_five_term.add(match)
                        pprint('  Blocks {} with frame {} match'
                               .format(', '.join([str(i) for i in r_compare[f]]), f), level='debug')
                        pprint('  Longest is {}'.format(match), level='debug')
                if match is None:
                    pprint('  No valid downstream exon found!', level='debug')
                    del new_ADJ[intron]
                    intron_discard.add(intron)
            ########################################################################
            else:
                pprint('Intron at {} is already flanked by exons. Continue.'.format(intron), level='debug')
                new_ADJ[intron] = [ l_adj, r_adj ]

        if new_ADJ == ADJ:
            pprint('No change. Stopping iterations.', level='debug')
            break

        pprint('Updating adjacency dictionary...', level='debug')
        ADJ = new_ADJ
        SITE = defaultdict(list)
        for intron in new_ADJ:
            chrom, left, right, strand = intron
            SITE[(chrom, left, strand)].append(intron)
            SITE[(chrom, right, strand)].append(intron)

        pprint('Removing internal exons without adjacent introns...', level='debug')
        temp = []
        for exon in exon_internal:
            chrom, left, right, strand = exon
            left_introns = SITE[(chrom, left-1, strand)]
            right_introns = SITE[(chrom, right+1, strand)]
            if len(left_introns) * len(right_introns) > 0:
                pprint('  Keeping exon', exon, level='debug')
                temp.append(exon)
            else:
                pprint('  Discarding exon', exon, level='debug')
                exon_discard.add(exon)
                # remove the exon from the adj dict
                for intron in left_introns + right_introns:
                    if intron in ADJ:
                        ADJ[intron][0].discard(exon)
                        ADJ[intron][1].discard(exon)
        exon_internal = temp

    pprint('\rChecking which terminal exons are in phase... Done!  ')
    pprint("  {:,} 5' terminal exons".format(len(exon_five_term)))
    pprint("  {:,} 3' terminal exons".format(len(exon_three_term)))

    print_as_gff3(single_intron, PREFIX + '.double.gff3') # DEBUG
    print_as_gff3(intron_discard, PREFIX + '.non_coding_introns.gff3') # DEBUG
    print_as_gff3(exon_discard, PREFIX + '.removed.gff3') # DEBUG

    exon_five_term = list(exon_five_term)
    exon_three_term = list(exon_three_term)

    return exon_internal, exon_five_term, exon_three_term


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
        ranges = sorted(flatten(((l-1, 1), (r+1, -1)) for l, r in ranges))
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

    return ends_f, ends_r, genes


def get_single_exons(index_f, index_r):
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
    introns, splice_f, splice_r, SITE = get_introns(args)
    ADJ = {i:[set(), set()] for i in introns}

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
    exon_internal, _ = intron_retention.filter(args, introns, exon_internal)

    # Try and resolve ambiguous phase for internal exons #
    ######################################################
    check_adjacency(exon_internal)
    num_unresolved = len(resolve_internal_frames(exon_internal))

    # Check which terminal exons are in frame #
    ###########################################
    exon_internal, exon_five_term, exon_three_term = terminal_exons_in_phase(exon_internal, exon_l_term_dict, exon_r_term_dict)
    # do one more round to use the terminal exons to determine the phase of internal exons
    if num_unresolved > 0:
        _ = resolve_internal_frames(exon_internal + exon_five_term + exon_three_term)

    # Get all single-exons #
    ########################
    ends_f, ends_r, genes = define_gene_ends(list(introns) + exon_internal + exon_five_term + exon_three_term)
    index_f = merge_dicts(index_f, ends_f)
    index_r = merge_dicts(index_r, ends_r)
    index_f = sort_indeces(index_f)
    index_r = sort_indeces(index_r)
    exon_single = get_single_exons(index_f, index_r)

    # Output the results #
    ######################
    pprint('Outputting results...', end='')
    exon_internal = sort_by_pos(exon_internal)
    exon_five_term = sort_by_pos(exon_five_term)
    exon_three_term = sort_by_pos(exon_three_term)
    exon_final = sort_by_pos(exon_internal + exon_five_term + exon_three_term)
    exon_single = sort_by_pos(exon_single)
    print_as_gff3(exon_internal, PREFIX + '.internal.gff3')
    print_as_gff3(exon_five_term, PREFIX + '.five_term.gff3')
    print_as_gff3(exon_three_term, PREFIX + '.three_term.gff3')
    print_as_gff3(exon_single, PREFIX + '.single.gff3')
    print_as_gff3(genes, PREFIX + '.gene_boundaries.gff3', kind='gene')
    print_as_gff3(exon_final, PREFIX + '.final.gff3')
    pprint('\rOutputting results... Done!')
    pprint('  {:,} multi-exon genes, {:,} exons found'.format(len(genes), len(exon_final)))
    pprint('    {:,} internal exons'.format(len(exon_internal)))
    pprint("    {:,} 5' terminal exons".format(len(exon_five_term)))
    pprint("    {:,} 3' terminal exons".format(len(exon_three_term)))
    pprint('  {:,} single exon genes'.format(len(exon_single)))

    # Report on some specific events #
    ##################################
    # pprint('Notes:')
    # # report any internal exons that don't have an adjacent exon
    # no_l = set()
    # no_r = set()
    # no_both = set()
    # for intron, adj in ADJ.items():
    #     chrom, start, end, strand = intron
    #     l_adj = adj[0]
    #     r_adj = adj[1]
    #     if len(l_adj) < 1 and len(r_adj) < 1: # intron has no exons on either side
    #         no_both.add(intron)
    #     elif len(l_adj) < 1: # intron has no exons on upstream side
    #         no_l = no_l | r_adj
    #     elif len(r_adj) < 1: # intron has no exons on downstream side
    #         no_r = no_r | l_adj
    # no_l = sort_by_pos(no_l)
    # no_r = sort_by_pos(no_r)
    # no_both = sort_by_pos(no_both)
    # print_as_gff3(no_l, PREFIX + '.no_upstream.gff3')
    # print_as_gff3(no_r, PREFIX + '.no_downstream.gff3')
    # print_as_gff3(no_both, PREFIX + '.alone_introns.gff3', kind='intron')
    # pprint('  {:,} introns have no adjacent upstream exon'.format( len(no_l)))
    # pprint('  {:,} introns have no adjacent downstream exon'.format( len(no_r)))
    # pprint('  {:,} introns have no adjacent exons at all'.format( len(no_both)))

    # Finish up #
    #############
    pprint('Finsihed!')
