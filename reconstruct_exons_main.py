#!/home2/mattdoug/python3/bin/python3
# Author: Matt Douglas
# Last updated: 2/6/2018

# PURPOSE: Use the translation blocks encoded in the genome and splice sites
#          identified from RNA-Seq data to reconstruct protein-coding exons.
# USAGE: reconstruct_exons.py [-q] [-d] [-p exons] [-l I:100..200] -x <PREFIX of index files> -i <path to introns>

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
            print(*args, **kwargs)
    elif level == 'progress':
        if not QUIET and not DEBUG:
            print(*args, **kwargs)
    else:
        if not QUIET:
            print(*args, **kwargs)


def fr(region):
    """Return a tuple of (chromosome, start position, end position, strand) in
    a more readable format.
    """
    if len(region) < 4:
        strand = ''
    else:
        strand = ' ({} strand)'.format(region[3])

    return '{}:{:,}-{:,}{}'.format(region[0], region[1], region[2], strand)


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
                if block.num_sites() < 1:
                    continue
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
                        pprint('    internal exon found at {}'.format(fr(exon)), level='debug')
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
                # join splice site to splice site
                for l in l_sites:
                    x = [r for r in r_sites if r > l]
                    for i, j in product([l], x):
                        exon = chrom, i+1, j, '-'
                        exon_internal.add(exon)
                        FRAME[exon].add(frame)
                        pprint('    internal exon found at {}'.format(fr(exon)), level='debug')
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
    pprint("  {:,} putative internal exons".format(len(exon_internal)))
    pprint("  {:,} putative 5' terminal exons".format(five_term_c))
    pprint("  {:,} putative 3' terminal exons".format(three_term_c))

    return exon_internal, exon_l_term_dict, exon_r_term_dict


def intron_retention_events(ret_introns, introns, all_exons):
    temp_ret_introns = defaultdict(list)
    temp_introns = defaultdict(list)
    temp_exons = defaultdict(list)
    discard = set()

    pprint('Filtering intron retention events...', end='')

    for i in ret_introns:
        chrom, start, end, strand = i
        temp_ret_introns[(chrom, strand)].append((start, end))
    for i in introns:
        chrom, start, end, strand = i
        temp_introns[(chrom, strand)].append((start, end))
    for i in all_exons:
        chrom, start, end, strand = i
        temp_exons[(chrom, strand)].append((start, end))

    for region, exons in temp_exons.items():
        exons = sorted(exons, key=lambda x: (x[0], x[1]))
        ret_i = sorted(temp_ret_introns[region], key=lambda x: (x[0], x[1]))
        not_i = sorted(temp_introns[region], key=lambda x: (x[0], x[1]))
        for e in exons:
            overlaps = False
            keep = False
            for i in not_i:
                if i[1] < e[0] + 1000: # if gone past the intron, goto next exon
                    break
                if e[0] < i[0] and i[1] < e[1]:
                    overlaps = True
                    if i in ret_i:
                        keep = True
            if overlaps and not keep:
                discard.add((region[0], e[0], e[1], region[1]))

    exons = all_exons - discard

    pprint('\rFiltering intron retention events... Done!')
    pprint('  Removed {:,} exons'.format(len(discard)))

    # DEBUG:
    for i in discard:
        print('   ', fr(i))

    return exons


def build_adj_dict(exon_list):
    pprint('Building adjacency dict...', level='debug')
    for exon in exon_list:
        chrom, left, right, strand = exon
        l_adj = chrom, left - 1, strand
        r_adj = chrom, right + 1, strand
        for intron in SITE[l_adj]:
            ADJ[intron][1].add(exon)
            pprint('  exon is right of intron {}'.format(fr(intron)), level='debug')
        for intron in SITE[r_adj]:
            ADJ[intron][0].add(exon)
            pprint('  exon is left of intron {}'.format(fr(intron)), level='debug')
    pprint('Building adjacency dict... Done!', level='debug')


def resolve_internal_frames(exon_internal, max_iterations=20):
    unresolved = []
    count = 0
    total = len(exon_internal)

    pprint('', level='debug')
    pprint('Identifying the correct phase for internal exons', end='')

    for n in range(max_iterations):
        pprint('\nStarting iteration {}:'.format(n+1), level='debug')
        new_FRAME = defaultdict(set)
        num_change = 0
        for exon in exon_internal:
            if len(FRAME[exon]) > 1:
                pprint('Iteration {}) Block {} has more than one frame ({}), skipping'.format(n+1, fr(exon), FRAME[exon]), level='debug')
                continue
            chrom, left, right, strand = exon
            (frame,) = FRAME[exon]
            upstream_introns = SITE[(chrom, left-1, strand)]
            downstream_introns = SITE[(chrom, right+1, strand)]
            pprint('Iteration {}) Using exon {} (frame {}) to resolve adjacent reading frames'.format(n+1, fr(exon), frame), level='debug')
            ####################################################################
            for intron in upstream_introns:
                pprint('  Looking for exons adjacent to upstream intron {}...'.format(fr(intron)), level='debug')
                for adj_exon in ADJ[intron][0]:
                    adj_frames = FRAME[adj_exon]
                    if len(adj_frames) != 1:
                        pprint('  Determining frame for exon {}, frames: {}'.format(fr(adj_exon), adj_frames), level='debug')
                        new_frame = shift_frame(frame, intron, strand+'<')
                        new_FRAME[adj_exon].add(new_frame)
            ####################################################################
            for intron in downstream_introns:
                pprint('  Looking for exons adjacent to downstream intron {}...'.format(fr(intron)), level='debug')
                for adj_exon in ADJ[intron][1]:
                    adj_frames = FRAME[adj_exon]
                    if len(adj_frames) != 1:
                        pprint('  Determining frame for exon {}, frames {}:'.format(fr(adj_exon), adj_frames), level='debug')
                        new_frame = shift_frame(frame, intron, strand+'>')
                        new_FRAME[adj_exon].add(new_frame)

        # update reading frames
        pprint('\nSummary of frame changes:', level='debug')
        for exon, new_frames in new_FRAME.items():
            old_frames = FRAME[exon]
            if old_frames != new_frames:
                pprint('  Block {} frames {} => {}'.format(exon, old_frames, new_frames), level='debug')
                FRAME[exon] = new_frames
                num_change += 1
        if num_change < 1:
            pprint('  None. Stopping at iteration {}'.format(n+1), level='debug')
            break

    # check for exons that still have more than one valid frame
    pprint('\nExons with unknown phase:', level='debug')
    for exon, frames in FRAME.items():
        if len(frames) > 1:
            unresolved.append(exon)
            pprint('  {}, frames: {}'.format(fr(exon), frames), level='debug')

    # DEBUG: Output internal exons with ambiguous frame
    unresolved = sort_by_pos(unresolved)
    print_as_gff3(unresolved, PREFIX + '.internal.multi_frame.gff3')

    pprint('\rIdentifying the correct phase for internal exons... Done!')
    pprint('  Resolved {:,} cases where an internal exons had >1 possible phase'.format(num_change))
    pprint('  The phase of {:,} ({}) internal exons could not be determined'
           .format(len(unresolved), percent(len(unresolved), len(exon_internal))))


def trim_false_internal_exons(exon_internal, exon_l_term_dict, exon_r_term_dict, max_iterations=3):
    utr_introns = set()
    removed_set = set()

    pprint('', level='debug')
    pprint('Trimming false internal exons', end='')

    for n in range(max_iterations):
        pprint('\nStarting iteration {}:'.format(n+1), level='debug')
        to_remove_exons = set()
        for exon in exon_internal:
            pprint('Checking internal exon {}, frame(s) {}...'.format(fr(exon), FRAME[exon]), level='debug')
            chrom, start, end, strand = exon
            upstream_introns = SITE[(chrom, start-1, strand)]
            downstream_introns = SITE[(chrom, end+1, strand)]
            validate_l = False
            validate_r = False
            ####################################################################
            for intron in upstream_introns:
                pprint('  Looking for exons adjacent to upstream intron {}...'.format(fr(intron)), level='debug')
                adj_internal_exons = ADJ[intron][0]
                adj_term_exons = exon_l_term_dict[intron[0], intron[1]-1, intron[3]]
                expect_frames = set([shift_frame(f, intron, strand+'<') for f in FRAME[exon]])
                pprint('    Expected frames =', expect_frames, level='debug')
                for adj_exon in adj_internal_exons:
                    if expect_frames & FRAME[adj_exon]:
                        pprint('    Internal exon {} is in frame (expect: {}, found: {})'.format(fr(adj_exon), expect_frames, FRAME[adj_exon]), level='debug')
                        validate_l = True
                        continue
                for adj_exon in adj_term_exons:
                    if expect_frames & FRAME[adj_exon]:
                        pprint('    Terminal exon {} is upstream, frame(s) {}'.format(fr(adj_exon), FRAME[adj_exon]), level='debug')
                        validate_l = True
                        continue
            ####################################################################
            for intron in downstream_introns:
                pprint('  Looking for exons adjacent to downstream intron {}...'.format(fr(intron)), level='debug')
                adj_internal_exons = ADJ[intron][1]
                adj_term_exons = exon_r_term_dict[intron[0], intron[2]+1, intron[3]]
                expect_frames = set([shift_frame(f, intron, strand+'>') for f in FRAME[exon]])
                for adj_exon in adj_internal_exons:
                    if expect_frames & FRAME[adj_exon]:
                        pprint('    Internal exon {} is in frame (expect: {}, found: {})'.format(fr(adj_exon), expect_frames, FRAME[adj_exon]), level='debug')
                        validate_r = True
                        continue
                for adj_exon in adj_term_exons:
                    if expect_frames & FRAME[adj_exon]:
                        pprint('    Terminal exon {} is downstream, frame(s) {}'.format(fr(adj_exon), FRAME[adj_exon]), level='debug')
                        validate_r = True
                        continue
            ####################################################################
            if validate_l and validate_r:
                pprint('  Valid!', level='debug')
                continue
            if validate_l:
                pprint('  No downstream exon! Removing.', level='debug')
                to_remove_exons.add(exon)
                for i in downstream_introns:
                    utr_introns.add(i)
            if validate_r:
                pprint('  No upstream exon! Removing.', level='debug')
                to_remove_exons.add(exon)
                for i in upstream_introns:
                    utr_introns.add(i)
            for intron in upstream_introns:
                pprint('    Removing exon from ADJ of intron {} (right)'.format(fr(intron)), level='debug')
                pprint('     ', ADJ[intron], level='debug')
                try:
                    ADJ[intron][1].remove(exon)
                except KeyError:
                    pass
                pprint('     ', ADJ[intron], level='debug')
            for intron in downstream_introns: # remove the exons from the adjaceny dict
                pprint('    Removing exon from ADJ of intron {} (left)'.format(fr(intron)), level='debug')
                pprint('     ', ADJ[intron], level='debug')
                try:
                    ADJ[intron][0].remove(exon)
                except KeyError:
                    pass
                pprint('     ', ADJ[intron], level='debug')

        if len(to_remove_exons) > 0:
            pprint('\nRemoved {:,} internal exons'.format(len(to_remove_exons)), level='debug')
            exon_internal = list(set(exon_internal) - to_remove_exons)
            removed_set = removed_set.intersection(to_remove_exons)
        else:
            pprint('\nNo more exons to remove. Stopping iterations.', level='debug')
            break

    pprint('\rTrimming back false-ends... Done!')
    pprint('  Removed {:,} false internal exons'.format(len(removed_set)))
    print_as_gff3(removed_set, PREFIX + '.trimmed_int_exons.gff3')

    return exon_internal, utr_introns


def terminal_exons_in_phase(exon_l_term_dict, exon_r_term_dict, utr_introns):
    # quick function to calculate the length of an exon
    length = lambda x: x[2] - x[1] + 1

    new_l_term = set()
    new_r_term = set()
    new_ADJ = defaultdict(list)
    five_term_c = 0
    three_term_c = 0
    exon_five_term = []
    exon_three_term = []
    # for keeping track of progress
    count = 0
    total = len(ADJ)

    single_intron = [] # DEBUG

    pprint('', level='debug')
    pprint('Checking which terminal exons are in phase', end='')

    for intron, adj in ADJ.items():
        pprint('\rChecking which terminal exons are in phase...', percent(count, total), end='', level='progress')
        count += 1
        chrom, l, r, strand = intron
        l_adj, r_adj = adj
        new_ADJ[intron] = [ set(), set() ]
        ########################################################################
        if intron in utr_introns:
            continue
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
                for f in FRAME[exon]:
                    if f not in l_compare:
                        l_compare[f] = set()
                    l_compare[f].add(exon)
            for exon in exon_r_term_dict[(chrom, r+1, strand)]:
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
            if len(results) > 0:
                final = max(results, key=lambda x: (x[0][2]-x[0][1] + x[1][2]-x[1][1]))
                new_l_term.add(final[0])
                new_r_term.add(final[1])
                five_term_c += 1
                three_term_c += 1
                pprint('  Longest ORF is:', level='debug')
                # add exons to the adjacency dictionary
                new_ADJ[intron][0].add(final[0])
                new_ADJ[intron][1].add(final[1])
                pprint('      {} is left of intron {}'.format(fr(final[0]), fr(intron)), level='debug')
                pprint('      {} is right of intron {}'.format(fr(final[1]), fr(intron)), level='debug')
                # DEBUG: print out any single-intron genes
                single_intron.append(final[0])
                single_intron.append(final[1])
            else:
                pprint('  No valid ORF found!', level='debug')
        ########################################################################
        # if no internal exon is upstream of the intron...
        elif len(l_adj) < 1:
            pprint('Looking upstream of intron at {}'.format(fr(intron)), level='debug')
            l_compare = {}
            expected = set()
            # determine what the expected frame(s) should be based on adjacent
            # internal exons
            for exon in r_adj:
                for f in FRAME[exon]:
                    new_f = shift_frame(f, intron, strand+'<')
                    pprint('  Downstream exon {}, frame {}'.format(exon, f), level='debug')
                    expected.add(new_f)
            pprint('  Expect adjacent frames to be one of: {}'.format(expected), level='debug')
            # only keep the longest exon(s) that match the expected frame(s)
            for exon in exon_l_term_dict[(chrom, l-1, strand)]:
                for f in FRAME[exon]:
                    if f not in l_compare:
                        l_compare[f] = set()
                    l_compare[f].add(exon)
            for f in expected:
                if f in l_compare:
                    match = max(l_compare[f], key=length)
                    new_l_term.add(match)
                    pprint('  Blocks {} with frame {} match'.format(', '.join([fr(i) for i in l_compare[f]]), f), level='debug')
                    pprint('  Longest is {}'.format(fr(match)), level='debug')
                    # add exons to the adjacency dictionary
                    new_ADJ[intron][0].add(match)
                    pprint('    {} is left of intron {}'.format(fr(match), fr(intron)), level='debug')
                    # count 5' and 3' terminal exons
                    if strand == '+':
                        five_term_c += 1
                    elif strand == '-':
                        three_term_c += 1
        ########################################################################
        # if no internal exon is downstream of the intron...
        elif len(r_adj) < 1:
            pprint('Looking downstream of intron at {}'.format(fr(intron)), level='debug')
            r_compare = {}
            expected = set()
            # determine what the expected frame(s) should be based on adjacent
            # internal exons
            for exon in l_adj:
                for f in FRAME[exon]:
                    new_f = shift_frame(f, intron, strand+'>')
                    expected.add(new_f)
            pprint('  Expect adjacent frames to be one of: {}'.format(expected), level='debug')
            # only keep the longest exon(s) that match the expected frame(s)
            for exon in exon_r_term_dict[(chrom, r+1, strand)]:
                for f in FRAME[exon]:
                    if f not in r_compare:
                        r_compare[f] = set()
                    r_compare[f].add(exon)
            for f in expected:
                if f in r_compare:
                    match = max(r_compare[f], key=length)
                    new_r_term.add(match)
                    pprint('  Blocks {} with frame {} match'.format(', '.join([fr(i) for i in r_compare[f]]), f), level='debug')
                    pprint('  Longest is {}'.format(fr(match)), level='debug')
                    # add exons to the adjacency dictionary
                    new_ADJ[intron][1].add(match)
                    pprint('    {} is right of intron {}'.format(fr(match), fr(intron)), level='debug')
                    # count 5' and 3' terminal exons
                    if strand == '+':
                        three_term_c += 1
                    elif strand == '-':
                        five_term_c += 1
        ########################################################################
        # otherwise, continue...
        else:
            pprint('Skipping intron at {}, is flanked by exons'.format(fr(intron)), level='debug')
            continue

    # update the adjacency dictionary
    pprint('Updating adjacency dictionary...', level='debug')
    for intron, adj in new_ADJ.items():
        for l in adj[0]:
            ADJ[intron][0].add(l)
        for r in adj[1]:
            ADJ[intron][1].add(r)

    # seperate out 5' terminal exons from 3' terminal exons
    for i in new_l_term:
        if i[3] == '+':
            exon_five_term.append(i)
        else:
            exon_three_term.append(i)
    for i in new_r_term:
        if i[3] == '+':
            exon_three_term.append(i)
        else:
            exon_five_term.append(i)

    pprint('\rChecking which terminal exons are in phase... Done!  ')
    pprint("  {:,} 5' terminal exons".format(five_term_c))
    pprint("  {:,} 3' terminal exons".format(three_term_c))

    print_as_gff3(single_intron, PREFIX + '.double.gff3') # DEBUG

    return exon_five_term, exon_three_term


def define_gene_ends(features):
    """Return gene boundaries (i.e. contigous regions of internal exons,
    terminal exons, and introns.
    """
    flatten = chain.from_iterable
    temp1   = defaultdict(list)
    temp2   = defaultdict(list)
    ends_f  = {}
    ends_r  = {}
    genes   = []

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

    # for DEBUGging: report gene boundaries
    if DEBUG:
        for d, s in ((ends_f, '+'), (ends_r, '-')):
            for chrom, frames in d.items():
                for pos, kind in frames[0]:
                    if kind == 2:
                        pprint('left side of gene at: {}:{:,} ({} strand)'.format(chrom, pos, s), level='debug')
                    elif kind == 3:
                        pprint('right side of gene at: {}:{:,} ({} strand)'.format(chrom, pos, s), level='debug')

    # DEBUG: Print out putative gene boundaries
    genes = sort_by_pos(genes)
    print_as_gff3(genes, PREFIX + '.gene_boundaries.gff3', kind='gene')

    pprint('\rDefining putative gene boundaries... Done!')
    pprint('  {:,} multi-exon genes found'.format( sum([len(i) for i in temp2.values()])))

    return ends_f, ends_r, genes


def get_single_exons(index_f, index_r):
    exon_single = set()
    count       = 0
    total       = len(index_f)*6 #calculate the number of positions to check

    pprint('', level='debug')
    pprint('Checking intergenic regions for single exons', end='')

    # scan the + strand of the genome
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
                elif kind == 1:
                    if len(track_list) > 0:
                        for t, x in track_list:
                            # if we're tracking from a start codon...
                            if x == 0:
                                exon = chrom, t, pos-1, '+'
                                exon_single.add(exon)
                                FRAME[exon].add(frame)
                                pprint('  single-exon exon found at {}'.format(fr(exon)), level='debug')
                    # reset
                    track_list = []
                # When a left splice site is reached...
                elif kind == 2:
                    track_list = []
                    intergenic = False
                # When a right splice site is reached...
                elif kind == 3:
                    track_list = []
                    intergenic = True

    # scan the - strand of the genome
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
                            pprint('  single-exon exon found at {}'.format(fr(exon)), level='debug')
                # When a stop codon is reached...
                elif kind == 1:
                    if intergenic:
                        track_list = [ (pos, 0) ]
                # When a left splice site is reached...
                elif kind == 2:
                    track_list = []
                    intergenic = False
                # When a right splice site is reached...
                elif kind == 3:
                    track_list = []
                    intergenic = True

    pprint('\rChecking intergenic regions for single exons... Done!    ')
    pprint('  {:,} single exons'.format(len(exon_single)))

    return exon_single


############
# Run loop #
############
def reconstruct_exons(index_f, index_r, introns, ret_introns, splice_f, splice_r, sd, args):
    global QUIET
    global DEBUG
    global PREFIX
    global ADJ
    global FRAME
    global NOTE
    global SITE
    QUIET = args.quiet
    DEBUG = args.debug
    PREFIX = args.prefix
    ADJ = {i:[set(), set()] for i in introns}
    FRAME = defaultdict(set)
    NOTE = defaultdict(set)
    SITE = sd

    ret_introns = [('II', 8053919, 8053975, '+')]

    pprint('Starting exon reconstruction', level='debug')

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
    exon_internal = intron_retention_events(ret_introns, introns, exon_internal) # TODO

    # Try and resolve ambiguous phase for internal exons #
    ######################################################
    build_adj_dict(exon_internal)
    resolve_internal_frames(exon_internal)

    # Check which terminal exons are in frame #
    ###########################################
    exon_internal, utr_introns = trim_false_internal_exons(exon_internal, exon_l_term_dict, exon_r_term_dict)
    exon_five_term, exon_three_term = terminal_exons_in_phase(exon_l_term_dict, exon_r_term_dict, utr_introns)
    # do one more round to use the terminal exons to determine the phase of internal exons
    resolve_internal_frames(list(exon_internal) + list(exon_five_term) + list(exon_three_term))

    # Get all single-exons #
    ########################
    ends_f, ends_r, genes = define_gene_ends(list(introns) + list(exon_internal) + exon_five_term + exon_three_term)
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
    exon_single = sort_by_pos(exon_single)
    print_as_gff3(exon_internal, PREFIX + '.internal.gff3')
    print_as_gff3(exon_five_term, PREFIX + '.five_term.gff3')
    print_as_gff3(exon_three_term, PREFIX + '.three_term.gff3')
    print_as_gff3(exon_single, PREFIX + '.single.gff3', NOTE)
    exon_all = sort_by_pos(exon_internal + exon_five_term + exon_three_term)
    print_as_gff3(exon_all, PREFIX + '.final.gff3')
    pprint('\rOutputting results... Done!')
    pprint('  {:,} internal exons'.format(len(exon_internal)))
    pprint("  {:,} 5' terminal exons".format(len(exon_five_term)))
    pprint("  {:,} 3' terminal exons".format(len(exon_three_term)))
    pprint('  {:,} single exon genes'.format(len(exon_single)))
    pprint('  {:,} exons total (not counting single-exon genes)'.format(len(exon_all)))

    # Report on some specific events #
    ##################################
    pprint('Notes:')
    # report any internal exons that don't have an adjacent exon
    no_l = set()
    no_r = set()
    no_both = set()
    for intron, adj in ADJ.items():
        chrom, start, end, strand = intron
        l_adj = adj[0]
        r_adj = adj[1]
        if len(l_adj) < 1 and len(r_adj) < 1: # intron has no exons on either side
            no_both.add(intron)
        elif len(l_adj) < 1: # intron has no exons on upstream side
            no_l = no_l | r_adj
        elif len(r_adj) < 1: # intron has no exons on downstream side
            no_r = no_r | l_adj
    no_l = sort_by_pos(no_l)
    no_r = sort_by_pos(no_r)
    no_both = sort_by_pos(no_both)
    print_as_gff3(no_l, PREFIX + '.no_upstream.gff3')
    print_as_gff3(no_r, PREFIX + '.no_downstream.gff3')
    print_as_gff3(no_both, PREFIX + '.alone_introns.gff3', kind='intron')
    pprint('  {:,} exons have no adjacent upstream exon'.format( len(no_l)))
    pprint('  {:,} exons have no adjacent downstream exon'.format( len(no_r)))
    pprint('  {:,} introns have no adjacent exons'.format( len(no_both)))

    # Finish up #
    #############
    pprint('Finsihed!')
