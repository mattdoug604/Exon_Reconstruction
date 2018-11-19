#!/usr/bin/python3
#
# An algorithm implemented in Python to identify putative protein-coding exons
# from a reference genome and a set of splice sites.
# Copyright (C) 2018 Matthew Douglas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# CONTACT INFORMATION:
# By email at 'mattdoug604@gmail.com'

from collections import defaultdict
from itertools import chain, product
import argparse, logging, re, sys
import pysam
from parse_genome_index import parse_index
from parse_intron_gff import get_introns
import filter_intron_retention

VERSION="0.18.3"

log_format = '%(message)s'
logging.basicConfig(format=log_format, level=logging.INFO)
log = logging.getLogger(__name__)

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

    def has_splice_site(self):
        return any((self.l_sites, self.r_sites))


#####################
# Utility functions #
#####################
def parse_region_str(string):
    """Convert a string in the format: 'I:1000..2000' or 'I_1000_2000' or
       'I:1000-2000' to a tuple.
    """
    if string is None:
        return None

    t = string.strip().replace(',', '').replace('..', '-')
    t = re.split(':|-|_', t)
    seq, start, end = t[0], int(t[1]), int(t[2])
    region = seq, start, end

    return region


def parse_commandline_args():
    parser = argparse.ArgumentParser(
        description='Reconstruct exons from RNA-Seq data.')
    requiredNamed = parser.add_argument_group('required input')
    # required input ###########################################################
    requiredNamed.add_argument(
        dest='index',
        type=str,
        help='prefix of the genome index files')
    requiredNamed.add_argument(
        dest='introns',
        type=str,
        help='intron in GFF3 format')
    requiredNamed.add_argument(
        dest='align',
        type=str,
        help='aligned reads in SAM or BAM format')
    # optional input ###########################################################
    parser.add_argument(
        '-o',
        dest='output',
        type=str,
        nargs='?',
        default='exons',
        help='prefix for the output files')
    parser.add_argument(
        '-r',
        dest='region',
        type=parse_region_str,
        nargs='?',
        help='limit search to genomic coordinates (e.g. I:1000-2000)')
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        default=False,
        help='do not print any information [incompatible w/ -d]')
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        default=False,
        help='print debug information')
    ############################################################################
    args = parser.parse_args()

    if not args.quiet:
        log.info('Index prefix  = {}'.format(args.index))
        log.info('Intron file   = {}'.format(args.introns))
        log.info('Aligned reads = {}'.format(args.align))
        log.info('Output prefix = {}'.format(args.output))
        if args.region:
            log.info('Search region = {}:{:,}-{:,}'.format(*args.region))
        log.info('')

    return args


def percent(val, total):
    if total == 0:
        return '0%'
    per = val * 100 / total
    per = '{0:.1f}%'.format(per)

    return per


def merge_dicts(x, y):
    """Merge to dicts (x, y) with the structure: key:[[], [], []] into a
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


def frame_shift(frame, intron, direction):
    """Calculate the frame a given exon uses based on the frame of the adjacent
    exon and the lenght of the intervening intron.
    """
    i_len = intron[2]-intron[1]+1
    i_mod = i_len % 3

    if direction in ('+<', '->'):
        new_frame = (frame - i_mod) % 3
    elif direction in ('+>', '-<'):
        new_frame = (frame + i_mod) % 3

    log.debug('    Calculating frame shift:')
    log.debug('      Direction is {}'.format(direction))
    log.debug('      Exon is in frame {}'.format(frame))
    log.debug('      Intron {} is {}bp long, modulo is {}'.format(intron, i_len, i_mod))
    log.debug('      New frame is {}'.format(new_frame))

    return new_frame


def sort_indeces(index_dict):
    """Take a dictionary of tuples, remove duplicates and sort by value."""
    log.debug('Removing duplicates and sorting by position...')

    for chrom, frames in index_dict.items():
        for n, f in enumerate(frames):
            index_dict[chrom][n] = sorted(set(f))

    log.debug('Removing duplicates and sorting by position... Done!')

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
    left_introns = SITE_R[(chrom, left-1, strand)]
    right_introns = SITE_L[(chrom, right+1, strand)]

    if len(left_introns) * len(right_introns) > 0:
        l_sum = sum(INTRON[i] for i in left_introns)
        r_sum = sum(INTRON[i] for i in right_introns)
        return min(l_sum, r_sum)
    elif len(left_introns) > 0:
        return sum(INTRON[i] for i in left_introns)
    elif len(right_introns) > 0:
        return sum(INTRON[i] for i in right_introns)
    else:
        return 0
        # raise KeyError('No adjacent intron(s) found for exon: {}:{}-{}{}'.format(*exon))
        # sys.exit(1)


def output_as_gff3(features, out_path, kind='CDS', mode='w'):
    """Take a list of features as a tuple (chromosome, start, end, strand),
    and print them as GFF3 formatted entries.
    """
    log.debug('Printing results to "{}"'.format(out_path))

    with open(out_path, mode) as f:
        if 'w' in mode: # print the header only if we're opening a new file
            print('##gff-version 3', file=f)
            print('##{} version_{}'.format(sys.argv[0], VERSION), file=f)
        for n, feat in enumerate(features):
            chrom, left, right, strand = feat
            if kind == 'CDS':
                score = calc_score(feat)
            elif kind == 'intron':
                score = INTRON[feat]
            else:
                score = '.'
            frame = ','.join(map(str, FRAME[feat]))
            notes = [';Frame='+frame if frame != '' else ''][0]
            line = chrom, '.', kind, left, right, score, strand, '.', 'ID={}{}{}'.format(kind, n+1, notes)
            print(*line, sep='\t', file=f)


#######################
# Exon Reconstruction #
#######################
def get_translation_blocks(index_f, index_r):
    """Return two dicts of TranslationBlock objects. A 'translation block' is
    defined as the region from the first base of the first codon after a stop
    codon to the last base of the next stop codon in the same reading frame.

    INPUT: Two dicts corresponding to the indexed positions of all stop codons,
    putative start codons, and splice sites on the forward and reverse strands
    of the genome, respectively. Keys in each dict are the chromosomes. The
    values are 3 lists corresponding to each reading frame:
    index_f = { chromsome: [ [frame0], [frame1], [frame2] ] }
    index_r = { chromsome: [ [frame0], [frame1], [frame2] ] }

    RETURN: Two dicts (blocks_f, blocks_r), structured the same as the input
    where each of the three lists (for each key) is a list of TranslationBlock
    objects in that reading frame.
    """
    sort_by_type = lambda x: sorted(x, key=lambda y: (y[0], y[1]))
    blocks_f = {}
    blocks_r = {}
    empty_blocks_f = {}
    empty_blocks_r = {}
    total_blk = 0
    empty_blk = 0
    count = 0
    total = len(index_f)*6 #calculate the number of positions to check

    log.info('\nIdentifying all translation blocks')

    for chrom, frames in index_f.items():
        blocks_f[chrom] = [ [], [], [] ]
        empty_blocks_f[chrom] = [ [], [], [] ]
        for frame, i in enumerate(frames):
            log.debug('\nScanning + strand in the > direction, reading frame {}'.format(frame))
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
                    if block.has_splice_site():
                        blocks_f[chrom][frame].append(block)
                    else:
                        empty_blocks_f[chrom][frame].append(block)
                        empty_blk += 1
                    log.debug('  {}:{:,}-{:,} (+ strand)'.format(chrom, block.start, block.end))
                    if len(block.l_sites) > 0:
                        log.debug("    5' splicesites: {}".format(' '.join(map(str, block.l_sites))))
                    if len(block.r_sites) > 0:
                        log.debug("    3' splicesites: {}".format(' '.join(map(str, block.r_sites))))
                    if len(block.s_sites) > 0:
                        log.debug("    start sites: {}".format(' '.join(map(str, block.s_sites))))
                    # reset
                    total_blk += 1
                    last = pos
                    block = TranslationBlock()
    ############################################################################
    for chrom, frames in index_r.items():
        blocks_r[chrom] = [ [], [], [] ]
        empty_blocks_r[chrom] = [ [], [], [] ]
        for frame, i in enumerate(frames):
            log.debug('\nScanning - strand in the > direction, reading frame {}'.format(frame))
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
                    if block.has_splice_site():
                        blocks_r[chrom][frame].append(block)
                    else:
                        empty_blocks_r[chrom][frame].append(block)
                        empty_blk += 1
                    log.debug('  {}:{:,}-{:,} (- strand)'.format(chrom, block.start, block.end))
                    if len(block.l_sites) > 0:
                        log.debug("    5' splicesites: {}".format(' '.join(map(str, block.l_sites))))
                    if len(block.r_sites) > 0:
                        log.debug("    3' splicesites: {}".format(' '.join(map(str, block.r_sites))))
                    if len(block.s_sites) > 0:
                        log.debug("    start sites: {}".format(' '.join(map(str, block.s_sites))))
                    # reset
                    total_blk += 1
                    last = pos
                    block = TranslationBlock()

    log.info('  {:,} translation blocks'.format(total_blk))
    log.info('  {:,} translation blocks contain one or more splice sites'.format(total_blk - empty_blk))

    # DEBUG: output all the translation blocks
    temp_f = merge_dicts(blocks_f, empty_blocks_f)
    temp_r = merge_dicts(blocks_r, empty_blocks_r)
    out_path = '{}.translation_blocks.gff3'.format(PREFIX)
    output_as_gff3([], out_path, mode='w') # create a new file
    for d, strand, suf in ((temp_f, '+', 'plus'), (temp_r, '-', 'minus')):
        temp = [ [], [], [] ]
        for chrom, frames in d.items():
            for frame, blocks in enumerate(frames):
                temp[frame] += [(chrom, i.start, i.end, strand) for i in blocks]
        for frame, blocks in enumerate(temp):
            blocks = sort_by_pos(blocks)
            output_as_gff3(blocks, out_path, kind=suf+str(frame), mode='a')

    return blocks_f, blocks_r


def get_putative_exons(blocks_f, blocks_r):
    """Search each translation block and return the genomic coordinates of all
    possible internal and terminal exons, which reading frame(s) each use, and
    which intron(s) the exons are adjacent to.

    INPUT: Two dicts of TranslationBlock objects on forward and reverse strands
    of the genome, respectively. Keys for each dict are the chromosomes.
    Values are three lists corresponding to each reading frame:
    blocks_f = { chromsome: [ [frame0], [frame1], [frame2] ] }
    blocks_r = { chromsome: [ [frame0], [frame1], [frame2] ] }

    RETURN:
    1) A set of "internal exons" as tuples in the format: (chromosome, start,
       end, strand).
    2) Two dicts of terminal exons (one for 5' terminal exons, one for 3'
       terminal) where the key is the adjacent intron(s) and the value is a set
       of exon positions as tuples (chromosome, start, end, strand).
    """
    exon_intrnl = set()
    exon_l_term_dict = defaultdict(set)
    exon_r_term_dict = defaultdict(set)
    interl_c = 0
    l_term_c = 0
    r_term_c = 0
    count = 0
    total = len(blocks_f)*6  #calculate the number of positions to check

    log.debug('')
    log.info('\nIdentifying putative exons positions')

    # scan the + strand of the genome
    for chrom, frames in blocks_f.items():
        for frame, i in enumerate(frames):
            log.debug('\nScanning + strand, reading frame {}'.format(frame))
            count += 1
            for block in i:
                l_sites = sorted(block.l_sites)
                r_sites = sorted(block.r_sites)
                s_sites = sorted(block.s_sites)
                log.debug('  Checking block: {}:{:,}-{:,} (+ strand)'.format(chrom, block.start, block.end))
                if len(l_sites) > 0: log.debug('  l_sites:'.format(l_sites))
                if len(r_sites) > 0: log.debug('  r_sites:'.format(r_sites))
                if len(s_sites) > 0: log.debug('  s_sites:'.format(s_sites))
                # join splice site to splice site
                for l in l_sites:
                    x = [r for r in r_sites if r > l]
                    for i, j in product([l], x):
                        exon = chrom, i, j-1, '+'
                        exon_intrnl.add(exon)
                        FRAME[exon].add(frame)
                        log.debug('    internal exon found at {}'.format(exon))
                # join start codons to 5' splice sites
                for s in s_sites:
                    x = [r for r in r_sites if r > s]
                    for i, j in product([s], x):
                        exon = chrom, i, j-1, '+'
                        r_adj = chrom, j-1, '+'
                        exon_l_term_dict[r_adj].add(exon)
                        FRAME[exon].add(frame)
                        l_term_c += 1
                        log.debug("    5' terminal exon found at {}".format(exon))
                # join 3' splice sites to stop codons
                for l in [s for s in l_sites if s < block.end]:
                    exon = chrom, l, block.end, '+'
                    l_adj = chrom, l, '+'
                    exon_r_term_dict[l_adj].add(exon)
                    FRAME[exon].add(frame)
                    r_term_c += 1
                    log.debug("    3' terminal exon found at {}".format(exon))
    ############################################################################
    # scan the - strand of the genome
    for chrom, frames in blocks_r.items():
        for frame, i in enumerate(frames):
            log.debug('\nScanning - strand, reading frame {}'.format(frame))
            count += 1
            for block in i:
                l_sites = sorted(block.l_sites)
                r_sites = sorted(block.r_sites)
                s_sites = sorted(block.s_sites)
                log.debug('  Checking block: {}:{:,}-{:,} (- strand)'.format(chrom, block.start, block.end))
                if len(l_sites) > 0: log.debug('  l_sites:'.format(l_sites))
                if len(r_sites) > 0: log.debug('  r_sites:'.format(r_sites))
                if len(s_sites) > 0: log.debug('  s_sites:'.format(s_sites))
                # join splice site to splice site
                for l in l_sites:
                    x = [r for r in r_sites if r > l]
                    for i, j in product([l], x):
                        exon = chrom, i+1, j, '-'
                        exon_intrnl.add(exon)
                        FRAME[exon].add(frame)
                        log.debug('    internal exon found at {}'.format(exon))
                # join stop codon to 3' splice sites
                for r in [s for s in r_sites if s > block.start]:
                    exon = chrom, block.start, r, '-'
                    r_adj = chrom, r, '-'
                    exon_l_term_dict[r_adj].add(exon)
                    FRAME[exon].add(frame)
                    r_term_c += 1
                    log.debug("    5' terminal exon found at {}".format(exon))
                # join 5' splice sites to start codons
                for l in l_sites:
                    x = [s for s in s_sites if s > l]
                    for i, j in product([l], x):
                        exon = chrom, i+1, j, '-'
                        l_adj = chrom, i+1, '-'
                        exon_r_term_dict[l_adj].add(exon)
                        FRAME[exon].add(frame)
                        l_term_c += 1
                        log.debug("    3' terminal exon found at {}".format(exon))

    log.info("  {:,} putative internal exons".format(len(exon_intrnl)))
    log.info("  {:,} putative 5' terminal exons".format(l_term_c))
    log.info("  {:,} putative 3' terminal exons".format(r_term_c))

    exon_intrnl = list(exon_intrnl)

    return exon_intrnl, exon_l_term_dict, exon_r_term_dict


def check_adjacency(exon_list):
    """Given a list of exons and a dictionary of splice sites and their
    corresponding introns, return a dictionary where each key is an intron, and
    the value is two lists of exons that are adjacent to the intron - the first
    list is exons at the 5' end, the second list is exons at the 3' end."""
    global ADJ

    log.debug('\nBuilding adjacency dict...')

    for exon in exon_list:
        chrom, left, right, strand = exon
        l_adj = chrom, left - 1, strand
        r_adj = chrom, right + 1, strand
        for intron in SITE_R[l_adj]:
            ADJ[intron][1].add(exon)
            log.debug('  exon {} is right of intron {}'.format(exon, intron))
        for intron in SITE_L[r_adj]:
            ADJ[intron][0].add(exon)
            log.debug('  exon {} is left of intron {}'.format(exon, intron))

    log.debug('Building adjacency dict... Done!')


def resolve_internal_frames(exon_list, max_iter=20):
    """Sometimes a putative internal exon may appear in more than one reading
    frame. Here we attempt to resolve any abiguous frames by iteratively
    checking the frame of adjacent exons, and using those to calulate what the
    frame should be."""
    unresolved = []
    count = 0
    total = len(exon_list)
    total_change = 0

    log.debug('')
    log.info('\nIdentifying the correct phase for internal exons')

    for n in range(max_iter):
        log.debug('\nStarting iteration {}:'.format(n+1))
        new_FRAME = defaultdict(set)
        num_change = 0
        for exon in exon_list:
            if len(FRAME[exon]) > 1:
                log.debug('Iteration {}) Block {} has more than one frame ({}), skipping'.format(n+1, exon, FRAME[exon]))
                continue
            chrom, left, right, strand = exon
            (frame,) = FRAME[exon]
            left_introns = SITE_R[(chrom, left-1, strand)]
            right_introns = SITE_L[(chrom, right+1, strand)]
            log.debug('Iteration {}) Using exon {} (frame {}) to resolve adjacent reading frames'.format(n+1, exon, frame))
            ####################################################################
            for intron in left_introns:
                log.debug('  Looking for exons adjacent to upstream intron {}...'.format(intron))
                for adj_exon in ADJ[intron][0]:
                    adj_frames = FRAME[adj_exon]
                    if len(adj_frames) != 1:
                        log.debug('  Determining frame for exon {}, frames: {}'.format(adj_exon, adj_frames))
                        new_frame = frame_shift(frame, intron, strand+'<')
                        if new_frame in adj_frames: # avoid conflicts where incompatible adjacent exons assign the wrong frame
                            new_FRAME[adj_exon].add(new_frame)
                        else:
                            log.debug('    Conflict: frame {} is not a valid choice {}'.format(new_frame, adj_frames))
            ####################################################################
            for intron in right_introns:
                log.debug('  Looking for exons adjacent to downstream intron {}...'.format(intron))
                for adj_exon in ADJ[intron][1]:
                    adj_frames = FRAME[adj_exon]
                    if len(adj_frames) != 1:
                        log.debug('  Determining frame for exon {}, frames {}:'.format(adj_exon, adj_frames))
                        new_frame = frame_shift(frame, intron, strand+'>')
                        if new_frame in FRAME[adj_exon]:
                            new_FRAME[adj_exon].add(new_frame)
                        else:
                            log.debug('    Conflict: frame {} is not a valid choice ({})'.format(new_frame, adj_frames))

        # update reading frames
        log.debug('\nSummary of frame changes:')
        for exon, new_frames in new_FRAME.items():
            old_frames = FRAME[exon]
            if old_frames != new_frames:
                log.debug('  Block {} frames {} => {}'.format(exon, old_frames, new_frames))
                FRAME[exon] = new_frames
                num_change += 1
        total_change += num_change
        if num_change < 1:
            log.debug('  None. Stopping at iteration {}'.format(n+1))
            break

    # check for exons that still have more than one valid frame
    log.debug('\nExons with unknown phase:')
    for exon, frames in FRAME.items():
        if len(frames) > 1:
            unresolved.append(exon)
            log.debug('  {}, frames: {}'.format(exon, frames))

    log.info('  Resolved {:,} case(s) where an internal exon had >1 possible phase'.format(total_change))
    log.info('  The phase of {:,} ({}) internal exons could not be determined'.format(len(unresolved), percent(len(unresolved), len(exon_list))))

    # DEBUG: Output internal exons with ambiguous frame
    output_as_gff3(unresolved, PREFIX+'.internal.multi_frame.gff3')

    return unresolved


def terminal_exons_in_phase(exon_intrnl, exon_l_term_dict, exon_r_term_dict, max_iter=3):
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
    global SITE_L
    global SITE_R
    # for keeping track of progress
    count = 0
    total = len(ADJ)
    # quick function to calculate the length of an exon
    length = lambda x: x[2] - x[1] + 1
    # for DEBUG-ing
    single_intron = []

    log.debug('')
    log.info('\nChecking which terminal exons are in phase')

    for n in range(max_iter):
        log.debug('\nStarting iteration {}:'.format(n+1))
        new_adj = {}
        for intron, adj in ADJ.items():
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
                log.debug('Looking on both sides of intron at {}'.format(intron))
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
                log.debug('  Left side:')
                for frame, exons in l_compare.items():
                    log.debug('    Frame {}) {}'.format(frame, ' '.join([str(i) for i in exons])))
                log.debug('  Right side:')
                for frame, exons in r_compare.items():
                    log.debug('    Frame {}) {}'.format(frame, ' '.join([str(i) for i in exons])))
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
                    log.debug('  Longest ORF is:')
                    log.debug('    {} is left of intron {}'.format(final[0], intron))
                    log.debug('    {} is right of intron {}'.format(final[1], intron))
                    # DEBUG: print out any single-intron genes
                    single_intron.append(final[0])
                    single_intron.append(final[1])
                else:
                    log.debug('  No valid ORF found! Delete intron {}'.format(intron))
                    del new_adj[intron]
                    intron_discard.add(intron)
            ####################################################################
            elif len(l_adj) < 1:
                log.debug('Looking upstream of intron at {}'.format(intron))
                new_adj[intron] = [ set(), r_adj ]
                # determine what the expected frame(s) should be based on
                # adjacent internal exons
                for exon in r_adj:
                    for f in FRAME[exon]:
                        new_f = frame_shift(f, intron, strand+'<')
                        log.debug('  Downstream exon {}, frame {}'.format(exon, f))
                        expected.add(new_f)
                log.debug('  Expect adjacent frames to be one of: {}'.format(expected))
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
                        log.debug('  Blocks {} with frame {} match'.format(', '.join([str(i) for i in l_compare[f]]), f))
                        log.debug('  Longest is {}'.format(match))
                if match is None:
                    log.debug('  No valid upstream exon found! Delete intron'.format(intron))
                    del new_adj[intron]
                    intron_discard.add(intron)
            ####################################################################
            elif len(r_adj) < 1:
                log.debug('Looking downstream of intron at {}'.format(intron))
                new_adj[intron] = [ l_adj, set() ]
                # determine what the expected frame(s) should be based on
                # adjacent internal exons
                for exon in l_adj:
                    for f in FRAME[exon]:
                        new_f = frame_shift(f, intron, strand+'>')
                        log.debug('  Upstream exon {}, frame {}'.format(exon, f))
                        expected.add(new_f)
                log.debug('  Expect adjacent frames to be one of: {}'.format(expected))
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
                        log.debug('  Blocks {} with frame {} match'.format(', '.join([str(i) for i in r_compare[f]]), f))
                        log.debug('  Longest is {}'.format(match))
                if match is None:
                    log.debug('  No valid downstream exon found! Delete intron'.format(intron))
                    del new_adj[intron]
                    intron_discard.add(intron)
            ####################################################################
            else:
                log.debug('Intron at {} is already flanked by exons. Continue.'.format(intron))
                new_adj[intron] = [ l_adj, r_adj ]

        if new_adj == ADJ:
            log.debug('No change. Stopping iterations.')
            break

        log.debug('Updating adjacency dictionary...')
        ADJ = new_adj
        SITE_L = defaultdict(list)
        SITE_R = defaultdict(list)
        for intron in new_adj:
            chrom, left, right, strand = intron
            SITE_L[(chrom, left, strand)].append(intron)
            SITE_R[(chrom, right, strand)].append(intron)

        # remove exons from the ADJ dictionary
        log.debug('Removing internal exons without adjacent introns...')
        temp = set()
        for exon in exon_intrnl:
            chrom, left, right, strand = exon
            left_introns = SITE_R[(chrom, left-1, strand)]
            right_introns = SITE_L[(chrom, right+1, strand)]
            if len(left_introns) * len(right_introns) > 0:
                log.debug('  Keeping exon {}'.format(exon))
                temp.add(exon)
            else:
                log.debug('  Discarding exon {}'.format(exon))
                exon_discard.add(exon)
                # remove the exon from the adj dict
                for intron in left_introns + right_introns:
                    if intron in ADJ:
                        ADJ[intron][0].discard(exon)
                        ADJ[intron][1].discard(exon)
        exon_intrnl = temp

        log.debug("Removing 5' terminal exons exons without adjacent introns...")
        temp = set()
        for exon in exon_five_term:
            chrom, left, right, strand = exon
            if strand == '+':
                right_introns = SITE_L[(chrom, right+1, strand)]
                if len(right_introns) > 0:
                    log.debug('  Keeping exon {}'.format(exon))
                    temp.add(exon)
                else:
                    log.debug('  Discarding exon {}'.format(exon))
                    exon_discard.add(exon)
                    # remove the exon from the adj dict
                    for intron in right_introns:
                        if intron in ADJ:
                            ADJ[intron][0].discard(exon)
            elif strand == '-':
                left_introns = SITE_R[(chrom, left-1, strand)]
                if len(left_introns) > 0:
                    log.debug('  Keeping exon {}'.format(exon))
                    temp.add(exon)
                else:
                    log.debug('  Discarding exon {}'.format(exon))
                    exon_discard.add(exon)
                    # remove the exon from the adj dict
                    for intron in left_introns:
                        if intron in ADJ:
                            ADJ[intron][1].discard(exon)
        exon_five_term = temp

        log.debug("Removing 3' terminal exons exons without adjacent introns...")
        temp = set()
        for exon in exon_three_term:
            chrom, left, right, strand = exon
            if strand == '+':
                left_introns = SITE_R[(chrom, left-1, strand)]
                if len(left_introns) > 0:
                    log.debug('  Keeping exon {}'.format(exon))
                    temp.add(exon)
                else:
                    log.debug('  Discarding exon {}'.format(exon))
                    exon_discard.add(exon)
                    # remove the exon from the adj dict
                    for intron in left_introns:
                        if intron in ADJ:
                            ADJ[intron][1].discard(exon)
            elif strand == '-':
                right_introns = SITE_L[(chrom, right+1, strand)]
                if len(right_introns) > 0:
                    log.debug('  Keeping exon {}'.format(exon))
                    temp.add(exon)
                else:
                    log.debug('  Discarding exon {}'.format(exon))
                    exon_discard.add(exon)
                    # remove the exon from the adj dict
                    for intron in right_introns:
                        if intron in ADJ:
                            ADJ[intron][0].discard(exon)
        exon_three_term = temp

    log.info("  {:,} 5' terminal exons".format(len(exon_five_term)))
    log.info("  {:,} 3' terminal exons".format(len(exon_three_term)))

    output_as_gff3(single_intron, PREFIX+'.double.gff3') # DEBUG
    output_as_gff3(intron_discard, PREFIX+'.non_coding_introns.gff3', kind='intron') # DEBUG
    output_as_gff3(exon_discard, PREFIX+'.terminal_removed.gff3') # DEBUG

    exon_intrnl = list(exon_intrnl)
    exon_five_term = list(exon_five_term)
    exon_three_term = list(exon_three_term)

    return exon_intrnl, exon_five_term, exon_three_term, intron_discard


def define_genes(introns, exon_intrnl, exon_l_term_dict, exon_r_term_dict):
    """Define coding gene boundaries (i.e. contigous regions of exons joined by
    introns).
    """
    flatten = chain.from_iterable
    temp = defaultdict(list)
    genes = []

    log.debug('')
    log.info('\nDefining gene boundaries')

    # merge sets of introns and internal exons
    for n, f in enumerate((introns, exon_intrnl, exon_l_term_dict, exon_r_term_dict)):
        for i in f:
            chrom, l, r, strand = i
            temp[(chrom, strand)].append((l, r, n))

    # find contiguous regions of exons joined by introns
    for region, ranges in temp.items():
        ranges = sorted(flatten(((l-1, n, 1), (r+1, n, -1)) for l, r, n in ranges))
        c = 0
        x, y = None, None
        for value, kind, label in ranges:
            if c == 0 and kind in (2, 3):
                x = value
            elif c != 0 and kind in (2, 3):
                y = value - 1
            c += label
            if c == 0:
                if all((x, y)):
                    genes.append((region[0], x+1, y, region[1]))
                x, y = None, None

    log.info('  {:,} multi-exon genes'.format(len(genes)))

    return genes


############
# Run loop #
############
if __name__ == '__main__':
    args = parse_commandline_args()
    QUIET = args.quiet
    DEBUG = args.debug
    PREFIX = args.output
    FRAME = defaultdict(set)

    # Parse the index files (positions of start and stop codons) #
    ##############################################################
    index_f, index_r = parse_index(args)

    # Get the positions of all the introns #
    ########################################
    INTRON, splice_f, splice_r, SITE_L, SITE_R = get_introns(args)
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
    exon_intrnl, exon_l_term_dict, exon_r_term_dict = get_putative_exons(blocks_f, blocks_r)
    exon_intrnl, _ = filter_intron_retention.filter(args, INTRON, exon_intrnl)
    check_adjacency(exon_intrnl)

    # Try and resolve ambiguous phase for internal exons #
    ######################################################
    num_unresolved = len(resolve_internal_frames(exon_intrnl))

    # Check which terminal exons are in frame #
    ###########################################
    exon_intrnl, exon_five_term, exon_three_term, intron_discard = terminal_exons_in_phase(exon_intrnl, exon_l_term_dict, exon_r_term_dict)
    for i in intron_discard:
        del INTRON[i]
    # do one more round to use the terminal exons to determine the phase of internal exons
    _ = resolve_internal_frames(exon_intrnl + exon_five_term + exon_three_term)

    # Define gene possible stuctures #
    #################################
    genes = define_genes(INTRON, exon_intrnl, exon_five_term, exon_three_term)

    # Output the results #
    ######################
    exon_intrnl = sort_by_pos(exon_intrnl)
    exon_five_term = sort_by_pos(exon_five_term)
    exon_three_term = sort_by_pos(exon_three_term)
    exon_final = sort_by_pos(exon_intrnl + exon_five_term + exon_three_term)
    output_as_gff3(genes, PREFIX+'.genes.gff3', kind='gene')
    output_as_gff3(exon_intrnl, PREFIX+'.internal.gff3')
    output_as_gff3(exon_five_term, PREFIX+'.five_term.gff3')
    output_as_gff3(exon_three_term, PREFIX+'.three_term.gff3')
    output_as_gff3(exon_final, PREFIX+'.final.gff3')
    log.info('\nSummary:'.format(len(exon_final)))
    log.info('{:,} exons found'.format(len(exon_final)))
    log.info('  {:,} internal exons'.format(len(exon_intrnl)))
    log.info("  {:,} 5' terminal exons".format(len(exon_five_term)))
    log.info("  {:,} 3' terminal exons".format(len(exon_three_term)))

    # Finish up #
    #############
    log.info('\nDone!')
