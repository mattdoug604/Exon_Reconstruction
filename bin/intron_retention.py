#!/home2/mattdoug/bin/python3
# Last updated: 17/11/2018
# Author: Matt Douglas

import logging
import pysam
from collections import defaultdict

log_format = '%(message)s'
logging.basicConfig(format=log_format, level=logging.DEBUG)
log = logging.getLogger(__name__)


def perc(count, total):
    if total == 0:
        return '0%'
    percent = count * 100 / total
    percent = '{0:.1f}%'.format(percent)

    return percent


def id_events(introns, exons):
    """Return and events where an exon completely overlaps an intron."""
    temp_introns = defaultdict(set)
    temp_exons = defaultdict(set)
    ret_intron_dict = defaultdict(set)

    # convert the sets to dicts so you dont have to iterate through the whole list
    for i in introns:
        chrom, start, end, strand = i
        temp_introns[(chrom, strand)].add((start, end))
    for e in exons:
        chrom, start, end, strand = e
        temp_exons[(chrom, strand)].add((start, end))

    p_total = len(temp_introns)

    log.info('\nIdentifying intron retention events...')
    for n, region in enumerate(temp_introns):
        introns = sorted(temp_introns[region])
        exons = sorted(temp_exons[region]) # sort by position so we can skip to the next introns if we search past the last one
        last = 0
        for i in introns:
            intron = region[0], i[0], i[1], region[1]
            for m, e in enumerate(exons[last:]):
                if e[0] < i[0] and i[1] < e[1]:
                    ret_intron_dict[intron].add(e)
                    log.debug('  exon {} overlaps intron {}'.format(e, intron))
                elif i[1] < e[0] - 1:
                    break
                elif e[1] < i[0]:
                    last = m
    log.info('  {:,} intron retention events identfied'.format(len(ret_intron_dict)))

    return ret_intron_dict


def filter_by_coverage(bamfile, ret_intron_dict, introns, min_depth, min_cov, min_ratio):
    """Return a list of exons that meet the minimum coverage requirements."""
    discard = set()
    f_total = len(ret_intron_dict)
    min_no_cov = 1 - min_cov

    log.info('\nFiltering events (min cov={}, min depth={:,}, min ratio={})...'.format(min_cov, min_depth, min_ratio))
    for n, f in enumerate(ret_intron_dict, 1):
        bp_no_cov = 0
        reads = set()
        fail = False
        score = introns[f]
        region = f[0], f[1]-1, f[2]
        length = f[2] - f[1] + 1
        if length <= 0:
            continue
        for pileupcolumn in bamfile.pileup(*region, truncate=True):
            depth = 0
            has_cov = False
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    reads.add(pileupread.alignment.query_name)
                    depth += 1
                    if depth >= min_depth:
                        has_cov = True
                        break
            if not has_cov:
                bp_no_cov += 1
            if (bp_no_cov / length) > min_no_cov: # goto next feature if exceeds minimum % of no coverage
                fail = True
                break
        ratio = (len(reads) / score if score > 0 else 1)
        if fail or ratio < min_ratio:
            for e in ret_intron_dict[f]:
                exon = f[0], e[0], e[1], f[3]
                discard.add(exon)
                log.debug('  Discrading exon {}'.format(exon))

    log.info('  Removed {:,} events'.format(len(discard)))

    return discard


def filter(args, introns, exons):
    global QUIET
    global DEBUG
    QUIET = args.quiet
    DEBUG = args.debug
    path = args.align
    min_depth = 2
    min_cov = 0.8
    min_ratio = 0

    if path.lower().endswith('.bam'):
        bamfile = pysam.AlignmentFile(path, 'rb')
    else:
        log.critical('Failure! Alignments must be a sorted and indexed BAM file!')
        sys.exit(1)

    ret_intron_dict = id_events(introns, exons)
    discard = filter_by_coverage(bamfile, ret_intron_dict, introns, min_depth, min_cov, min_ratio)
    exons = list(set(exons) - discard)

    return exons, discard
