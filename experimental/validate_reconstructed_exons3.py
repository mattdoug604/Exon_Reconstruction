#!/usr/local/bin/python3
# Last updated: 16/5/2018

import sys
from collections import defaultdict

def parse_gff3(path, kind=None):
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.split('\t')
                if kind is None or line[2] == kind:
                    yield line[0], int(line[3]), int(line[4]), line[6]


def overlaps_gene(exon, gene_dict):
    seqid, left, right, strand = exon
    key = seqid, strand

    for gene in gene_dict[key]:
        gleft, gright = gene
        if left <= gright and gleft <= right:
            return True
    return False


def extends_gene(exon, gene):
    return False


def cmp_e_intrn(results, exon):
    seqid, left, right, strand = exon
    key = seqid, strand
    matchlvl = 0

    if left in r_left[key]:
        matchlvl += 1
        if right in r_left[key][left]:
            matchlvl += 1
    if right in r_right[key]:
        matchlvl += 1

    if matchlvl == 0:
        if overlaps_gene(exon, gene):
            matchlvl = 0
        elif overlaps_gene(exon, gene_nc):
            matchlvl = 5
        elif extends_gene(exon, gene):
            matchlvl = 6
        else:
            matchlvl = 7

    results[matchlvl].add(exon)

    return results


def cmp_e_term(results, exon):
    seqid, left, right, strand = exon
    key = seqid, strand
    matchlvl = 0

    if left in r_left[key] and right in r_right[key]:
        matchlvl = 3
    elif left in r_left[key]:
        matchlvl = 4
    elif right in r_right[key]:
        matchlvl = 4

    if matchlvl == 0:
        if overlaps_gene(exon, gene):
            matchlvl = 0
        elif overlaps_gene(exon, gene_nc):
            matchlvl = 5
        elif extends_gene(exon, gene):
            matchlvl = 6
        else:
            matchlvl = 7

    results[matchlvl].add(exon)

    return results


def print_gff3(exon_list, path, kind='CDS'):
    n = 1

    with open(path, 'w') as out:
        print('##gff-version 3', file=out)
        for exon in exon_list:
            seqid, left, right, strand = exon
            line = seqid, '.', kind, left, right, '.', strand, '.', 'ID={}'.format(n)
            n += 1
            print(*line, sep='\t', file=out)


if __name__ == '__main__':
    path_intron = '/home2/mattdoug/Thesis/intron_database/v2/intron_database.min5.gff3'
    path_exon = '/home2/mattdoug/Thesis/exon_database/v2/exon_database.min5.gff3'
    path_gene = '/home2/mattdoug/Thesis/reference/c_elegans.PRJNA13758.WS250.protein_coding.gff3'
    path_gene_nc = '/home2/mattdoug/Thesis/reference/c_elegans.PRJNA13758.WS250.non_coding.gff3'

    ########################################################
    # Part 1: Determine which exons are internal, terminal #
    ########################################################
    i_left = defaultdict(set)
    i_right = defaultdict(set)
    r_intron_left = defaultdict(set)
    r_intron_right = defaultdict(set)
    r_left = defaultdict(dict)
    r_right = defaultdict(set)
    gene = defaultdict(set)
    gene_nc = defaultdict(set)
    result_internal = {i:set() for i in range(7)}
    result_terminal = {i:set() for i in range(7)}
    # results[0] == neither end matches (two novel splice sites)
    # results[1] == one end matches (alternative splice site)
    # results[2] == both ends match, but not the same exon (intron retention)
    # results[3] == both ends match the same exon (exact match)
    # results[4] == terminal matches terminal at intron end (counted as exact)
    # results[5] == terminal does not match terminal at intron end (novel terminal)
    # results[5] == overlaps non-coding gene
    # results[6] == extends coding gene
    # results[7] == does not overlap any known gene

    # Parse the RNA-Seq intron database
    print('Looking for "intron" in {}'.format(path_intron))
    for n, x in enumerate(set([i for i in parse_gff3(path_intron, kind='intron')]), 1):
        seqid, left, right, strand = x
        k = seqid, strand
        i_left[k].add(left)
        i_right[k].add(right)
    print('  {:,} introns found!'.format(n))

    # Parse intron features from the WormBase gene models
    print('Looking for "intron" in {}'.format(path_gene))
    for n, x in enumerate(set([i for i in parse_gff3(path_gene, kind='intron')]), 1):
        seqid, left, right, strand = x
        k = seqid, strand
        r_intron_left[k].add(left)
        r_intron_right[k].add(right)
    print('  {:,} introns found!'.format(n))

    # Parse CDS features from the WormBase gene models
    print('Looking for "CDS" in {}'.format(path_gene))
    for n, x in enumerate(set([i for i in parse_gff3(path_gene, kind='CDS')]), 1):
        seqid, left, right, strand = x
        k = seqid, strand
        if left not in r_left[k]:
            r_left[k][left] = set()
        r_left[k][left].add(right)
        r_right[k].add(right)
        if (left-1 in r_intron_right[k]) is not (right+1 in r_intron_left[k]): # introns on one side only
            terminal.add(left)
            terminal.add(right)
    print('  {:,} coding exons found!'.format(n))

    # Parse coding genes from the WormBase gene models
    print('Looking for "gene" in {}'.format(path_gene))
    for n, x in enumerate(set([i for i in parse_gff3(path_gene, kind='gene')]), 1):
        seqid, left, right, strand = x
        k = seqid, strand
        gene[k].add((left, right))
    print('  {:,} coding genes found!'.format(n))

    # Parse non-coding genes from the WormBase gene models
    print('Looking for "gene" (non-coding genes) in {}'.format(path_gene_nc))
    for n, x in enumerate(set([i for i in parse_gff3(path_gene_nc, kind='gene')]), 1):
        seqid, left, right, strand = x
        k = seqid, strand
        gene_nc[k].add((left, right))
    print('  {:,} non-coding genes found!'.format(n))

    # Parse RNA-Seq exon database
    print('Comparing exons to WormBase, from {}'.format(path_exon))
    for n, exon in enumerate(set([i for i in parse_gff3(path_exon)]), 1):
        seqid, left, right, strand = exon
        k = seqid, strand
        # determine which exons are terminal based on how many introns are adjacent to them
        if left-1 in i_right[k] and right+1 in i_left[k]: # introns on both sides
            results = cmp_e_intrn(results, exon)
        elif left-1 in i_right[k]:  # introns on the left side only
            results = cmp_e_term(results, exon)
        elif right+1 in i_left[k]:  # introns on the right side only
            results = cmp_e_term(results, exon)
        else:                                     # no introns on either side
            results = cmp_e_intrn(results, exon)
    print('  {:,}\tExons total'.format(n))

    #################################
    # Part 2: Print out the results #
    #################################
    # results[0] == neither end matches (two novel splice sites)
    # results[1] == one end matches (alternative splice site)
    # results[2] == both ends match, but not the same exon (intron retention)
    # results[3] == both ends match the same exon (exact match)
    # results[4] == terminal exon matches at intron end (counted as exact)
    # results[5] == overlaps non-coding gene
    # results[6] == does not overlap any known gene
    results0 = sorted(results[0], key=lambda x: (x[0], x[1], x[2]))
    results1 = sorted(results[1], key=lambda x: (x[0], x[1], x[2]))
    results2 = sorted(results[2], key=lambda x: (x[0], x[1], x[2]))
    results3 = sorted(results[3], key=lambda x: (x[0], x[1], x[2]))
    results4 = sorted(results[4], key=lambda x: (x[0], x[1], x[2]))
    results5 = sorted(results[5], key=lambda x: (x[0], x[1], x[2]))
    results6 = sorted(results[6], key=lambda x: (x[0], x[1], x[2]))
    results7 = sorted(results[7], key=lambda x: (x[0], x[1], x[2]))

    print('  {:,}\tNeither end matches (both splice sites novel)'.format(len(results0)))
    print('  {:,}\tOne end matches (one novel splice site)'.format(len(results1)))
    print('  {:,}\tBoth ends match, but different exons (intron retention)'.format(len(results2)))
    print('  {:,}\tMatches reference exon'.format(len(results3)))
    print('    {:,}\tTerminal exon differs at the non-intron end (counts as exact)'.format(len(results4)))
    print('  {:,}\tOverlaps non-coding gene'.format(len(results5)))
    print('  {:,}\tDoes not overlap a known gene'.format(len(results6)))

    print_gff3(results0, 'validate.results0.gff3')
    print_gff3(results1, 'validate.results1.gff3')
    print_gff3(results2, 'validate.results2.gff3')
    print_gff3(results3, 'validate.results3.gff3')
    print_gff3(results4, 'validate.results4.gff3')
    print_gff3(results5, 'validate.results5.gff3')
    print_gff3(results6, 'validate.results6.gff3')
    print_gff3(results7, 'validate.results7.gff3')
