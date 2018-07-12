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


def match_intrl(qry, ref, result):
    qry_d = defaultdict(set)
    ref_d = defaultdict(set)
    result = {'qry_only':[], 'both':[], 'ref_only':[]}

    for i in qry:
        seqid, left, right, strand = i
        qry_d[(seqid, strand)].add((left, right))
    for i in ref:
        seqid, left, right, strand = i
        ref_d[(seqid, strand)].add((left, right))

    for k, v in qry_d.items():
        seqid, strand = k
        result['qry_only'] += [(seqid, i[0], i[1], strand) for i in v - ref_d[k]]
        result['both'] += [(seqid, i[0], i[1], strand) for i in v.intersection(ref_d[k])]
        result['ref_only'] += [(seqid, i[0], i[1], strand) for i in ref_d[k] - v]

    return result


def match_lterm(qry, ref, fresult, rresult):
    qry_d = defaultdict(set)
    ref_d = defaultdict(dict)
    ref_match = set()

    for i in qry:
        seqid, left, right, strand = i
        qry_d[(seqid, strand)].add((left, right))
    for i in ref:
        seqid, left, right, strand = i
        if right not in ref_d[(seqid, strand)]:
            ref_d[(seqid, strand)][right] = set()
        ref_d[(seqid, strand)][right].add(left)

    for k, v in qry_d.items():
        seqid, strand = k
        if strand == '+':
            result = fresult
        else:
            result = rresult
        for i in v:
            if i[1] in ref_d[k]:
                if i[0] in ref_d[k][i[1]]:
                    result['both_exact'].append((seqid, i[0], i[1], strand))
                else:
                    result['both'].append((seqid, i[0], i[1], strand))
                for j in ref_d[k][i[1]]:
                    ref_match.add((seqid, j, i[1], strand))
            else:
                result['qry_only'].append((seqid, i[0], i[1], strand))

    for i in [x for x in ref - ref_match if x[3] == '+' ]:
        fresult['ref_only'].append(i)
    for i in [x for x in ref - ref_match if x[3] != '+' ]:
        rresult['ref_only'].append(i)

    return fresult, rresult


def match_rterm(qry, ref, fresult, rresult):
    qry_d = defaultdict(set)
    ref_d = defaultdict(dict)
    ref_match = set()

    for i in qry:
        seqid, left, right, strand = i
        qry_d[(seqid, strand)].add((left, right))
    for i in ref:
        seqid, left, right, strand = i
        if left not in ref_d[(seqid, strand)]:
            ref_d[(seqid, strand)][left] = set()
        ref_d[(seqid, strand)][left].add(right)

    for k, v in qry_d.items():
        seqid, strand = k
        if strand == '+':
            result = rresult
        else:
            result = fresult
        for i in v:
            if i[0] in ref_d[k]:
                if i[1] in ref_d[k][i[0]]:
                    result['both_exact'].append((seqid, i[0], i[1], strand))
                else:
                    result['both'].append((seqid, i[0], i[1], strand))
                for j in ref_d[k][i[0]]:
                    ref_match.add((seqid, i[0], j, strand))
            else:
                result['qry_only'].append((seqid, i[0], i[1], strand))

    for i in [x for x in ref - ref_match if x[3] == '+' ]:
        rresult['ref_only'].append(i)
    for i in [x for x in ref - ref_match if x[3] != '+' ]:
        fresult['ref_only'].append(i)

    return fresult, rresult


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
    gene = defaultdict(set)
    gene_nc = defaultdict(set)

    qry_intron_l = defaultdict(set)
    qry_intron_r = defaultdict(set)
    ref_intron_l = defaultdict(set)
    ref_intron_r = defaultdict(set)

    ref_intrl = set()
    ref_lterm = set()
    ref_rterm = set()
    ref_singl = set()
    qry_intrl = set()
    qry_lterm = set()
    qry_rterm = set()
    qry_singl = set()

    result_intrl = {'qry_only':[], 'both':[], 'ref_only':[]}
    result_fterm = {'qry_only':[], 'both':[], 'both_exact':[], 'ref_only':[]}
    result_rterm = {'qry_only':[], 'both':[], 'both_exact':[], 'ref_only':[]}
    result_singl = {'qry_only':[], 'both':[], 'ref_only':[]}

    # # Parse coding genes from the WormBase gene models
    # print('Looking for "gene" in {}'.format(path_gene))
    # for n, x in enumerate(set([i for i in parse_gff3(path_gene, kind='gene')]), 1):
    #     seqid, left, right, strand = x
    #     gene[(seqid, strand)].add((left, right))
    # print('  {:,} coding genes found!'.format(n))
    #
    # # Parse non-coding genes from the WormBase gene models
    # print('Looking for "gene" (non-coding genes) in {}'.format(path_gene_nc))
    # for n, x in enumerate(set([i for i in parse_gff3(path_gene_nc, kind='gene')]), 1):
    #     seqid, left, right, strand = x
    #     gene_nc[(seqid, strand)].add((left, right))
    # print('  {:,} non-coding genes found!'.format(n))

    # Parse intron features from the WormBase gene models
    print('Looking for "intron" in {}'.format(path_gene))
    for n, x in enumerate(set([i for i in parse_gff3(path_gene, kind='intron')]), 1):
        seqid, left, right, strand = x
        ref_intron_l[(seqid, strand)].add(left)
        ref_intron_r[(seqid, strand)].add(right)
    print('  {:,} introns found!'.format(n))

    # Parse CDS features from the WormBase gene models
    print('Looking for "CDS" in {}'.format(path_gene))
    for n, x in enumerate(set([i for i in parse_gff3(path_gene, kind='CDS')]), 1):
        seqid, left, right, strand = x
        k = seqid, strand
        if left-1 in ref_intron_r[k] and right+1 in ref_intron_l[k]: #intron on both sides
            ref_intrl.add(x)
        elif left-1 in ref_intron_r[k]:  #intron on left side only
            ref_rterm.add(x)
        elif right+1 in ref_intron_l[k]: #intron on right side only
            ref_lterm.add(x)
        else:                          #no intron on either side
            ref_singl.add(x)
    print('  {:,} coding exons found!'.format(n))
    print('    {:,} internal'.format(len(ref_intrl)))
    print('    {:,} left terminal'.format(len(ref_lterm)))
    print('    {:,} right terminal!'.format(len(ref_rterm)))
    print('    {:,} single'.format(len(ref_singl)))

    # Parse the RNA-Seq intron database
    print('Looking for "intron" in {}'.format(path_intron))
    for n, x in enumerate(set([i for i in parse_gff3(path_intron, kind='intron')]), 1):
        seqid, left, right, strand = x
        qry_intron_l[(seqid, strand)].add(left)
        qry_intron_r[(seqid, strand)].add(right)
    print('  {:,} introns found!'.format(n))

    # Parse RNA-Seq exon database
    print('Looking for "CDS" in {}'.format(path_exon))
    for n, x in enumerate(set([i for i in parse_gff3(path_exon)]), 1):
        seqid, left, right, strand = x
        k = seqid, strand
        # determine which exons are terminal based on how many introns are adjacent to them
        if left-1 in qry_intron_r[k] and right+1 in qry_intron_l[k]: #intron on both sides
            qry_intrl.add(x)
        elif left-1 in qry_intron_r[k]:     #intron on left side only
            qry_rterm.add(x)
        elif right+1 in qry_intron_l[k]:    #intron on right side only
            qry_lterm.add(x)
        else:                             #no intron on either side
            qry_singl.add(x)
    print('  {:,} exons found!'.format(n))
    print('    {:,} internal'.format(len(qry_intrl)))
    print('    {:,} left terminal'.format(len(qry_lterm)))
    print('    {:,} right terminal!'.format(len(qry_rterm)))
    print('    {:,} single'.format(len(qry_singl)))

    print_gff3(ref_intrl, 'category.ref.internal.gff3')
    print_gff3(ref_lterm, 'category.ref.lterminal.gff3')
    print_gff3(ref_rterm, 'category.ref.rterminal.gff3')
    print_gff3(ref_singl, 'category.ref.single.gff3')
    print_gff3(qry_intrl, 'category.qry.internal.gff3')
    print_gff3(qry_lterm, 'category.qry.lterminal.gff3')
    print_gff3(qry_rterm, 'category.qry.rterminal.gff3')
    print_gff3(qry_singl, 'category.qry.single.gff3')

    print('Comparing features...')
    result_intrl = match_intrl(qry_intrl, ref_intrl, result_intrl)
    result_fterm, result_rterm = match_lterm(qry_lterm, ref_lterm, result_fterm, result_rterm)
    result_fterm, result_rterm = match_rterm(qry_rterm, ref_rterm, result_fterm, result_rterm)
    result_singl = match_intrl(qry_singl, ref_singl, result_singl)

    print('Type', 'Qry.Only', 'Both', 'Both.Exact', 'Ref.Only', sep='\t')
    print('Internal\t{:,}\t{:,}\t.\t{:,}'.format(len(result_intrl['qry_only']),
                                                 len(result_intrl['both']),
                                                 len(result_intrl['ref_only'])))
    print('5.terminal\t{:,}\t{:,}\t{:,}\t{:,}'.format(len(result_fterm['qry_only']),
                                                      len(result_fterm['both']),
                                                      len(result_fterm['both_exact']),
                                                      len(result_fterm['ref_only'])))
    print('3.terminal\t{:,}\t{:,}\t{:,}\t{:,}'.format(len(result_rterm['qry_only']),
                                                      len(result_rterm['both']),
                                                      len(result_rterm['both_exact']),
                                                      len(result_rterm['ref_only'])))
    print('Single\t{:,}\t{:,}\t.\t{:,}'.format(len(result_singl['qry_only']),
                                               len(result_singl['both']),
                                               len(result_singl['ref_only'])))

    print_gff3(result_intrl['qry_only'], 'result.internal.qry_only.gff3')
    print_gff3(result_intrl['both'], 'result.internal.both.gff3')
    print_gff3(result_intrl['ref_only'], 'result.internal.ref_only.gff3')
    print_gff3(result_fterm['qry_only'], 'result.fterm.qry_only.gff3')
    print_gff3(result_fterm['both'], 'result.fterm.both.gff3')
    print_gff3(result_fterm['both_exact'], 'result.fterm.both_exact.gff3')
    print_gff3(result_fterm['ref_only'], 'result.fterm.ref_only.gff3')
    print_gff3(result_rterm['qry_only'], 'result.rterm.qry_only.gff3')
    print_gff3(result_rterm['both'], 'result.rterm.both.gff3')
    print_gff3(result_rterm['both_exact'], 'result.rterm.both_exact.gff3')
    print_gff3(result_rterm['ref_only'], 'result.rterm.ref_only.gff3')
    print_gff3(result_singl['qry_only'], 'result.single.qry_only.gff3')
    print_gff3(result_singl['both'], 'result.single.both.gff3')
    print_gff3(result_singl['ref_only'], 'result.single.ref_only.gff3')
