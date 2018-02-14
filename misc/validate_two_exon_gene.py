#!/usr/local/bin/python3

import sys

def parse_ref_gff3(path):
    unsorted = set()
    exon_set = set()
    left_set = set()
    right_set = set()
    last = None

    with open(path, 'r') as f:
        for line in f:
            if line[0] != '#':
                i = line.split('\t')
                exon = i[0], int(i[3]), int(i[4]), i[6]
                unsorted.add(exon)

    sorted_exons = sorted(unsorted, key=lambda x: (x[0], x[3], int(x[1]), int(x[2])))

    for n, exon in enumerate(sorted_exons):
        left = exon[0], exon[2], exon[3]
        right = exon[0], exon[1], exon[3]
        if n % 2 == 0:
            exon_set.add(exon)
            left_set.add(left)
            last = exon
        else:
            try:
                assert exon[0] == last[0] and exon[1] > last[2]
            except AssertionError:
                print('AssertionError in reference file!')
                print(exon)
                print(last)
                sys.exit(1)
            exon_set.add(exon)
            right_set.add(right)

    return exon_set, left_set, right_set


def parse_qry_gff3(path):
    # unsorted = set()
    unsorted = []
    gene_set = set()
    temp = [None, None]
    last = None

    with open(path, 'r') as f:
        for line in f:
            if line[0] != '#':
                i = line.split('\t')
                exon = i[0], int(i[3]), int(i[4]), i[6]
                # unsorted.add(exon)
                unsorted.append(exon)

    #sorted_exons = sorted(unsorted, key=lambda x: (x[0], x[3], int(x[1]), int(x[2])))
    sorted_exons = unsorted

    for n, exon in enumerate(sorted_exons):
        left = exon[0], exon[2], exon[3]
        right = exon[0], exon[1], exon[3]
        if n % 2 == 0:
            temp[0] = exon
            last = exon
        else:
            try:
                assert exon[0] == last[0] and exon[1] > last[2]
            except AssertionError:
                print('AssertionError in query file!')
                print(exon)
                print(last)
                sys.exit(1)
            temp[1] = exon
            temp = tuple(temp)
            gene_set.add(temp)
            temp = [None, None]

    return gene_set


def perc(count, total):
    if total == 0:
        return '0%'
    percent = count * 100 / total
    percent = '{0:.1f}%'.format(percent)

    return percent


def print_as_gff3(gene_set, out_path):
    with open(out_path, 'w') as f:
        print('##gff-version 3', file=f)
        for exon_tuple in gene_set:
            for exon in exon_tuple:
                chrom, left, right, strand = exon
                line = chrom, '.', 'exon', left, right, 0, strand, '.', '.'
                print(*line, sep='\t', file=f)


if __name__ == '__main__':
    ref = sys.argv[1]
    qry = sys.argv[2]
    both = set()
    three_diff = set()
    five_diff = set()
    overlap = set()
    neither = set()

    ref_set, left_set, right_set = parse_ref_gff3(ref)
    qry_set = parse_qry_gff3(qry)

    for gene in qry_set:
        left, right = gene
        if left in ref_set and right in ref_set:
            both.add(gene)
        elif left in ref_set:
            if left[3] == '+':
                three_diff.add(gene)
            else:
                five_diff.add(gene)
        elif right in ref_set:
            if left[3] == '+':
                five_diff.add(gene)
            else:
                three_diff.add(gene)
        else:
            if (left[0], left[2], left[3]) in left_set:
                overlap.add(gene)
            elif (right[0], right[1], right[3]) in right_set:
                overlap.add(gene)
            else:
                neither.add(gene)

    print('Reference genes\t{}\t '.format( len(ref_set)//2 ))
    print('Query genes\t{}\t '.format( len(qry_set) ))
    print('Both ends match\t{}\t{}'.format( len(both), perc(len(both), len(qry_set)) ))
    print("5' end differs\t{}\t{}".format( len(five_diff), perc(len(five_diff), len(qry_set)) ))
    print("3' end differs\t{}\t{}".format( len(three_diff), perc(len(three_diff), len(qry_set)) ))
    print('Neither end match\t{}\t{}'.format( len(overlap), perc(len(overlap), len(qry_set)) ))
    print('Totally novel gene\t{}\t{}'.format( len(neither), perc(len(neither), len(qry_set)) ))

    print_as_gff3(both, 'two_exons.both_match.gff3')
    print_as_gff3(three_diff, 'two_exons.three_diff.gff3')
    print_as_gff3(five_diff, 'two_exons.five_diff.gff3')
    print_as_gff3(overlap, 'two_exons.both_diff.gff3')
    print_as_gff3(neither, 'two_exons.no_overlap.gff3')
