#!/usr/local/bin/python3
# Last updated: 26/4/2018

import sys
from collections import defaultdict

def parse_gff3(path):
    feature = defaultdict(set)

    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip()
                i = line.split('\t')
                seqid, left, right, strand = i[0], i[3], i[4], i[6]
                feature[(seqid, strand)].add((left, right))

    return feature


def print_gff3(feature, path, kind='intron'):
    n = 1

    with open(path, 'w') as out:
        print('##gff-version 3', file=out)
        for region, feats in feature.items():
            seqid, strand = region
            for f in feats:
                left, right = f
                line = seqid, '.', kind, left, right, '.', strand, '.', 'ID={}'.format(n)
                n += 1
                print(*line, sep='\t', file=out)


def compare_features(qry_orig, ref, oargs):
    qry = qry_orig.copy()
    exact = defaultdict(set)
    five_match = defaultdict(set)
    three_match = defaultdict(set)
    overlap = defaultdict(set)
    no_overlap = defaultdict(set)
    c_exact = 0
    c_five = 0
    c_three = 0
    c_overlap = 0
    c_no_overlap = 0

    print('Query {:,}'.format(sum(len(i) for i in qry.values())), file=sys.stderr)
    print('Reference {:,}'.format(sum(len(i) for i in ref.values())), file=sys.stderr)

    # check for exact matches
    print('Checking for exact matches...', end=' ', file=sys.stderr)
    for region, feats in qry.items():
        exact[region] = feats.intersection(ref[region])
        qry[region] = feats - ref[region]
        c_exact += len(exact[region])

    print('{:,}'.format(c_exact), file=sys.stderr)
    path = '.'.join(( *oargs, 'exact', 'gff3' ))
    print_gff3(exact, path)

    # check for matches at one end
    print('Checking for partial matches...', end=' ', file=sys.stderr)
    for region, feats in qry.items():
        strand = region[1]
        l_qry = {i[0]:i for i in feats}
        r_qry = {i[1]:i for i in feats}
        l_ref = {i[0]:i for i in ref[region]}
        r_ref = {i[1]:i for i in ref[region]}
        l_match = set(l_qry).intersection(set(l_ref))
        r_match = set(r_qry).intersection(set(r_ref))
        l_match_feats = set([l_qry[i] for i in l_match])
        r_match_feats = set([r_qry[i] for i in r_match])
        if strand == '+':
            five_match[region] = l_match_feats
            three_match[region] = r_match_feats
        elif strand == '-':
            five_match[region] = r_match_feats
            three_match[region] = l_match_feats
        qry[region] = feats - l_match_feats - r_match_feats
        c_five += len(five_match[region])
        c_three += len(three_match[region])

    print("{:,} 5' matches, {:,} 3' matches".format(c_five, c_three), file=sys.stderr)
    path = '.'.join(( *oargs, '5only', 'gff3' ))
    print_gff3(five_match, path)
    path = '.'.join(( *oargs, '3only', 'gff3' ))
    print_gff3(three_match, path)

    # check for feature that overlap, but differ at both ends
    print('Checking for overlaps...', end=' ', file=sys.stderr)
    for region, feats in qry.items():
        for q in feats:
            for r in ref[region]:
                if q[0] <= r[0] <= q[1] or r[0] <= q[0] <= r[1]:
                    overlap[region].add(q)
        qry[region] = feats - overlap[region]
        c_overlap += len(overlap[region])

    print('{:,}'.format(c_overlap), file=sys.stderr)
    path = '.'.join(( *oargs, 'overlap', 'gff3' ))
    print_gff3(overlap, path)

    # feature that do not overlap
    c_no_overlap = sum(len(i) for i in qry.values())

    print('Checking for non-overlapping... {:,}'.format(c_no_overlap), file=sys.stderr)
    path = '.'.join(( *oargs, 'no_overlap', 'gff3' ))
    print_gff3(qry, path)

    return exact, five_match, three_match, overlap, qry


if __name__ == '__main__':
    qry_prefix = sys.argv[1]
    ref_prefix = '/home2/mattdoug/Thesis/reference/categorized_exons_by_type/coding'
    ref_nc_prefix = '/home2/mattdoug/Thesis/reference/categorized_exons_by_type/non_coding'
    out_prefix = 'validate'

    categories = ('total',
                  'exact',
                  'diff_five_end',
                  'diff_three_end',
                  'overlap',
                  'no_overlap',
                  'total ref',
                  'missed')

    # paths to query files
    path = '.'.join((qry_prefix, 'internal', 'gff3'))
    internal = parse_gff3(path)
    path = '.'.join((qry_prefix, 'five_term', 'gff3'))
    five_term = parse_gff3(path)
    path = '.'.join((qry_prefix, 'three_term', 'gff3'))
    three_term = parse_gff3(path)

    # paths to reference files
    path = '.'.join((ref_prefix, 'internal', 'gff3'))
    ref_internal = parse_gff3(path)
    path = '.'.join((ref_prefix, 'five_term', 'gff3'))
    ref_five_term = parse_gff3(path)
    path = '.'.join((ref_prefix, 'three_term', 'gff3'))
    ref_three_term = parse_gff3(path)

    header = "Internal", "5'Terminal", "3'Terminal"
    for n, i in enumerate((internal, five_term, three_term)):
        print(header[n])
        print(header[n], file=sys.stderr)
        print('Compare.To', 'Exact', 'Diff.5end', 'Diff.3end', 'Overlap', 'No.Overlap', sep='\t')
        # compare i to internal
        compare_to = "Internal"
        results = compare_features(i, ref_internal, [out_prefix, header[n], compare_to])
        line = [compare_to] + [sum([len(v) for v in d.values()]) for d in results]
        print(*line, sep='\t')
        # compare i to 5' terminal
        compare_to = "5'Terminal"
        results = compare_features(i, ref_five_term, [out_prefix, header[n], compare_to])
        line = [compare_to] + [sum([len(v) for v in d.values()]) for d in results]
        print(*line, sep='\t')
        # compare i to 3' terminal
        compare_to = "3'Terminal"
        results = compare_features(i, ref_three_term, [out_prefix, header[n], compare_to])
        line = [compare_to] + [sum([len(v) for v in d.values()]) for d in results]
        print(*line, sep='\t')
