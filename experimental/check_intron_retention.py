#!/usr/local/bin/python3
# Last updated: 20/4/2018

import argparse
import sys
from collections import defaultdict

def parse_commandline_arguments():
    parser = argparse.ArgumentParser(description='Compare a set of feats agains a set of exons and look for feat retention events.')
    parser.add_argument('-i',
                        type=str,
                        nargs='?',
                        help='a GFF3 file of feats')
    parser.add_argument('-e',
                        type=str,
                        nargs='?',
                        help="a GFF3 file of exons")
    args = parser.parse_args()

    if None in (args.i, args.e):
        parser.print_help()
        sys.exit(1)

    return args.i, args.e


def parse_gff3(path):
    region_dict = defaultdict(list)
    site_dict = defaultdict(list)

    with open(path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            line = line.split('\t')
            chrom = line[0]
            left = int(line[3])
            right = int(line[4])
            try:
                score = int(line[5])
            except ValueError:
                score = 0
            strand = line[6]
            region_dict[(chrom, strand)].append((left, right, score))
            site_dict[(chrom, left, strand)].append((right, score))
            site_dict[(chrom, right, strand)].append((left, score))

    return region_dict, site_dict


def retention_scores(e_dict, i_dict, i_site):
    ret_dict = {}

    total = sum([len(i) for i in i_dict.values()])
    count = 0

    for r, introns in i_dict.items():
        for i in introns:
            count += 1
            print('\r{:,}/{:,}'.format(count, total), end='', file=sys.stderr)
            intron = r[0], i[0], i[1], r[1]
            for e in e_dict[r]:
                if e[0] < i[0] and i[1] < e[1]:
                    exon = r[0], e[0], e[1], r[1]
                    if exon not in ret_dict:
                        ret_dict[exon] = [0, i[2], 0]
                    else:
                        ret_dict[exon][1] += i[2]
    print('', file=sys.stderr)

    for exon, scores in ret_dict.items():
        chrom, left, right, strand = exon
        for i in i_site[(chrom, left-1, strand)]:
            scores[0] += i[1]
        for i in i_site[(chrom, right+1, strand)]:
            scores[2] += i[1]

    return ret_dict


def print_as_gff3(feature_dict, score_dict, path, kind):
    with open(path, 'w') as f:
        print('##gff-version 3', file=f)
        for r, features in feature_dict.items():
            chrom, strand = r
            for i in features:
                start, end = i
                try:
                    score = score_dict[(chrom, start, end, strand)]
                except KeyError:
                    score = '.'
                line = chrom, '.', kind, start, end, score, strand, '.', '.'
                print(*line, sep='\t', file=f)


if __name__ == '__main__':
    intron_path, exon_path = parse_commandline_arguments()
    i_dict, i_site = parse_gff3(intron_path)
    e_dict = parse_gff3(exon_path)[0]
    ret_dict = retention_scores(e_dict, i_dict, i_site)

    for exon, scores in ret_dict.items():
        min_adj = min(scores[0], scores[2])
        score_adj = min_adj - scores[1]
        if min_adj > 0:
            print('{}:{}-{}{}\t{}\t{}'.format(*exon, ','.join(map(str, scores)), score_adj))
