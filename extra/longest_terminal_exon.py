#!/usr/local/bin/python

from __future__ import division, print_function
import pysam, sys

def compare_l(file_A, file_B):
    l_sites = {}

    for f, x in ((file_A, 0), (file_B, 1)):
        for line in f:
            if line[0] != '#':
                i = line.split('\t')
                chrom, start, end, strand = i[0], int(i[3]), int(i[4]), i[6]
                try:
                    l_sites[(chrom, end, strand)][x].append(start)
                except KeyError:
                    l_sites[(chrom, end, strand)] = [ [], [] ]
                    l_sites[(chrom, end, strand)][x].append(start)

    for site, starts in l_sites.items():
        if len(starts[0]) * len(starts[1]) > 0:
            longest_A = max(starts[0]) # reference
            longest_B = max(starts[1]) # reconstructed
            exon = site[0], longest_B, site[1], site[2]
            if longest_A > longest_B:
                yield exon


def compare_r(file_A, file_B):
    r_sites = {}

    for f, x in ((file_A, 0), (file_B, 1)):
        for line in f:
            if line[0] != '#':
                i = line.split('\t')
                chrom, start, end, strand = i[0], int(i[3]), int(i[4]), i[6]
                try:
                    r_sites[(chrom, start, strand)][x].append(end)
                except KeyError:
                    r_sites[(chrom, start, strand)] = [ [], [] ]
                    r_sites[(chrom, start, strand)][x].append(end)

    for site, ends in r_sites.items():
        if len(ends[0]) * len(ends[1]) > 0:
            longest_A = max(ends[0]) # reference
            longest_B = max(ends[1]) # reconstructed
            exon = site[0], site[1], longest_B, site[2]
            if longest_A < longest_B:
                yield exon


def compute_cov(exon, samfile):
    bp_cov = 0
    bp_tot = exon[2] - exon[1] + 1

    for col in samfile.pileup(exon[0], exon[1], exon[2]):
        if exon[1] <= col.pos <= exon[2]:
            bp_cov += 1

    return bp_cov / bp_tot


if __name__ == '__main__':
    l_file_A = open(sys.argv[1], 'r') # reference, l_term
    l_file_B = open(sys.argv[2], 'r') # reconstructed, l_term
    r_file_A = open(sys.argv[3], 'r') # reference, r_term
    r_file_B = open(sys.argv[4], 'r') # reconstructed, r_term
    samfile = pysam.AlignmentFile(sys.argv[5], 'rb')

    for i in compare_l(l_file_A, l_file_B):
        c = compute_cov(i, samfile)
        print('{}:{:,}-{:,} ({} strand)\t{}'.format(i[0], i[1], i[2], i[3], c))

    for i in compare_r(r_file_A, r_file_B):
        c = compute_cov(i, samfile)
        print('{}:{:,}-{:,} ({} strand)\t{}'.format(i[0], i[1], i[2], i[3], c))

    l_file_A.close()
    l_file_B.close()
    r_file_A.close()
    r_file_B.close()
