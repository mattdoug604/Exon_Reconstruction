#!/usr/local/bin/python
#Author: Matt Douglas
#Last updated: 15/2/2018

from __future__ import print_function
import sys
import pysam
from Bio import SeqIO
from collections import defaultdict

def find_codons(fasta):
    stop_dict = defaultdict(list)
    start_dict = defaultdict(list)

    for record in fasta:
        seq = record.seq
        chrom = record.name
        seq_len = len(seq)
        for strand, nuc in [(1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                trans = str(nuc[frame:].translate())
                trans_len = len(trans)
                # get all the stop codons
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find('*', aa_start)
                    if aa_end == -1:
                        break
                    if strand == 1:
                        end = frame + aa_end*3 + 3
                    elif strand == -1:
                        end = seq_len - frame - aa_start*3
                    aa_start = aa_end+1
                    if end < 1:
                        end = 1
                    stop_dict[(chrom, strand)].append(end)
                # get all the start codons
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find('M', aa_start)
                    if aa_end == -1:
                        break
                    if strand == 1:
                        start = frame + (aa_start - 1)*3 + 1
                    elif strand == -1:
                        start = seq_len - frame - aa_end*3
                    aa_start = aa_end+1
                    if start < 1:
                        start = 1
                    start_dict[(chrom, strand)].append(start)

    return stop_dict, start_dict


def check_coverage(pos_dict, bam):
    cov_dict = {}

    for cs, pos_list in pos_dict.items():
        for pos in pos_list:
            region = cs[0], pos-1, pos
            for pileupcolumn in bam.pileup(*region, truncate=True):
                depth = len(pileupcolumn.pileups)
                cov_dict[region] = depth

    return cov_dict


if __name__ == '__main__':
    fasta = SeqIO.parse(sys.argv[1], 'fasta')
    stop_dict, start_dict = find_codons(fasta)

    bam = pysam.AlignmentFile(sys.argv[2], 'rb')
    cov_dict = check_coverage(start_dict, bam)

    for region, depth in cov_dict.items():
        print('{}:{}-{}\t{}'.format(region[0], region[1], region[2], depth))
