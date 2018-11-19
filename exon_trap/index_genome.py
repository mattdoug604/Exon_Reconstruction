#!/usr/bin/python3

import sys
from Bio import SeqIO
from collections import defaultdict

if __name__ == '__main__':
    fasta = SeqIO.parse(sys.argv[1], 'fasta')
    prefix = (sys.argv[2] if len(sys.argv) > 1 else 'exons.index')
    stop_dict = defaultdict(list)
    start_dict = defaultdict(list)

    for record in fasta:
        seq = record.seq
        chrom = record.name
        seq_len = len(seq)
        for strand, nuc in [(1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                print('Searching {}, strand {}, frame {}'.format(chrom, strand, frame))
                trans = str(nuc[frame:].translate())
                trans_len = len(trans)
                s_dict = {chrom:[]}
                m_dict = {chrom:[]}
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
                    if end < 0:
                        end = 0
                    s_dict[chrom].append(end)
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
                    if start < 0:
                        start = 0
                    m_dict[chrom].append(start)
                # add both dicts
                stop_dict[(strand, frame)].append(s_dict)
                start_dict[(strand, frame)].append(m_dict)

    for n, i in enumerate(((1,0),(1,1),(1,2),(-1,0),(-1,1),(-1,2)), 1):
        with open('{}.{}'.format(prefix, n), 'w') as f:
            stops = stop_dict[i]
            starts = start_dict[i]
            for d in stops:
                for chrom, pos in d.items():
                    line = chrom, 3, ','.join(map(str, pos))
                    print(*line, sep='\t', file=f)
            for d in starts:
                for chrom, pos in d.items():
                    line = chrom, 0, ','.join(map(str, pos))
                    print(*line, sep='\t', file=f)
