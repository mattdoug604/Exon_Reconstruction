#!/usr/bin/python3

import gzip
import sys
from Bio import SeqIO
from collections import defaultdict
from pathlib import Path

VERSION='0.2'

def open_file(input):
    is_gzipped = False
    in_type = None
    ext_list = [x.lower() for x in Path(input).suffixes]
    ext = ext_list.pop()
    # check if FASTx file is gzipped
    if ext in ('.gz', '.gzip'):
        is_gzipped = True
        ext = ext_list.pop()
    # check if file is FASTA or FASTQ
    if ext in ('.fasta', '.fa'):
        in_type = 'fasta'
    elif ext in ('.fastq', '.fq'):
        in_type = 'fastq'
    else:
        print('Could not determine if input is FASTA or FASTQ!', file=sys.stderr)
        raise(TypeError)

    if is_gzipped:
        with gzip.open(input, 'rt') as f:
            for line in SeqIO.parse(f, in_type):
                yield line
    else:
        for line in SeqIO.parse(input, in_type):
            yield line


if __name__ == '__main__':
    input = sys.argv[1]
    output = (sys.argv[2] if len(sys.argv) > 2 else 'exons.index')
    x = 0
    y = len([x for x in open_file(input)])*6 #6 reading frames per DNA sequence
    fastx = open_file(input)
    stop_dict = defaultdict(list)
    start_dict = defaultdict(list)

    for record in fastx:
        seq = record.seq
        seqid = record.name
        seqlen = len(seq)
        for strand, nuc in [(1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                x += 1
                print('[{:,}/{:,}] Indexing sequence: {}, strand {}, frame {}'.format(x, y, seqid, strand, frame))
                trim_seq = nuc[frame:][:len(nuc[frame:]) // 3 * 3] # trim sequence so its length is divisible by 3
                trans = str(trim_seq.translate())
                trans_len = len(trans)
                s_dict = {seqid:[]}
                m_dict = {seqid:[]}
                # get all the stop codons
                aa_start = trans.find('*') + 1
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find('*', aa_start)
                    if aa_end == -1:
                        break
                    if strand == 1:
                        end = frame + aa_end*3 + 3
                    elif strand == -1:
                        end = seqlen - frame - aa_start*3
                    aa_start = aa_end + 1
                    s_dict[seqid].append(end)
                # get all the start codons
                aa_start = trans.find('M') + 1
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find('M', aa_start)
                    if aa_end == -1:
                        break
                    if strand == 1:
                        start = frame + (aa_start - 1)*3 + 1
                    elif strand == -1:
                        start = seqlen - frame - aa_end*3
                    aa_start = aa_end + 1
                    m_dict[seqid].append(start)
                # add both dicts
                stop_dict[(strand, frame)].append(s_dict)
                start_dict[(strand, frame)].append(m_dict)

    for n, i in enumerate(((1,0),(1,1),(1,2),(-1,0),(-1,1),(-1,2)), 1):
        with open('{}.{}'.format(output, n), 'w') as f:
            stops = stop_dict[i]
            starts = start_dict[i]
            print('#version: {}'.format(VERSION), file=f)
            for d in starts:
                for seqid, pos in d.items():
                    line = seqid, 0, ','.join(map(str, pos))
                    print(*line, sep='\t', file=f)
            for d in stops:
                for seqid, pos in d.items():
                    line = seqid, 3, ','.join(map(str, pos))
                    print(*line, sep='\t', file=f)

    print('Wrote index files with the prefix: {}'.format(output))
