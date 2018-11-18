#!/usr/local/bin/python
# Last updated: 31/1/2018

from __future__ import division, print_function
import pysam
import sys

def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def parse_gff3(path):
    """Parse a GFF3 file and return all the feature coordinates."""
    features = set()

    with open(path, 'r') as f:
        for line in f:
            if line[0] != '#':
                i = line.split('\t')
                feature = i[0], int(i[3]), int(i[4]), int(i[5]), i[6] # chrom, start, end, support, strand
                features.add(feature)

    features = list(features)

    return features


def check_coverage(features):
    cov_dict = {}
    f_total = len(features)

    for n, f in enumerate(features, 1):
        # eprint('\r{:,}/{:,} features checked...'.format(n, f_total), end='')
        region = f[0], f[1]-1, f[2]
        temp = set()
        bp_cov = 0
        for pileupcolumn in bamfile.pileup(*region, truncate=True):
            is_cov = False
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    is_cov = True
                    if pileupread.alignment.is_read1:
                        name = pileupread.alignment.query_name + '/1'
                    elif pileupread.alignment.is_read2:
                        name = pileupread.alignment.query_name + '/2'
                    else:
                        name = pileupread.alignment.query_name
                    temp.add(name)
            if is_cov:
                bp_cov += 1
        cov_dict[f] = (len(temp), bp_cov)
    # eprint('\r{:,}/{:,} features checked!    '.format(n, f_total))

    return cov_dict


if __name__ == '__main__':
    gff_path = sys.argv[1]
    bam_path = sys.argv[2]

    bamfile = pysam.AlignmentFile(bam_path, 'rb')
    features = parse_gff3(gff_path)
    #features = [('I', 1839193, 1839250, 11538, '+')]
    #features = [('I', 390352, 390400, 15, '-')]
    #features = [('I', 484859, 484906, 31, '-')]
    # eprint('{:,} features total'.format(len(features)))

    cov_dict = check_coverage(features)
    sorted_by_pos = sorted(cov_dict, key=lambda x: (x[0], x[1], x[2]))

    print('chrom', 'start', 'end', 'strand', 'support_for_retention', 'support_for_intron', 'ratio', 'bp_cov', sep='\t')
    for intron in sorted_by_pos:
        chrom, start, end, support, strand = intron
        reads, bp_cov = cov_dict[intron]
        per_cov = bp_cov / (end - start + 1)
        try:
            ratio = reads / support
        except ZeroDivisionError:
            ratio = 1
        line = chrom, start, end, strand, reads, support, round(ratio, 4), round(per_cov, 4)
        print(*line, sep='\t')
