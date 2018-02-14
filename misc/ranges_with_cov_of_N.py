#!/usr/local/bin/python3

import sys
from collections import defaultdict

if __name__ == '__main__':
    cov_dict = defaultdict(list)

    infile = sys.argv[1]
    thresh = float(sys.argv[2])

    temp = [None, None]
    prev = 0
    with open(infile, 'r') as f:
        for line in f:
            i = line.strip().split('\t')
            chrom, start, end, cov = i[0], int(i[1]), int(i[2]), float(i[3])
            if (start != prev) or (cov < thresh):
                if temp[0] is not None:
                    temp[1] = prev
                    cov_dict[chrom].append(temp)
                    temp = [None, None]
            elif cov >= thresh:
                if temp[0] is None:
                    temp[0] = start
            prev = end
        if temp[0] is not None:
            temp[1] = end
            cov_dict[chrom].append(temp)

    for chrom, ranges in cov_dict.items():
        for r in ranges:
            print('{}:{}-{}'.format(chrom, r[0], r[1]))
