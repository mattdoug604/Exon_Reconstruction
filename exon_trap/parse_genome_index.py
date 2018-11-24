#!/usr/bin/python3

import logging, sys

log_format = '%(message)s'
logging.basicConfig(format=log_format)
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

def parse_index(args):
    global QUIET
    global DEBUG
    QUIET = args.quiet
    DEBUG = args.debug
    PREFIX = args.output
    index_path = args.index
    region = args.region
    index_f = {}
    index_r = {}
    start_c = 0
    end_c = 0

    log.info('Parsing genome index...')
    for i in range(6):
        index = index_path + '.' + str(i+1)
        with open(index, 'r') as f:
            for line in f:
                if line[0] == '#' or len(line) < 3: #if no start/stop codons were found, the 3rd column may be empty
                    continue
                line = line.strip().split('\t')
                seqid, kind, pos_list = line
                seqid = str(seqid)
                kind  = int(kind)
                # if the feature is outside the specified region, skip it
                if region is not None and seqid != region[0]:
                    continue
                pos_list = list(map(int, pos_list.split(',')))
                if region is not None:
                    temp = []
                    for n in pos_list:
                        if region[1] <= n <= region[2]:
                            temp.append(n)
                    pos_list = temp
                # add the features to the dictionaries
                d = [index_f if i < 3 else index_r][0] # 1 == "+", -1 == "-"
                frame = [i-3 if i >= 3 else i][0]
                try:
                    d[seqid][frame] += [(x, kind) for x in pos_list]
                except KeyError:
                    d[seqid] = [ [], [], [] ]
                    d[seqid][frame] += [(x, kind) for x in pos_list]
                # update the counters
                if kind == 0:
                    start_c += len(pos_list)
                elif kind == 3:
                    end_c += len(pos_list)

    log.info('  {:,} putative start codons'.format(start_c))
    log.info('  {:,} stop codons'.format(end_c))

    return index_f, index_r
