import sys
from collections import defaultdict
import networkx as nx

def pprint(*args, **kwargs):
    level = kwargs.pop('level', {'level':None})
    if level == 'debug':
        if DEBUG:
            print(*args, **kwargs)
    elif level == 'progress':
        if True not in (QUIET, NOPROG, DEBUG):
            print(*args, **kwargs)
    else:
        if not QUIET:
            print(*args, **kwargs)


def percent(val, e_total):
    if e_total == 0:
        return '0%'
    per = val * 100 / e_total
    per = '{0:.1f}%'.format(per)

    return per


def format_as_gff3(features, kind, name='', parent=''):
    """Take a list of features as a tuple (chromosome, start, end, strand),
    and print them as GFF3 formatted entries.
    """
    for f in features:
        chrom, left, right, strand = f
        if kind == 'gene':
            line = chrom, '.', kind, left, right, '.', strand, '.', 'ID={}'.format(name)
        elif kind == 'mRNA':
            # line = chrom, '.', kind, left, right, '.', strand, '.', 'ID={};Parent={}'.format(name, parent)
            line = chrom, '.', kind, left, right, '.', strand, '.', 'ID={}'.format(name)
        elif kind == 'CDS':
            line = chrom, '.', kind, left, right, '.', strand, '.', 'Parent={}'.format(parent)
        yield '\t'.join(map(str, line))


def frame_shift(frame, intron, direction):
    """Calculate the frame a given exon uses based on the frame of the adjacent
    exon and the lenght of the intervening intron.
    """
    i_len = intron[2]-intron[1]+1
    i_mod = i_len % 3

    if direction in ('+<', '->'):
        new_frame = (frame - i_mod) % 3
    elif direction in ('+>', '-<'):
        new_frame = (frame + i_mod) % 3

    pprint('    Calculating frame shift:', level='debug')
    pprint('      Direction is {}'.format(direction), level='debug')
    pprint('      Exon is in frame {}'.format(frame), level='debug')
    pprint('      Intron {} is {}bp long, modulo is {}'.format(intron, i_len, i_mod), level='debug')
    pprint('      New frame is {}'.format(new_frame), level='debug')

    return new_frame


def build_graph(adj, frame):
    graph = nx.DiGraph()

    pprint('Building graph...', level='debug')
    for intron, sides in adj.items():
        l, r = sides
        pprint('  Looking adjacent to intron', intron, level='debug')
        pprint('    Upstream exons:', l, level='debug')
        pprint('    Downstream exons:', r, level='debug')
        for i in l:
            strand = i[3]
            for j in r:
                j_frame = frame[j]
                x_frame = set([frame_shift(x, intron, strand+'>') for x in frame[i]])
                if len(j_frame.intersection(x_frame)) > 0:
                    graph.add_edge(i, j)
                    pprint('    Add edge: {} > {} ({})'.format(i, j, intron), level='debug')

    pprint('  Graph:', level='debug')
    for node, adj in graph.edges():
        pprint('    {} > {}'.format(node, adj), level='debug')

    return graph


def find_all_paths(graph, s, d):
    path  = []
    paths = []
    queue = [(s, d, path)]

    while queue:
        s, d, path = queue.pop()
        path = path + [s]
        if s == d:
            paths.append(path)
        for node in set(graph[s]).difference(path):
            queue.append((node, d, path))

    return paths


def run(args, exon_intrnl, exon_five_term, exon_three_term, genes, adj, frame):
    global QUIET
    global NOPROG
    global DEBUG
    QUIET = args.quiet
    NOPROG = args.noprog
    DEBUG = args.debug
    prefix = args.prefix
    out_path = prefix+'.transcripts.gff3'
    exclude_path = prefix+'.excluded_exons.gff3'
    n = 0
    gc = 0
    t_total = 0
    d_total = 0
    gtotal = len(genes)
    e_total = sum([len(i) for i in (exon_intrnl, exon_five_term, exon_three_term)])
    skip_c = 0

    outfile = open(out_path, 'w')
    excludefile = open(exclude_path, 'w')
    print('##gff-version 3', file=outfile)
    print('##gff-version 3', file=excludefile)

    graph = build_graph(adj, frame)

    pprint('\nAssembling transcripts...', level='debug')
    for gene, features in genes.items():
        pprint('\rAssembling transcripts...', percent(gc, gtotal), end='', level='progress')
        gc += 1
        tc = 1
        c = 0
        all_exons = set()
        for i in (1, 2, 3):
            all_exons = all_exons.union(set([j for j in features[i]]))
        included = set()
        root = features[2]
        exit = features[3]
        pprint("Checking gene {}:".format(gene), level='debug')
        pprint("  {:,} 5' terminal exon(s)".format(len(features[2])), level='debug')
        pprint("  {:,} internal exon(s)".format(len(features[1])), level='debug')
        pprint("  {:,} 3' terminal exon(s)".format(len(features[3])), level='debug')
        pprint("  {:,} introns(s)".format(len(features[0])), level='debug')
        pairs = [(x, y) for y in exit for x in root]
        g_id = 'TEST{}'.format(gc)
        for line in format_as_gff3([gene], 'gene', name=g_id): print(line, file=outfile)
        # skip if there are lots of exons, otherwise script takes forever
        if len(features[1]) > 300:
            skip_c += 1
            continue
        # check for cycles
        for x in root:
            try:
                cycle = nx.find_cycle(graph, x)
                print('[WARNING] cycle found from node {}: {}'.format(x, cycle))
                sys.exit(1)
            except nx.exception.NetworkXNoCycle:
                pass
        for x, y in pairs:
            pprint('  {} to {}'.format(x, y), level='debug')
            path_list = find_all_paths(graph, x, y)
            included = included.union(set([j for i in path_list for j in i]))
            for path in path_list:
                transcript = [(path[0][0], path[0][1], path[-1][2], path[0][3])]
                t_id = 'TEST{}.{}'.format(gc, tc)
                tc += 1
                for line in format_as_gff3(transcript, 'mRNA', name=t_id, parent=g_id):
                    print(line, file=outfile)
                    t_total += 1
                for line in format_as_gff3(path, 'CDS', parent=t_id):
                    print(line, file=outfile)
                pprint(' ', ' > '.join(['{0}'.format(i) for i in path]), level='debug')
        excluded = all_exons - included
        d_total += len(excluded)
        pprint('{:,} exons not included in gene {}'.format(d_total, gene), level='debug')
        for line in format_as_gff3(excluded, 'CDS'):
            print(line, file=excludefile)

    pprint('', level='debug')
    pprint('\rAssembling transcripts... Done!')
    pprint('  {:,} transcripts assembled!'.format(t_total))
    pprint('  {:,}/{:,} ({}) exons could not be placed into a coding transcript'.format(d_total, e_total, percent(d_total, e_total)))
    pprint('  {:,} genes skipped (exceeded maximum depth)'.format(skip_c))
    outfile.close()
    excludefile.close()
