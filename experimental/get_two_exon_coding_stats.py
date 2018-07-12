#!/usr/local/bin/python3

import re, sys
from collections import defaultdict

attr_gene_name  = re.compile('(?<=sequence_name=)[^;]+(?=(;|))')
attr_biotype    = re.compile('(?<=biotype=)[^;]+(?=(;|))')
attr_trans_name = re.compile('(?<=ID=)[^;]+(?=(;|))')
attr_parent     = re.compile('(?<=Parent=)[^;]+(?=(;|))')
attr_status     = re.compile('(?<=prediction_status=)[^;]+(?=(;|))')

class CDS(object):
    __slots__ = ['chrom', 'start', 'end', 'strand', 'phase', 'status']
    def __init__(self, chrom, start, end, strand, phase, status):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.phase = phase
        self.status = status.lower()

    def __len__(self):
        return self.end - self.start + 1

    def coords(self):
        return self.chrom, self.start, self.end, self.strand


def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def fr(region):
    """Print a tuple of (chromosome, start position, end position, strand) into
    a readable format.
    """
    return '{}:{}-{}{}'.format(region[0], region[1], region[2], region[3])


def parent_gene_name(transcript_id):
    gene_name = '.'.join(transcript_id.split('.')[:2])
    if gene_name[-1].isalpha():
        gene_name = gene_name[:-1]  #Y71H2AM.19a.1 <- ".1" indicates different UTRs

    return gene_name


def transcript_id_names(attrs):
    transcript_ids = attr_parent.search(attrs).group(0)
    transcript_ids = transcript_ids.split(',')
    transcript_ids = [i.split(':')[1] for i in transcript_ids]

    return transcript_ids


def parse_wormbase_transcripts(GFF3_file):
    with open(GFF3_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line[0] != '#':
                i = line.split('\t')
                kind = i[2]
                pos = i[0], int(i[3]), int(i[4]), i[6]
                phase = [int(i[7]) if i[7] != '.' else '.'][0]
                try:
                    if kind == 'gene':
                        biotype = attr_biotype.search(i[8]).group(0)
                        seq_name = attr_gene_name.search(i[8]).group(0)
                        gene_type_dict[seq_name] = biotype
                        gene_pos_dict[seq_name] = pos
                    elif 'RNA' in kind:
                        transcript_id = attr_trans_name.search(i[8]).group(0).split(':')[1]
                        gene_name = parent_gene_name(transcript_id)
                        gene_dict[gene_name].append(transcript_id)
                    elif kind == 'intron':
                        for i in transcript_id_names(i[8]):
                            trans_intron_dict[i].append(pos)
                    elif kind == 'CDS':
                        try:
                            status = attr_status.search(i[8]).group(0)
                        except AttributeError:
                            status = 'n/a'
                        cds = CDS(pos[0], pos[1], pos[2], pos[3], phase, status)
                        for i in transcript_id_names(i[8]):
                            trans_cds_dict[i].append(cds)
                except AttributeError:
                    eprint('[AttributeError] Could not parse line:')
                    eprint(line)
                    sys.exit(1)


def feats_per_gene(transcripts, feature_dict):
    feats_per_gene = set()

    for t in transcripts:
        for i in feature_dict[t]:
            feats_per_gene.add(i)

    return feats_per_gene


def uniquify(cds):
    cds_dict = {}

    for i in cds:
        cds_dict[i.coords()] = i
    unique_cds = list(cds_dict.values())
    unique_cds = sorted(unique_cds, key=lambda x: (x.start))

    return unique_cds


def format_as_gff3(intron_list):
    yield '##gff-version 3'
    for n, i in enumerate(intron_list):
        chrom, left, right, strand = i
        line = chrom, '.', 'intron', left, right, '.', strand, '.', 'ID={}{}'.format('intron', n+1)
        yield '\t'.join(map(str, line))


if __name__ == '__main__':
    gene_dict = defaultdict(list)
    gene_type_dict = {}
    gene_pos_dict = {}
    trans_cds_dict = defaultdict(list)
    trans_intron_dict = defaultdict(list)
    intron_list = []
    gff3_path = '/home2/mattdoug/Thesis/reference/c_elegans.PRJNA13758.WS250.protein_coding.gff3'
    #gff3_path = '/home2/mattdoug/scripts/two_exon_genes/test.gff3'

    parse_wormbase_transcripts(gff3_path)

    # skip non-coding transcripts or transcripts with > 2 coding exons
    for gene, transcripts in gene_dict.items():
        biotype = gene_type_dict[gene]
        num_introns = len(feats_per_gene(transcripts, trans_intron_dict))
        if biotype != 'protein_coding':
            continue
        if num_introns !=  1:
            continue
        for t in transcripts:
            cds_list = uniquify(trans_cds_dict[t])
            if len(cds_list) == 2:
                ######################### DO STUFF #########################
                left = cds_list[0]
                right = cds_list[1]
                intron = left.chrom, (left.end + 1), (right.start - 1), left.strand
                intron_list.append(intron)
                line = fr(left.coords()), fr(right.coords()), left.phase, right.phase, left.status, right.status, fr(intron)
                print(*line, sep='\t')

    with open('two_exon_introns.gff3', 'w') as f:
        for i in format_as_gff3(intron_list):
            print(i, file=f)
