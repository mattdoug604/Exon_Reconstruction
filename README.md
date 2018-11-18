# Identify putative protein-coding exons from a reference genome and a set of splice sites.

This program takes a set of aligned reads (BAM format), introns (GFF3 format), and an index of "translation blocks" encoded by the reference genome:

Requires Python 3 and the following modules:
- BioPython
- pysam

## Before running:
The main script requires that your genome (or scaffold) sequence be indexed. Run 'exon_trap/generate_index.py' to index your desired genome sequence (in FASTA format). Note: This only has to be done once per genome.

## To run:
exon_trap [options] <prefix/of/index/files> <introns.gff3> <sorted_alignments.bam>
