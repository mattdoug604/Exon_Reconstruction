# Exon_Trap

This is a program for identifying putative protein-coding exons. Conceptually similar to the [experimental technique of the same name](https://en.wikipedia.org/wiki/Exon_trapping), it takes a set of known or predicted splice sites and a refernece genome to define exons that make up protein-coding genes (exons that do not contain a premature stop codon). 

Input:
- aligned reads (BAM format)
- introns (GFF3 format)
- index of "translation blocks" encoded by the reference genome

Requires Python 3 and the following modules:
- BioPython
- pysam

## Before running:
The main script requires that your genome sequence be indexed. Run 'exon_trap/generate_index.py' to index your desired genome sequence (in FASTA format). Note: This only has to be done once per genome.

## To run:
exon_trap [options] <prefix/of/index/files> <introns.gff3> <sorted_alignments.bam>
