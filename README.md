# Reconstruct protein coding exons from RNA-Seq data

This program takes a set of aligned reads (BAM format) introns (GFF3 format) format and an index of "translation blocks" from a genome sequence (FASTA format):

Depends on Python3 modules:
- BioPython
- pysam

## Before running:
The main script requires that your genome (or scaffold) sequence be indexed first. Run 'reconstruct_exons_index.py' to index your desired genome sequence (in FASTA format). Note: This only has to be done once per genome.

## To run:
reconstruct_exons.py -x path/to/index/files -i introns.gff3 -a sorted_alignments.bam

## Future work:
- Support additional input formats (GTF, BED, SAM, FASTQ)
- Support additional output formats (GTF, BED)
- Remove need to pre-generate intron file
- Reduce memory usage by not storing whole index in memory at once
