#!/home2/mattdoug/python3/bin/python3
# Last updated: 14/2/2018
# Author: Matt Douglas

def parse_arguments():
    import argparse
    import re
    from sys import exit

    parser = argparse.ArgumentParser(
                        description='Reconstruct exons from RNA-Seq data.')
    parser.add_argument('-x',
                        dest='index',
                        type=str,
                        nargs='?',
                        help='prefix of the genome index files')
    parser.add_argument('-i',
                        dest='introns',
                        type=str,
                        nargs='?',
                        help='intron in GFF3 format')
    parser.add_argument('-a',
                        dest='bam',
                        type=str,
                        nargs='?',
                        help='aligned reads in SAM/BAM format')
    parser.add_argument('-r',
                        dest='region',
                        type=str,
                        nargs='?',
                        help='limit search to genomic coordinates (e.g. I:1000-2000)')
    parser.add_argument('-p',
                        dest='prefix',
                        type=str,
                        nargs='?',
                        default='exons', help='prefix for the output files')
    parser.add_argument('-q',
                        dest='quiet',
                        action='store_true',
                        default=False,
                        help='do not print progress information [incompatible w/ -d]')
    parser.add_argument('-d',
                        dest='debug',
                        action='store_true',
                        default=False,
                        help='print debug information')
    args = parser.parse_args()

    # the required parameters are the index and the aligned reads
    if not args.index or not args.introns:
        parser.print_help()
        exit(1)

    # convert a string (e.g. "I:1000..2000" or "I_1000_2000") to a tuple
    if args.region is not None:
        temp = args.region.strip().replace(',', '').replace('..', '-')
        temp = re.split(':|-|_', temp)
        chrom, start, end = temp[0], int(temp[1]), int(temp[2])
        args.region = chrom, start, end
        region_str = '{}:{:,}-{:,}'.format(*args.region)
    else:
        region_str = 'All'

    # report arguments from the command line
    if not args.quiet:
        print('Index prefix  =', args.index)
        print('Intron file   =', args.introns)
        print('Output prefix =', args.prefix)
        print('Search region =', region_str)
        print('')

    return args
