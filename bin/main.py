#!/home2/mattdoug/bin/python3
# Last updated: 17/11/2018
# Author: Matt Douglas

import argparse
import os
import re
import logging

log_format = '%(message)s'
logging.basicConfig(format=log_format, level=logging.INFO)
log = logging.getLogger(__name__)


def parse_region_str(string):
    """Convert a string in the format: 'I:1000..2000' or 'I_1000_2000' or
       'I:1000-2000' to a tuple.
    """
    if string is None:
        return None

    t = string.strip().replace(',', '').replace('..', '-')
    t = re.split(':|-|_', t)
    seq, start, end = t[0], int(t[1]), int(t[2])
    region = seq, start, end

    return region


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Reconstruct exons from RNA-Seq data.')
    requiredNamed = parser.add_argument_group('required input')
    # required input ###########################################################
    requiredNamed.add_argument(
        dest='index',
        type=str,
        help='prefix of the genome index files')
    requiredNamed.add_argument(
        dest='introns',
        type=str,
        help='intron in GFF3 format')
    requiredNamed.add_argument(
        dest='align',
        type=str,
        help='aligned reads in SAM or BAM format')
    # optional input ###########################################################
    parser.add_argument(
        '-o',
        dest='output',
        type=str,
        nargs='?',
        default='exons',
        help='prefix for the output files')
    parser.add_argument(
        '-r',
        dest='region',
        type=parse_region_str,
        nargs='?',
        help='limit search to genomic coordinates (e.g. I:1000-2000)')
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        default=False,
        help='do not print any information [incompatible w/ -d]')
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        default=False,
        help='print debug information')
    ############################################################################
    args = parser.parse_args()

    if not args.quiet:
        log.info('Index prefix  = {}'.format(args.index))
        log.info('Intron file   = {}'.format(args.introns))
        log.info('Aligned reads = {}'.format(args.align))
        log.info('Output prefix = {}'.format(args.output))
        if args.region:
            log.info('Search region = {}:{:,}-{:,}'.format(*args.region))
        log.info('')

    return args


if __name__ == '__main__':
    log.info('run "reconstruct_exons.py"')
