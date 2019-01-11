#!/usr/bin/env python3
"""
This script groups fasta reads by their length, and output them as new files.

Usage:

python group-fasta-readlen.py reads.fasta
"""

import argparse
from math import trunc


class Read:
    def __init__(self, id, read):
        self.id = id
        self.read = read


def cli_main():
    """Command line entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'file', type=argparse.FileType('r'), help='fasta file to be grouped')

    args = parser.parse_args()

    with args.file as file:
        it = iter(file)
        reads = {}
        for line in it:
            read_str = next(it)
            # round down to nearest thousand
            rounded_length = round(len(read_str), -3)
            read = Read(line, read_str)
            try:
                reads[rounded_length].add(read)
            except KeyError:
                reads[rounded_length] = {read}

    for k, v in reads.items():
        fn = '{}k_reads.fasta'.format(trunc(k / 1000))
        with open(fn, 'a') as f:
            for read in v:
                f.write(read.id)
                f.write(read.read)


if __name__ == '__main__':
    cli_main()
