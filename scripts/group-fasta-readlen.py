#!/usr/bin/env python3
"""
This script groups fasta reads by their length, and output them as new files.

Usage:

python group-fasta-readlen.py reads.fasta
"""

import argparse
from math import trunc


def parse_fasta(handle):
    """Iterate over Fasta records as string tuples.

    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in parse_fasta(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')

    https://github.com/biopython/biopython/blob/master/Bio/SeqIO/FastaIO.py
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    # This matches the previous implementation where .readline() was used
    for line in handle:
        if line[0] == '>':
            title = line[1:].rstrip()
            break
    else:  # no break encountered
        return  # Premature end of file, or just empty?

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == '>':
            yield title, ''.join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, ''.join(lines).replace(" ", "").replace("\r", "")


def cli_main():
    """Command line entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'file', type=argparse.FileType('r'), help='fasta file to be grouped')

    args = parser.parse_args()

    with args.file as file:
        reads = {}
        for header, seq in parse_fasta(file):
            # round down to nearest thousand
            rounded_length = round(len(seq), -3)
            read = (header, seq)
            try:
                reads[rounded_length].add(read)
            except KeyError:
                reads[rounded_length] = {read}

    for k, v in reads.items():
        fn = '{}k_reads.fasta'.format(trunc(k / 1000))
        with open(fn, 'a') as f:
            for read in v:
                f.write('>' + read[0] + '\n')
                f.write(read[1] + '\n')


if __name__ == '__main__':
    cli_main()
