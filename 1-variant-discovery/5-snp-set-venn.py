#!/usr/bin/env python3
"""Venn diagram of the SNP sets shared between variant databases.

Every input file holds one variant identifier per line. The script draws a
Venn diagram of all input sets and, for each input file, writes out the
identifiers that are private to that file.

Usage:
    python 5-snp-set-venn.py 1KGP3.SNPs.txt ChinaMAP.SNPs.txt \
        CPGDP6K.SNPs.txt dbsnp156.SNPs.txt WBBC.SNPs.txt

Output:
    venn_diagram.pdf
    <input basename>_unique.txt   for every input file
"""

import argparse

import matplotlib.pyplot as plt
from geneview import venn


def read_files_to_sets(file_paths):
    """Read every input file into a set of lines, keyed by the file name."""
    data = {}
    for file_path in file_paths:
        try:
            with open(file_path, 'r') as file:
                lines = file.read().splitlines()
                file_name = file_path.split('/')[-1].replace('.txt', '')
                data[file_name] = set(lines)
        except FileNotFoundError:
            print(f"file {file_path} not find")
    return data


def write_unique_lines(data):
    """Write the entries private to each input set to <name>_unique.txt."""
    for file_name, line_set in data.items():
        other_sets = [s for name, s in data.items() if name != file_name]
        combined_other = set().union(*other_sets)
        unique_lines = line_set - combined_other
        output_file = f"{file_name}_unique.txt"
        with open(output_file, 'w') as outfile:
            for line in unique_lines:
                outfile.write(line + '\n')


def main():
    parser = argparse.ArgumentParser(description='plot')
    parser.add_argument('files', nargs='+', help='input files')
    args = parser.parse_args()

    data = read_files_to_sets(args.files)

    if data:
        ax = venn(data)

        plt.rcParams['figure.dpi'] = 300
        plt.savefig('venn_diagram.pdf', format='pdf', bbox_inches='tight')

        write_unique_lines(data)


if __name__ == "__main__":
    main()
