#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gisaid_MI_features.py
    Get genomic features for MI positions.

@author: Bill Thompson
@license: GPL 3
@copyright: 2021_11_21
"""
import sys
import argparse
import pandas as pd
from Bio import SeqIO

def GetArgs():
    def ParseArgs(parser):
        class Parser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

        parser = Parser(description='Get features for MI positions.')
        parser.add_argument('-i', '--input_file',
                            required = True,
                            help = 'Input MI file (required). A CSV file created by gisaid_MI.py',
                            type = str)
        parser.add_argument('-o', '--output_file',
                            required = True,
                            help = 'Output CSV file (required).',
                            type = str)
        parser.add_argument('-m', '--mutation_file',
                            required = True,
                            help = 'Mutation CSV file (required). A CSV file created by gisaid_mutations.py',
                            type = str)
        parser.add_argument('-c', '--MI_cutoff',
                            required = False,
                            help = 'The lowest mutual information value to consider (default = 0.5).',
                            default = 0.5,
                            type = float)

        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def main():
    args = GetArgs()
    input_file = args.input_file
    output_file = args.output_file
    mutation_file = args.mutation_file
    mi_cutoff = args.MI_cutoff

    df_mi = pd.read_csv(input_file)
    df_mutations = pd.read_csv(mutation_file)

    df = pd.DataFrame(columns = list(df_mutations.columns))

    for idx, row in df_mi.iterrows():
        pos1 = int(row['Position_1'])
        pos2 = int(row['Position_2'])
        MI = row['MI']
        if MI < mi_cutoff:
            break

        ref_row = df_mutations.loc[pos1 == df_mutations['reference positions']]
        df = df.append(ref_row)

        ref_row = df_mutations.loc[pos2 == df_mutations['reference positions']]
        df = df.append(ref_row)

    df = df.drop_duplicates()
    df.to_csv(output_file, index = False)

if __name__ == "__main__":
    main()

