
"""
gisaid_MI.py
    Calculate mutual information among aligned columns.

@author: Bill Thompson
@license: GPL 3
@copyright: 2021_11_19
"""
import sys
import argparse
from collections import Counter
import pandas as pd
import numpy as np
from gisaid import gisaid_ref_pos_to_alignment, \
    gisaid_read_alignment_file, \
    gisaid_get_col_range, \
    gisaid_get_columns_variation, \
    gisaid_get_varying_columns

def GetArgs():
    def ParseArgs(parser):
        class Parser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

        parser = Parser(description='Calculate mutual information among aligned columns.')
        parser.add_argument('-i', '--input_file',
                            required = True,
                            help = 'Input MSA file (required). Fasta file headers must have been reformatted with gisaid_reformat_fasta_headers.py',
                            type = str)
        parser.add_argument('-o', '--output_file',
                            required = True,
                            help = 'Output CSV file (required).',
                            type = str)
        parser.add_argument('-q', '--min_column_quality',
                            required = False,
                            help = 'Minimum quality to be accepted. Quality is the percent of A, C, G, T in column. (default = 0)',
                            default = 0.9,
                            type = float)
        parser.add_argument('-c', '--consensus_cutoff',
                            required = False,
                            help = 'The upper value for variation in an alignment column. (default = 0.98)',
                            default = 0.98,
                            type = float)
        parser.add_argument('-p', '--pseudocounts',
                            required = False,
                            help = 'Pseudocounts for MI calculations. (default = 0.05)',
                            default = 0.05,
                            type = float)
        parser.add_argument('-n', '--num_processors',
                            required = False,
                            help = 'Number of processors to use for parallel calculations. (default = 8)',
                            default = 8,
                            type = int)

        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def MI(col1, col2, pseudo = 0.05):
    """
    MI - calculate the mutual information between two columns of a 
         multisequence alignment
    
    arguments:
    col1, col2 - two columns from the alignment
    pseudo - pseudocounts
    
    returns:
    the mutual information in bits
    
    requires:
    len(col1) == len(col2)
    pseudo >= 0
    np.log2
    
    notes:
        This method is faster than using contingency_matrix from sklearn
    """
    nts = [b'A', b'C', b'G', b'T']
    
    c1 = {b'A': pseudo, b'C': pseudo, b'G': pseudo, b'T': pseudo}
    c2 = {b'A': pseudo, b'C': pseudo, b'G': pseudo, b'T': pseudo}
    p = {b'A': {b'A': pseudo, b'C': pseudo, b'G': pseudo, b'T': pseudo},
          b'C': {b'A': pseudo, b'C': pseudo, b'G': pseudo, b'T': pseudo},
          b'G': {b'A': pseudo, b'C': pseudo, b'G': pseudo, b'T': pseudo},
          b'T': {b'A': pseudo, b'C': pseudo, b'G': pseudo, b'T': pseudo}}
    
    N = 0
    for i in range(len(col1)):
        if col1[i] in nts and col2[i] in nts:
            c1[col1[i]] += 1
            c2[col2[i]] += 1
            p[col1[i]][col2[i]] += 1
            N += 1
            
    mi = 0
    for nt1 in nts:
        count1 = c1[nt1]
        for nt2 in nts:
            count2 = c2[nt2]
            mi += (p[nt1][nt2]/N) * np.log2(N * p[nt1][nt2] / (count1 * count2))
            
    return mi

def MI_table(variant_cols, pos_map, pseudo = 0.05):
    """
    MI_table - generate a table of mutual information values from all paires of
               columns in a  dictionary of aligned nucleotides
    
    arguments:
    variant_cols - a dictionary of columns. 
                   key = aligned column number, value = column data
    pos_map - a dictionary, key = aligned column position, value = reference coulum  position
    pseudo - pseudocounts
    
    returns:
    a dictionary
         'Position_1' - a list of genome positions
         'Position_2' - a list of genome positions
         'MI' - list of mutual information between Position_1 and Position_2
    
    requires:
    pseudo >= 0
    """
    var_cols = list(variant_cols.keys())
    
    mi_tab = {'Position_1': [], 'Position_2': [], 'MI': []}
    for i in range(len(var_cols)-1):
        for j in range(i+1, len(var_cols)):
            mi = MI(variant_cols[var_cols[i]][1],
                    variant_cols[var_cols[j]][1],
                    pseudo = pseudo)
            mi_tab['Position_1'].append(pos_map[var_cols[i]]+1)
            mi_tab['Position_2'].append(pos_map[var_cols[j]]+1)
            mi_tab['MI'].append(mi)
            
    return mi_tab

def main():
    args = GetArgs()
    align_file = args.input_file
    mutual_info_csv = args.output_file
    min_col_quality = args.min_column_quality
    consensus_cutoff = args.consensus_cutoff
    pseudo = args.pseudocounts
    num_processes = args.num_processors
    
    ref_id = 'EPI_ISL_402124'
    
    # get the data
    align = gisaid_read_alignment_file(align_file)

    # map alignment to the reference sequence
    pos_map = gisaid_ref_pos_to_alignment(align, ref_id)

    variation = gisaid_get_columns_variation(align, num_processes = num_processes)
    start, end = gisaid_get_col_range(variation, align.seq_count, align.align_length,
                                      min_col_quality)

   # get varying columns
    variant_cols = gisaid_get_varying_columns(align, 
                                              variation,
                                              consensus_cutoff = consensus_cutoff,
                                              start = start, end = end)

    # mutual information
    mi_tab = MI_table(variant_cols, pos_map, pseudo)
    df_mi = pd.DataFrame(mi_tab).sort_values(by = 'MI', ascending = False)
    df_mi.to_csv(mutual_info_csv, index=False)

if __name__ == "__main__":
    main()
