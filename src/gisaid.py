#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gisaid.py
Functions for dealing with GISAID msa and fasta files

@author: Bill Thompson
@license: GPL 3
@copyright: 2021_02_05
"""
import subprocess
import os
import psutil
from collections import Counter
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from dateutil import parser
from datetime import datetime
import numpy as np
import multiprocessing as mp
from multiprocessing import shared_memory

class GISAID_Data_Record(object):
    """
    GISAID_Data_Record - a class to hold fasta header data
    """
    def __init__(self, line, alignment_row):
        """
        Initialize the object

        Parameters
        ----------
        line : str
            A FASTA header from a BioSeq object
        alignment_row : int
            The sequence position in the FASTA file.

        Requires
        --------
        0 < alignment_row < number of sequences in filr
        """
        line = line.strip()
        seq_id, description = line.split(' ', maxsplit = 1)
        self._seq_id = seq_id[1:]
        
        description_list = description.split('|')

        # attributes
        self._date = gisaid_format_date(description_list[2])
        self._lineage = description_list[3]
        self._region = description_list[4]
        self._country = description_list[5]
        self._state = description_list[6]
        if len(description_list) > 7:
            self._city = description_list[7]
        else:
            self._city = ''

        self._alignment_row = alignment_row

    # getters for attributes
    @property
    def seq_id(self):
        return self._seq_id
        
    @property
    def date(self):
        return self._date
        
    @property
    def region(self):
        return self._region
        
    @property
    def lineage(self):
        return self._lineage
        
    @property
    def country(self):
        return self._country
    
    @property
    def state(self):
        return self._state

    @property
    def city(self):
        return self._city

    @property
    def alignment_row(self):
        return self._alignment_row
    
class GISAIDData(object):
    """
    Class for storing GISAID MSA data

    Sequence data is stored in shared memmory as a numby array of bytes.
    """
    def __init__(self, filename):
        """
        Initialize object

        Parameters
        ----------
        filename : str
            The path to the FASTA file

        Requires
        --------
        All sequences must be the same length.
        """
        def get_msa_size():
            """
            Get the sequence length and number of seqeunces in the file.

            Returns
            -------
            tuple
                seq_count - the number of sequences
                align_length - the length of the sequences
            """
            res = subprocess.check_output(['/usr/bin/grep', '>', filename]).split(b'\n')
            res.remove(b'')
            rec = next(SeqIO.parse(filename, 'fasta'))
            seq_count = len(res)
            align_length = len(rec.seq)
            
            return (seq_count, align_length)
            
        self._seq_count, self._align_length = get_msa_size()

        self._dtype = '|S1'  # For shared memory processing
        self._order = 'F'

        # data is stored in shared memory so we can use multiple processes
        self._shm = shared_memory.SharedMemory(create = True, 
                                               size = self._seq_count * self._align_length)
        self._align_array = np.ndarray((self._seq_count, self._align_length), 
                                  dtype = self._dtype, 
                                  order = self._order,
                                  buffer = self._shm.buf)
        
        self._header_info = dict()  # storage for FASTA headers

        # read the file                            
        count = 0
        curr_line = ''
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if line[0] == '>':               
                    rec = GISAID_Data_Record(line, count)
                    self._header_info[rec.seq_id] = rec
                    if len(curr_line) > 0:
                        self._align_array[count, :] = list(curr_line)
                        count += 1
                    curr_line = ''
                elif len(line) > 0:
                    curr_line += line
                    
        if len(curr_line) > 0:
            self._align_array[count, :] = list(curr_line)
            
        self._align_array.flags.writeable = False  # data is read only

    # getters for attributes
    #         
    @property
    def header_info(self):
        return self._header_info  # dictionary of headers. Key is seq_id, value is a GISAID_Data_Record object
    
    @property
    def seq_count(self):
        return self._seq_count  # number of sequences
    
    @property
    def align_length(self):     # length of each sequence
        return self._align_length
    
    @property
    def shared_memory(self):
        return self._shm        # memory buffer
    
    @property
    def dtype(self):
        return self._dtype      # for shared memory
    
    @property
    def order(self):
        return self._order      # for shared memory
    
    @property
    def align_array(self):
        return self._align_array  # sequence data, a shared numpy array of bytes
    
    def __del__(self):
        self._shm.close()
        self._shm.unlink()
    
def gisaid_read_msa(msa_file):
    """
    gisaid_read_msa - read a Multiple Sequence Alignment file

    arguments:
    msa_file - path to Multiple Sequence Alignment file

    returns:
    a BioPython align object
    """
    align = AlignIO.read(msa_file, 'fasta')
    
    return align

def gisaid_format_date(date):
    """
    format_date - a helper function. Dates from GISAID fasta files are soometimes improper.
    
    arguments:
    date - a string containing a date
    
    returns:
    a datetime date
    """
    try:
        dt = parser.parse(date)
    except ValueError:
        try:
            dt = datetime.strptime(date, '%Y-%m-00')
        except ValueError:
            dt = datetime.strptime(date, '%Y-00-00')
        
    return dt

def gisaid_read_alignment_file(filename):
    """
    Read a collection of aligned sequences in FASTA format

    arguments:
        filename - the name of the FASTA file
        
    returns:
        a GISAIDData object
    """
    align = GISAIDData(filename)

    return align

def gisaid_get_col_range(col_variation, seq_count, align_length,
                         min_quality = 0.99):
    """

    Parameters
    ----------
    align : a BioPython Bio.Align.MultipleSeqAlignment object
    min_quality : a float
        the minimum column quality
        quality is the percent of A, C, G, T in column

    Returns
    -------
    start, end : ints
        the start and ending positions of positions with sufficient quality

    """    
    start = 0
    for col in range(align_length):
        c = col_variation[col].most_common(4)
        qual = sum([x[1] for x in c if x[0] in [b'A', b'C', b'G', b'T']])/seq_count
        if qual > min_quality:
            start = col
            break
        
    end = 0
    for col in range(align_length-1, -1, -1):
        c = col_variation[col].most_common(4)
        qual = sum([x[1] for x in c if x[0] in [b'A', b'C', b'G', b'T']])/seq_count
        if qual > min_quality:
            end = col
            break

    return start, end

def print_memory_use():
    """
    Print percent of memory used.
    From https://stackoverflow.com/questions/276052/how-to-get-current-cpu-and-ram-usage-in-python

    Returns
    -------
    None.

    """    
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2.**30  # memory use in GB...I think
    print('memory use:', memoryUse)
    
def gisaid_ref_pos_to_alignment(align_data, ref_id):
    """
    Generate a map of aligned columns to reference genome positions.

    Args:
        align_data : a GISAIDData object
            Aligned data read by gisaid_read_alignment_file.

        ref_id : str
            The sequence ID of the reference sequence.

    Returns:
        a dictionary: 
            key is the aligned position, value is position in the reference sequence.
            -1 indicates a gap in the reference.
    """
    shm = shared_memory.SharedMemory(name = align_data.shared_memory.name)
    align = np.ndarray((align_data.seq_count, align_data.align_length), 
                        dtype = align_data.dtype, 
                        order = align_data.order,
                        buffer = shm.buf)
    
    seq = align[align_data.header_info[ref_id].alignment_row, :]  # reference sequence data
    pos_map = dict()
    ref_pos = 0    # 0 based
    for pos in range(len(seq)):
        if seq[pos] in [b'A', b'C', b'G', b'T']:
            pos_map[pos] = ref_pos
            ref_pos += 1
        else:
            pos_map[pos] = -1

    return pos_map

def get_col_variation(shm_name, shape, col, 
                      dtype = '|S1', order = 'F'):
    """
    Get count of each nucleotide in each column.

    Parameters
    ----------
    shm_name : str
        The id of the shared memory block.
    shape : a tuple of ints.
        (number of sequences, alignmnet width)
    col : int
        The column number.

    Returns
    -------
    TYPE
        A count of each nucleotide in the column.

    """
    shm = shared_memory.SharedMemory(name = shm_name)
    align = np.ndarray(shape, 
                        dtype = dtype, 
                        order = order,
                        buffer = shm.buf)
    
    return (col, Counter(align[:, col]))

def gisaid_get_columns_variation(alignment, 
                                 num_processes = 10):
    """
    Get counts of each nucleotide type in each column of MSA.

    Parameters
    ----------
    shm_name : str
        The assigned name of the shared memory block.
    shape : a tuple of ints
        (number of sequences, alignment width)
    num_processes : int
        The number of subprocesses to start.

    Returns
    -------
    column_variation : a list
        A list counts for each column of the MSA.

    """
    pool = mp.Pool(processes = num_processes)
    column_variation = pool.starmap(get_col_variation, [(alignment.shared_memory.name, (alignment.seq_count, alignment.align_length), col)
                                                        for col in range(alignment.align_length)])
    
    col_vary = dict()
    for col, var in column_variation:
        col_vary[col] = var
        
    return col_vary

def gisaid_get_varying_columns(align, 
                               variation,
                               consensus_cutoff = 1.0,
                               start = 130, end = 29840):
    """
    get_varying_columns - find columns in alignment that vary by more than a
                          given amount

    arguments:
    align - a Bio.Align object produced by Bio.AlignIO
    consensus_cutoff - upper value for variation in an alignment column
    start, end - starting and ending columns in the alignment

    returns:
    a dictionary - key = a column number. Thse colummns have identity <= cutoff
                   value = tuple (percent of variation, 
                                  a list containing the aligned column nucleotides)

    requires:
    consensus_cutoff, start, end >= 0
    start, end with alignment column bounds
    """                               
    variant_cols = dict()
    for col in range(start, end):
        # c = Counter(align[:, col]).most_common(1)
        c = variation[col]
        common = c.most_common(1)
        if common[0][0] in [b'A', b'C', b'G', b'T']:
            denom = sum([c[k] for k in [b'A', b'C', b'G', b'T', b'-']])
            pct = common[0][1] / denom
            if pct < consensus_cutoff:
                variant_cols[col] = (pct, list(align.align_array[:, col]))

    return variant_cols

