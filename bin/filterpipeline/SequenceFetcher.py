# Version              Date            Developer                   Comments
#--------------------------------------------------------------------------
#     0.1                 ?                    ?        initial development

import os
from pyfaidx import Fasta
from Bio.Seq import Seq
import pandas as pd

class SequenceFetcher():

    def __init__(self, fasta_file):

        self.fasta_file = fasta_file
        if not os.path.exists(self.fasta_file):
            raise ValueError('Fasta file does not exist.')
        if not os.path.exists(self.fasta_file + '.fai'):
            raise ValueError('Fasta index file does not exist. Please ensure the .fai file is in the same directory as the Fasta file.')

    def load_fasta(self):

        if os.path.exists(self.fasta_file) & (os.path.exists(self.fasta_file + '.fai')):
            hg = Fasta(self.fasta_file)
        self.hg = hg
        return True

    def fetch_seq(self, chromosome, start, end, reverse=True, complement=True):

        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("Start and End coordinates must be integers.")
        seq = self.hg[chromosome-1][start-1:end]
        seq = seq.seq
        if reverse is True and complement is False:
            return seq[::-1]
        if reverse is False and complement is True:
            seq = Seq(seq)
            return seq.complement()
        if reverse is True and complement is True:
            seq = Seq(seq)
            return seq.reverse_complement()
        return seq


    def fetch_seq_df(self, chromosome, start, end, reverse=False, complement=False):

        if chromosome == 'X': chromosome = 23
        if chromosome == 'Y': chromosome = 24
        chromosome = int(chromosome)

        seq = self.fetch_seq(chromosome, start, end, reverse, complement)
        df = pd.DataFrame({'pos': list(range(1, len(seq)+1)), 'seq': list(seq)})
        if reverse is False:
            df['genomic_pos'] = list(range(start, end+1))
        else:
            df['genomic_pos'] = list(range(start, end+1))[::-1]
        return df