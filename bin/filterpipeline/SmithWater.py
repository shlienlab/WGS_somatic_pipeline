# Version              Date            Developer                   Comments
#--------------------------------------------------------------------------
#     0.1                 ?                    ?        initial development

from Bio.Emboss.Applications import WaterCommandline
import re 
from Bio import AlignIO
import pandas as pd
import numpy as np
import SequenceFetcher
import sys

class SmithWater():

    def __init__(self, reference_fasta):

        self.reffetch = SequenceFetcher.SequenceFetcher(fasta_file=reference_fasta)
        self.reffetch.load_fasta()

    def align_to_reference(self, query_name, sequence, reference_chromosome, reference_start, reference_end, reverse=False, complement=False, gapopen=10.0, gapextend=1.0):

        if reference_chromosome == 'X': reference_chromosome = 23
        if reference_chromosome == 'Y': reference_chromosome = 24

        reference_sequence_df = self.reffetch.fetch_seq_df(int(reference_chromosome), int(reference_start), int(reference_end), reverse, complement)
        reference_sequence = ''.join(reference_sequence_df.seq)
        water_cline = WaterCommandline(
            asequence = 'asis:{s}'.format(s=sequence),
            bsequence = 'asis:{s}'.format(s=reference_sequence),
            gapextend = gapextend,
            gapopen = gapopen,
            stdout = True,
            auto = True,
            filter = True,
        )
        stdout, stderr = water_cline()
        alignment_df = self.parse_stdout(stdout, query_name, sequence, reference_sequence)
        (rfd, alignment_df) = self.resolve_reference_position(reference_sequence_df, alignment_df)
        return alignment_df

    def resolve_reference_position(self, rfd, alignment_df):
        
        rfd['sequence_mapped'] = np.nan
        rfd.ix[((rfd.pos) >= alignment_df.ref_start.item()) & ((rfd.pos) <= alignment_df.ref_end.item()), 'sequence_mapped'] = True
        alignment_df['ref_genomic_start'] = rfd[rfd.sequence_mapped == True].genomic_pos.min()
        alignment_df['ref_genomic_end'] = rfd[rfd.sequence_mapped == True].genomic_pos.max() 
        return (rfd, alignment_df)

    def parse_stdout(self, water_stdout, read_id, sequence, reference_sequence, a_only=True, b_only=False):

        out = re.split('\n|#', water_stdout)
        
        # Get the sequence start and stop
        asis = []
        scores = []
        for line in out:
            if 'asis' in line:
                asis.append(line)
            if 'Similarity' in line or 'Identity' in line:
                scores.append(line)
        asis = asis[-2:]
        
        alignment = []
        for i, s in enumerate(asis):
            seq = s.split()
            alignment.append( [int(seq[1]), int(seq[3]), seq[2]] )

        alignment_df = pd.DataFrame(alignment, columns=['start', 'end', 'aligned_sequence'])
        alignment_df['query_name'] = read_id
        alignment_df.ix[0, 'sequence'] = sequence
        alignment_df.ix[1, 'sequence'] = reference_sequence
        alignment_df['ref_start'] = alignment_df.start.shift(-1).astype(float)
        alignment_df['ref_end'] = alignment_df.end.shift(-1).astype(float)
        
        matches = np.nan
        possible = np.nan
        percent = np.nan
        score = np.nan
        for score in scores:
            s = [x.strip().lower() for x in score.split(':')]
            result = s[1].split()
            matches = float(result[0].split('/')[0])
            possible = float(result[0].split('/')[1])
            for q in result:
                if '%' in q:
                    percent = float(q.strip('(').strip(')').strip('%').strip())
            if s[0] == 'similarity':
                alignment_df.ix[:, 'similarity_matches'] = matches
                alignment_df.ix[:, 'similarity_sites'] = possible
                alignment_df.ix[:, 'similarity_percent'] = percent
            elif s[0] == 'identity':
                alignment_df.ix[:, 'identity_matches'] = matches
                alignment_df.ix[:, 'identity_sites'] = possible
                alignment_df.ix[:, 'identity_percent'] = percent
        if a_only is True:
            return alignment_df.ix[[0], :]
        if b_only is True:
            return alignment_df.ix[[1], :]
        else:
            return alignment_df
