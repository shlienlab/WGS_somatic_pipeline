### filter_clips.py #############################################################################
# A python script that loads a mutation dataframe (from vcf) and its associated BAM file and creates a dataframe where all hard/soft clippings are flagged

### HISTORY #######################################################################################
# Version           Date            Developer                                Comments
#------------------------------------------------------------------------------------
#    0.01     2016-12-08              chrissy                     initial development
#    0.02     2017-04-12             rdeborja    Rscript now copied to exec directory
#    0.03     2019-07-31        Drew Thompson                   move to pipeline repo,
#                                                         bring in external functions
#    0.04     2019-12-09  Lisa-Monique Edward           allow for chromosomes with no 
#                                                                 unfiltered variants
#    1.00     2020-04-21  Lisa-Monique Edward              translate from R to Python 


import pandas as pd
import numpy as np
import os
import pysam
import sys
import re
from itertools import chain
import argparse

PYSAM_ZERO_BASED = 1
cigar_df = pd.DataFrame({
    'operation':list('MIDNSHP=X'),
    'as_int':list(range(9)),
    'name':["alignment_matches", "insertions", "deletions", "skipped_regions", "soft_clips", "hard_clips", "padding", "sequence_matches", "sequence_mismatches"]
})
cigar_df = cigar_df.set_index('operation')

def main(bam, input_file, output_file):
    muts = pd.read_csv(input_file, sep = "\t")
    bam_file = pysam.AlignmentFile(bam, "rb")
    
    filter_clips(muts, bam_file, operation = "H")
    filter_clips(muts, bam_file, offsets = 30)

    bam_file.close()
    
    muts.to_csv(output_file, sep="\t", index=False)

# prev get_cigar_mutations
def filter_clips(muts, bam_file, operation = "S", offsets = 0):
    if len(muts.index) > 0:
        output = pd.Series([0] * len(muts.index))
        # check if reads overlapping position of each mutation contains clipped reads
        for index,row in muts.iterrows():
            chrom = row.annovar_chr
            mut_pos = row.annovar_start
            alt = row.annovar_alt
            
            # get list of reads overlapping position
            reads = bam_file.fetch(contig=str(chrom), start=mut_pos-offsets-PYSAM_ZERO_BASED, end=mut_pos+offsets+1-PYSAM_ZERO_BASED)
            
            clip_count = 0
            for read in reads:
                # check if reads overlapping position of mutation contains clipped reads
                clips = get_clips(read, operation)
                if isinstance(clips, pd.DataFrame):
                    # go through each positive, check where the mutation occurs (from md_snv_position)
                    # and compare against the given position, to see if they are at the same position
                    pos_vector = count_clips(read, chrom, mut_pos, operation, clips)
                    
                    # if position in right place or indel, keep true
                    if isinstance(pos_vector, pd.DataFrame):
                        if len(pos_vector.loc[(pos_vector.position == mut_pos) & (pos_vector.mut == alt)]) > 0:
                            clip_count += 1
                    elif pos_vector == "indel":
                        clip_count += 1
            
            output[index] = clip_count
        
        # append column to initial data frame
        muts['in_' + cigar_df.loc[operation]['name']] = output
    else:
        muts['in_' + cigar_df.loc[operation]['name']] = np.nan

# prev get_cigar_mutations
def get_clips(read, operation):
    start = read.reference_start + PYSAM_ZERO_BASED
    end = start + read.query_length
    # get cigar string as tuples of (operation, length)
    cigars = pd.DataFrame(read.cigartuples, columns = ['operation', 'length'])
    
    # check if read contains clip of interest
    if cigar_df.loc[operation]['as_int'] in cigars.operation.values:
        if operation == "H":
            # for hard clips, return entire read as clipped region
            return pd.DataFrame([[start, end]], columns=['start', 'end'])
        else:
            # for soft clips, return positions of clipped regions
            # determine positions for each operation in cigar string
            cigars = cigars.assign(start = [start + sum([0] + cigars.length.values.tolist()[:i]) for i in range(len(cigars.length))])
            cigars['end'] = cigars.start + cigars.length
            return cigars[cigars.operation == cigar_df.loc[operation]['as_int']].reset_index()[['start','end']] 
    else:
        return False

# prev md_snv_position
def count_clips(read, chrom, mut_pos, operation, clips):
    if read.has_tag("MD"):
        pos_vector = pd.DataFrame(columns=["mut", "position"])
        md = read.get_tag("MD")
        
        # if no mutation in md, unmutated
        if not any((c in md) for c in "TCAG"):
            return "unmutated"
        
        # if insertion in md, indel
        if "^" in md:
            return "indel"
        
        # separate ints and strings from md
        md_l = list(filter(None, re.split(r'(\d+)', md)))
        alpha = [ e for e in md_l if not e.isnumeric() ]
        num = [ int(e) for e in md_l if e.isnumeric() ]
        
        # if deletion in md, indel
        md_len = len(list(chain.from_iterable([list(e) for e in alpha]))) + sum(num)
        if md_len != read.query_alignment_length:
            return "indel"
        
        start = read.reference_start + PYSAM_ZERO_BASED
        
        if operation != "S":
        # for hard clips, determine positions for each mutation in md string
            # trim starting matched positions
            if md_l[0].isnumeric():
                start = num[0] + start
                d = num.pop(0)
            
            for index,value in pd.Series(alpha).items():
                for i,m in enumerate(list(value)):
                    # enumerate potential adjacent mismatches
                    if index == 0:
                        pos = start + i
                    else:
                        pos = pos_vector.iloc[index - 1]['position'] + num[index - 1] + i
                    # get the alt base at the determined position using the read sequence
                    alt = read.query_sequence[pos-read.reference_start-PYSAM_ZERO_BASED]
                    
                    pos_vector = pos_vector.append(pd.DataFrame([[alt, pos]], columns=['mut', 'position']), ignore_index=True)
        else:
            # for soft clips, return positions of bases in soft clip
            # get bases for all soft clips (beginning and end) in read
            for index,clip in clips.iterrows():
                pos_vector = pos_vector.append(pd.DataFrame({'mut': list(read.query_sequence[clip.start-start:clip.end-start]), 'position': list(range(clip.start, clip.end))}))
            
            # if soft clip at the beginning of read, adjust starting base
            if read.cigartuples[0][0] == cigar_df.loc[operation]['as_int']:
                pos_vector['position'] = pos_vector.position - read.cigartuples[0][1]
        
        if len(pos_vector.index > 0):
            return pos_vector
    return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tumor_bam', required=True)
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--output_file', required=True)

    args = parser.parse_args()

    main(
        bam=args.tumor_bam,
        input_file=args.input_file,
        output_file=args.output_file
    )