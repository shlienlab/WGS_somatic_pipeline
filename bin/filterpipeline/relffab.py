# Version              Date            Developer                   Comments
#--------------------------------------------------------------------------
#     0.1                 ?                    ?        initial development

import pysam
import SequenceFetcher
import pandas as pd
import numpy as np
import warnings
import SmithWater
import os
import yaml
import logging
import datetime


class baffler():

    def __init__(self, reference_fasta, log=None, log_level=None):

        self.PYSAM_ZERO_BASED = 1

        self.sequence_alignment = {
            'gapopen': 10.0,
            'gapextend': 1.0,
            'identity_percent': 80.0,
            'similarity_percent': 80.0,
            'length_cutoff': 5
        }

        self.ratio_filter = {
            'clipped_reads': 0.8,
            'aligned_reads': 0.8
        }
        self.BAF = {
            'extend_soft_clipped_side': 20,
            'extend_mapped_side': 50,
            'max_junction_diff': 10,
        }
        self.max_coverage = {
            'coverage': 5000,
            'window': 1000
        }


        self.baffler_columns = ['clip_pos_breakpoint_1', 'aligned_breakpoint_2', 'passed_ratio_filter', 'reads_supporting_translocation', 'COV_TYPE_1', 'BAF_1_y', 'AVG_COV_1' ,'COV_STR_1' ,'clip_pos_breakpoint_2','aligned_breakpoint_1','COV_TYPE_2','BAF_2_y','AVG_COV_2','COV_STR_2']

        self.smithwater = SmithWater.SmithWater(reference_fasta)

        cigar_labels = [[0, 'M'], [1, 'I'], [2,'D'], [3,'N'], [4, 'S'], [5,'H'], [6,'P'], [7,'='], [8,'X']]
        self.cigar_df = pd.DataFrame.from_records(cigar_labels, columns=['label', 'cigar_str'])
        self.soft_clip_direction = {'5': 'left', '3': 'right'}

        if log is None: log = './baffler.tmp.'  +  str(datetime.datetime.now()).replace(' ','') + np.random.randint(0, 1000000) + '.log'
        if log_level is None: log_level = 'DEBUG'
        log_id = 'baffler.' + str(np.random.randint(0, 1000))
        self.logger = self._setup_log(log_id=log_id, logfile=log, level=log_level)

    def _setup_log(self, log_id, logfile, level=logging.INFO, mode='w'):

        logging.basicConfig(level=level)
        logger = logging.getLogger(log_id)
        logger.setLevel(level)
        handler = logging.FileHandler(logfile, mode=mode)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        return logger

    def configure(self, config, class_name='baffler'):

        self.logger.info('Configuring baffler.')

        with open(config, 'r') as f:
            options = yaml.safe_load(f)
            options = options.get(class_name, {})
            PYSAM_ZERO_BASED = options.get('PYSAM_ZERO_BASED', 1)
        
        self.PYSAM_ZERO_BASED = PYSAM_ZERO_BASED
        for option in options:
            if isinstance(options[option], dict):
                default = getattr(self, option)
                nested = options[option]
                for opt in nested.keys():
                    default_arg = default.get(opt, 'MISSING')
                    if default_arg is 'MISSING':
                        warnings.warn('Argument {arg} specified in config is not an allowed property.'.format(arg=opt))
                        continue
                    if not isinstance(nested[opt], type(default_arg)) and default_arg is not None:
                        warnings.warn('Argument {arg} specified in config is not the same type as the default: {default} = {d}, default type: {t}, passed type: {tp}'.format(arg=opt, default=opt, d=default_arg, t=type(default_arg), tp=type(nested[opt])))
                        continue
                    default[opt] = nested[opt]
                setattr(self, option, default)
            else:
                default_arg = getattr(self, option, 'MISSING')
                if default_arg is 'MISSING':
                        warnings.warn('Argument {arg} specified in config is not an allowed property.'.format(arg=opt))
                        continue
                if not isinstance(options[option], type(default_arg)) and default_arg is not None:
                    warnings.warn('Argument {arg} specified in config is not the same type as the default: {default} = {d}, default type: {t}, passed type: {tp}'.format(arg=opt, default=opt, d=default_arg, t=type(default_arg), tp=type(options[option])))
                    continue
                setattr(self, option, options[option])

        self.logger.info('Successfully configured baffler.')

    def clean_chromosomes(self, chromosomes):

        return [str(int(x)) if x in list(range(0, 23)) else str(x) for x in chromosomes ]

    def run_Baffler(self, tumor_tab, tumor_bam, merge=False, output_tab=None, passed_tumor_df=False, limited=True):

        tb_obj = pysam.AlignmentFile(tumor_bam, 'rb')

        if passed_tumor_df is True:
            df_tumor = tumor_tab.copy()
        else:
            df_tumor = pd.read_csv(tumor_tab, sep="\t", index_col=False)

        df_tumor = df_tumor.reset_index(drop=True)
        df_tumor.Chromosome1 = self.clean_chromosomes(df_tumor.Chromosome1)
        df_tumor.Chromosome2 = self.clean_chromosomes(df_tumor.Chromosome2)
        df_tumor = self.flag_excessive_coverage(df_tumor, tb_obj)

        temp = []
        for idx, x in df_tumor.iterrows():
            placeholder_record = x.to_frame().transpose().copy()
            if (x.excessive_coverage == True):
                self.logger.debug("Excessive coverage found for record index: {idx}, {chr1}:{pos1};{chr2}:{pos2}".format(idx=idx, chr1=x.Chromosome1, pos1=x.Position1, chr2=x.Chromosome2, pos2=x.Position2))
                placeholder_record = self.fill_placeholder(placeholder_record, 'EC')
                temp.append(placeholder_record)
            else:
                try:
                    # Try to calculate the translocation baf, otherwise fill the dataframe with a placeholder
                    result = self.calculate_translocation_baf(x.Chromosome1, int(x.Position1), x.Chromosome2, int(x.Position2), x.CT, tb_obj)
                    if result is not None:
                        temp.append(result)                    
                    else:
                        placeholder_record = self.fill_placeholder(placeholder_record, np.nan)
                        temp.append(placeholder_record)
                except (KeyboardInterrupt, SystemExit):
                    raise
                except:
                    self.logger.exception("Could not calculate translocation BAF for record index: {idx}, {chr1}:{pos1};{chr2}:{pos2}".format(idx=idx, chr1=x.Chromosome1, pos1=x.Position1, chr2=x.Chromosome2, pos2=x.Position2))
                    placeholder_record = self.fill_placeholder(placeholder_record, 'ERROR')
                    temp.append(placeholder_record)                    
        
        results = pd.concat(temp)
        results.Chromosome1 = self.clean_chromosomes(results.Chromosome1)
        results.Chromosome2 = self.clean_chromosomes(results.Chromosome2)

        if merge is True:
            if limited is True: results = self.clean_results(results)
            results = pd.merge(df_tumor, results, how='left', on=['Chromosome1', 'Chromosome2', 'Position1', 'Position2', 'CT'])
        
        if output_tab is not None:
            results.to_csv(output_tab, sep="\t", index=False)

        return results

    def fill_placeholder(self, placeholder, comment):

        for column in self.baffler_columns:
            placeholder[column] = comment
        return placeholder

    def flag_excessive_coverage(self, df_tumor, tb_obj):

        df_tumor['junction1_coverage'] = df_tumor.apply(lambda x: self._fetch_coverage( x.Chromosome1, x.Position1 - self.max_coverage['window'] / 2, x.Position1 + self.max_coverage['window'] / 2, tb_obj)[1] , axis=1)
        df_tumor['junction2_coverage'] = df_tumor.apply(lambda x: self._fetch_coverage( x.Chromosome2, x.Position2 - self.max_coverage['window'] / 2, x.Position2 + self.max_coverage['window'] / 2, tb_obj)[1] , axis=1)
        df_tumor['excessive_coverage'] = False
        df_tumor.ix[(df_tumor.junction1_coverage > self.max_coverage['coverage'] ) | (df_tumor.junction2_coverage > self.max_coverage['coverage']), 'excessive_coverage'] = True
        return df_tumor

    def clean_results(self, results):

        return results[['CT','Chromosome1', 'Position1', 'clip_pos_breakpoint_1', 'aligned_breakpoint_2', 'passed_ratio_filter', 'reads_supporting_translocation', 'COV_TYPE_1', 'BAF_1', 'AVG_COV_1', 'COV_STR_1', 'Chromosome2', 'Position2', 'clip_pos_breakpoint_2', 'aligned_breakpoint_1', 'COV_TYPE_2' ,'BAF_2', 'AVG_COV_2', 'COV_STR_2']]

    def apply_ratio_filter(self, new_junctions):


        new_junctions['passed_ratio_filter'] = False
        new_junctions.ix[
            ( (new_junctions.aligned_breakpoint_reads_1 * 1.0 / (new_junctions.aligned_differently_1 + new_junctions.aligned_breakpoint_reads_1) ) >= self.ratio_filter['aligned_reads'] )
            & ( (new_junctions.aligned_breakpoint_reads_2 * 1.0 / (new_junctions.aligned_differently_2 + new_junctions.aligned_breakpoint_reads_2) ) >= self.ratio_filter['aligned_reads'] )
            & ( (new_junctions.clip_pos_breakpoint_reads_1 * 1.0 / (new_junctions.clipped_differently_1 + new_junctions.clip_pos_breakpoint_reads_1) ) >= self.ratio_filter['clipped_reads'] )
            & ( (new_junctions.clip_pos_breakpoint_reads_2 * 1.0 / (new_junctions.clipped_differently_2 + new_junctions.clip_pos_breakpoint_reads_2) ) >= self.ratio_filter['clipped_reads'] )
        , 'passed_ratio_filter'] = True
        return new_junctions

    def get_soft_clipped_reads(self, junction_prime, chromosome, position, start, stop, bam_obj):
        
        splits = []

        if chromosome != 'X' and chromosome != 'Y': chromosome = int(chromosome)

        for r in bam_obj.fetch(str(chromosome), int(start - self.PYSAM_ZERO_BASED), int(stop - self.PYSAM_ZERO_BASED)):
            if r.cigartuples is None: continue

            cigar = pd.DataFrame(r.cigartuples, columns=['label','bases'])
            cigar = pd.merge(cigar, self.cigar_df, how='left', on='label')
            exploded_cigar = self.explode_cigar(r, cigar)
            exploded_cigar = self.explode_sequence_onto_vector(r, exploded_cigar)
            clipped_bases = self.clip_at_clipping_start(cigar, exploded_cigar, junction_prime, clip_type='S')
            
            if clipped_bases is None or len(clipped_bases) == 0: continue # no soft clip
            splits.append([chromosome, position, junction_prime, r.query_name, r.cigarstring, ''.join(clipped_bases.base), clipped_bases.clip_pos.head(1).item()])

        if len(splits) < 1:
            splits = [[chromosome, position, junction_prime, np.nan, np.nan, np.nan, np.nan]]

        splits_df = pd.DataFrame(splits, columns=['chromosome', 'position', 'junction_prime', 'query_name', 'cigar', 'seq', 'clip_pos'])
        return splits_df

    # def cigar_to_df(self, cigartuples):
        
    #     cigar = pd.DataFrame(cigartuples, columns=['label','bases'])
    #     cigar = pd.merge(cigar, self.cigar_df, how='left', on='label')
    #     return cigar

    def explode_cigar(self, read, cigar):

        start = read.reference_start - read.query_alignment_start + self.PYSAM_ZERO_BASED

        positions = []
        cigar_letters = []
        for idx, row in cigar.iterrows():
            cigar_str = [row.cigar_str] * row.bases
            pos = list(range(start, start + row.bases))
            positions.extend(pos)
            cigar_letters.extend(cigar_str)
            start = start + row.bases
        df = pd.DataFrame({'pos': positions, 'cig': cigar_letters})
        return df

    def clip_at_clipping_start(self, cigar, exploded_cigar, junction_prime, clip_type='S'):

        cig = exploded_cigar.copy()
        cigar['cigar_order']= ((cigar.cigar_str != cigar.cigar_str.shift(-1))).shift(1).fillna(False).cumsum()
        cig['cigar_order']= ((cig.cig != cig.cig.shift(-1))).shift(1).fillna(False).cumsum()

        if junction_prime == '5':
            # 5 Prime junction should have cigar like 50S100M***
            soft_clip = cigar[
                (cigar.cigar_str == clip_type)
                & (cigar.shift(-1).cigar_str == 'M')
            ]
            if len(soft_clip) > 0:
                soft_clip = soft_clip[soft_clip.cigar_order == soft_clip.cigar_order.min()].cigar_order

        if junction_prime == '3':
            # 3 prime junction should have cigar like ***100M50S
            soft_clip = cigar[
                (cigar.cigar_str == clip_type)
                & (cigar.shift(1).cigar_str == 'M')
            ]
            if len(soft_clip) > 0:
                soft_clip = soft_clip[soft_clip.cigar_order == soft_clip.cigar_order.max()].cigar_order

        if len(soft_clip) != 1:
            return None
        else:
            soft_clip = soft_clip.item()

        out = cig[cig.cigar_order == soft_clip]
        out = out.drop('cigar_order', 1)

        if junction_prime == '5':
            out['clip_pos'] = out.pos.max() + 1
        elif junction_prime == '3':
            out['clip_pos'] = out.pos.min() - 1
        return out

    def get_direction_for_soft_clips(self, prime):

        return self.soft_clip_direction[prime]

    def map_soft_clipped_to_reference(self, CT, prime, chromosome, position, junction_window, ref_chromosome, ref_position, ref_window, bam_obj):

        prime1, prime2 = CT.split('to')

        ref_window_left, ref_window_right = ref_window
        reference_start = ref_position - ref_window_left
        reference_end = ref_position + ref_window_right

        soft_clip_direction = self.get_direction_for_soft_clips(prime)
        soft_clipped = self.get_soft_clipped_reads(junction_prime=prime, chromosome=chromosome, position=position, start=position-junction_window, stop=position+junction_window, bam_obj = bam_obj)

        if prime1 == prime2: revcomp = (True, True)
        if prime1 != prime2: revcomp = (False, False)

        alignments = []
        for idx, row in soft_clipped.iterrows():
            if pd.isnull(row.seq): continue
            alignments.append(
                self.smithwater.align_to_reference(
                    query_name=row.query_name,
                    sequence=row.seq,
                    reference_chromosome=ref_chromosome,
                    reference_start=reference_start,
                    reference_end=reference_end,
                    reverse=revcomp[0],
                    complement=revcomp[1],
                    gapopen=self.sequence_alignment['gapopen'],
                    gapextend=self.sequence_alignment['gapextend']
                )
            )

        if len(alignments) > 0:
            alignments_df = pd.concat(alignments)
            soft_clipped = pd.merge(soft_clipped, alignments_df, how='left', on='query_name')

        return soft_clipped

    def calculate_translocation_baf(self, chromosome1, position1, chromosome2, position2, CT, bam_obj):

        position1 = int(position1)
        position2 = int(position2)

        first, second = CT.split('to')
        reads_first = self.map_soft_clipped_to_reference(
            CT=CT, 
            prime=first,
            chromosome=chromosome1,
            position=position1,
            junction_window=10,
            ref_chromosome=chromosome2,
            ref_position=position2,
            ref_window= self.build_ref_window(second),
            bam_obj=bam_obj
        )
        reads_second = self.map_soft_clipped_to_reference(
            CT=CT,
            prime=second,
            chromosome=chromosome2,
            position=position2,
            junction_window=10,
            ref_chromosome=chromosome1,
            ref_position =position1,
            ref_window= self.build_ref_window(first),
            bam_obj=bam_obj
        )

        if 'identity_percent' not in reads_first.columns and 'identity_percent' not in reads_second.columns:
            return None
        if len(reads_first.columns) > len(reads_second.columns):
            # Second junction has no soft clipped reads
            soft_clipped_reads = pd.concat([reads_first, reads_second]).reindex_axis(reads_first.columns, axis=1)
            reads_second = soft_clipped_reads[soft_clipped_reads.chromosome.isin(reads_second.chromosome)].reset_index(drop=True)
        elif len(reads_first.columns) < len(reads_second.columns):
            # First junction has no soft clipped reads
            soft_clipped_reads = pd.concat([reads_first, reads_second]).reindex_axis(reads_second.columns, axis=1)
            reads_first = soft_clipped_reads[soft_clipped_reads.chromosome.isin(reads_first.chromosome)].reset_index(drop=True)

        # Add a column to designate which side of the translocation comes first, so we can decouple them after we put them together        
        reads_first['event_order'] = 1
        reads_second['event_order'] = 2

        reads_first = self.apply_mapping_filters(reads_first)
        reads_first = self.mark_aligned_breakpoints(prime_other=second, soft_clipped=reads_first)
        reads_second = self.apply_mapping_filters(reads_second)
        reads_second = self.mark_aligned_breakpoints(prime_other=first, soft_clipped=reads_second)

        passed_filters_first = self.get_reads_passing_filters(reads_first)
        passed_filters_second = self.get_reads_passing_filters(reads_second)

        if passed_filters_first is None: passed_filters_first = reads_first
        if passed_filters_second is None: passed_filters_second = reads_second

        soft_clipped_reads = pd.concat([passed_filters_first, passed_filters_second])
        resolved_junctions = self.define_new_junctions(soft_clipped_reads)
        resolved_junctions = resolved_junctions.sort_values(by=['event_order'], ascending=True)

        if len(resolved_junctions) > 2:
            warnings.warn("More than two junctions resolved for {chr1}:{pos1};{chr2}:{pos2}. Skipping...".format(chr1=chromosome1, pos1=position1, chr2=chromosome2, pos2=position2))
            return None

        resolved_junctions = self.swap_aligned_breakpoints(resolved_junctions)

        resolved_junctions = self.calculate_baf(resolved_junctions, soft_clipped_reads, bam_obj)
        junctions_out = self.merge_junction_output(resolved_junctions)
        junctions_out = self.apply_ratio_filter(junctions_out)
        junctions_out['CT'] = CT

        return junctions_out


    def build_ref_window(self, prime):

        if prime == '5': return (self.BAF['extend_soft_clipped_side'], self.BAF['extend_mapped_side'])
        if prime == '3': return (self.BAF['extend_mapped_side'], self.BAF['extend_soft_clipped_side'])

    def get_reads_passing_filters(self, soft_clipped):

        putative_junctions = soft_clipped.copy()
        putative_junctions = putative_junctions[
            (putative_junctions.identity_filter == True) & 
            (putative_junctions.similarity_filter == True) & 
            (putative_junctions.length_filter == True)
        ]
        if len(putative_junctions) < 1: return None
        return putative_junctions

    def mark_aligned_breakpoints(self, prime_other, soft_clipped):

        soft_clipped['aligned'] = np.nan
        if prime_other == '5':
            soft_clipped['aligned'] = soft_clipped.ref_genomic_start
        elif prime_other == '3':
            soft_clipped['aligned'] = soft_clipped.ref_genomic_end

        return soft_clipped

    def define_new_junctions(self, putative_junctions):

        best_junctions = putative_junctions.groupby(['chromosome', 'position', 'event_order']).agg(
            {'clip_pos': {
                            'breakpoint': lambda x:x.value_counts().index[0] if all(pd.notnull(x)) else None,
                            'breakpoint_reads': lambda x:x.value_counts().values[0] if all(pd.notnull(x)) else None,
                            'other': lambda x:x.value_counts().values[1:].sum() if all(pd.notnull(x)) else None
                        },
             'aligned': {
                            'breakpoint': lambda x:x.value_counts().index[0] if all(pd.notnull(x)) else None,
                            'breakpoint_reads': lambda x:x.value_counts().values[0] if all(pd.notnull(x)) else None,
                            'other': lambda x:x.value_counts().values[1:].sum() if all(pd.notnull(x)) else None
                        },
            'query_name': 'count'
            }
        ).reset_index()
        best_junctions.columns = ['_'.join(x).strip('_') for x in best_junctions.columns]

        reads_supporting_translocation = len(
            putative_junctions[putative_junctions.clip_pos.isin(best_junctions.clip_pos_breakpoint)
                               & putative_junctions.aligned.isin(best_junctions.aligned_breakpoint)]
        )
        best_junctions['reads_supporting_translocation'] = reads_supporting_translocation

        best_junctions = best_junctions.rename(
            columns = {
                'query_name_count': 'total_split_reads',
                'clip_pos_other': 'clipped_differently',
                'aligned_other': 'aligned_differently'
            }
        )

        return best_junctions

    def swap_aligned_breakpoints(self, new_junctions):

        new_junctions['aligned_breakpoint_from_other'] = list(new_junctions.aligned_breakpoint)[::-1]
        return new_junctions

    def apply_mapping_filters(self, mappings_df):

        if self.sequence_alignment['identity_percent'] is not None:
            mappings_df['identity_filter'] = False
            mappings_df.ix[(mappings_df.identity_percent > self.sequence_alignment['identity_percent']), 'identity_filter'] = True
        if self.sequence_alignment['similarity_percent'] is not None:
            mappings_df['similarity_filter'] = False
            mappings_df.ix[(mappings_df.identity_percent > self.sequence_alignment['similarity_percent']), 'similarity_filter'] = True
        if self.sequence_alignment['length_cutoff'] is not None:
            mappings_df['length_filter'] = False
            mappings_df.ix[[len(x) >= self.sequence_alignment['length_cutoff'] for x in mappings_df.sequence if pd.notnull(x)], 'length_filter'] = True

        return mappings_df


    def merge_junction_output(self, new_junctions):

        new_junctions['swap_chromosome'] = list(new_junctions.chromosome)[::-1]
        new_junctions['swap_position'] = list(new_junctions.position)[::-1]
        junctions_out = pd.merge(new_junctions, new_junctions, how='left', left_on=['chromosome', 'position'], right_on=['swap_chromosome', 'swap_position'], suffixes=['_1', '_2'])
        junctions_out = junctions_out.head(1)
        junctions_out = junctions_out.drop(['event_order_1', 'event_order_2', 'swap_chromosome_1', 'swap_position_1', 'swap_chromosome_2', 'swap_position_2', 'aligned_breakpoint_from_other_1', 'aligned_breakpoint_from_other_2', 'reads_supporting_translocation_2'], 1)
        junctions_out = junctions_out.rename(columns={'chromosome_1': 'Chromosome1', 'chromosome_2': 'Chromosome2', 'position_1': 'Position1', 'position_2': 'Position2', 'reads_supporting_translocation_1': 'reads_supporting_translocation'})

        return junctions_out

    def calculate_baf(self, new_junctions, soft_clipped_reads, bam_obj):

        translocation_query_names = soft_clipped_reads[
                                soft_clipped_reads.clip_pos.isin(new_junctions.clip_pos_breakpoint)
                               & soft_clipped_reads.aligned.isin(new_junctions.aligned_breakpoint)].query_name
        for idx, row in new_junctions.iterrows():

            position1 = row.clip_pos_breakpoint
            position2 = row.aligned_breakpoint_from_other

            if pd.isnull(position1) or pd.isnull(position2) or (np.abs(position1 - position2) > self.BAF['max_junction_diff']):
                coverage_str, coverage = self._fetch_coverage(row.chromosome, row.position, row.position+1, bam_obj, reads_to_skip=translocation_query_names)
                coverage_type = 'TIER2'
            else:
                coverage_str, coverage = self._fetch_coverage(row.chromosome,  np.min([position1, position2] ), np.max([position1, position2] ), bam_obj, reads_to_skip=translocation_query_names)
                coverage_type = 'TIER1'
            new_junctions.ix[idx, 'AVG_COV'] = coverage
            new_junctions.ix[idx, 'COV_STR'] = coverage_str
            new_junctions.ix[idx, 'COV_TYPE'] = coverage_type

        new_junctions['BAF'] = np.round((new_junctions.reads_supporting_translocation * 1.0 / new_junctions.AVG_COV ).values, 4)
        return new_junctions


    def explode_sequence_onto_vector(self, read, exploded_cigar):

        sequence_out = exploded_cigar.copy()
        sequence_iterator = iter(list(read.seq))
        temp = []
        for idx, row in exploded_cigar.iterrows():
            if row.cig == 'H' or row.cig == 'D':
                temp.append(np.nan)
            else:
                try:
                    base = next(sequence_iterator)
                except StopIteration:
                    break
                temp.append(base)
        sequence_out['base'] = temp
        return sequence_out

    def _fetch_coverage(self, chromosome, position1, position2, bam_obj, reads_to_skip=[]):

        chromosome = self.clean_chromosomes([chromosome])[0]
        position1 = int(position1)
        position2 = int(position2)

        reads_to_skip = list(reads_to_skip)

        SAME_POS = False
        if position1 == position2:
            SAME_POS = True
            position2 += 1
        
        coverages = []
        for pileup in bam_obj.pileup(str(chromosome), int(position1 - self.PYSAM_ZERO_BASED), int(position2 - self.PYSAM_ZERO_BASED)):
            if ((pileup.pos + self.PYSAM_ZERO_BASED) >= position1) and ( ((pileup.pos + self.PYSAM_ZERO_BASED <= position2) and SAME_POS == False) or ((pileup.pos + self.PYSAM_ZERO_BASED < position2) and SAME_POS == True)):
                coverage_counter = len(reads_to_skip)
                if coverage_counter > 0:
                    for pileupread in pileup.pileups:
                        if pileupread.alignment.query_name in reads_to_skip: continue
                        coverage_counter += 1
                    coverages.append(coverage_counter)
                else:
                    coverages.append(pileup.n)
        return (str(coverages), np.round(np.mean(coverages), 3))


