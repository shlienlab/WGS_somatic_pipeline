# Version              Date            Developer                   Comments
#--------------------------------------------------------------------------
#     0.1                 ?                    ?        initial development

import sys
import os
import pysam
import csv
import pandas as pd
import numpy as np
import yaml
import warnings

class sFilter():

    def __init__(self):

        self.PYSAM_ZERO_BASED = 1

        tier1_filter = {'tumor_reads_min': 4, 'normal_reads_max': 0}
        self.tier1_filter = tier1_filter

        tier2_filter = {'tumor_reads_min': 7, 'normal_reads_max': 0}
        self.tier2_filter = tier2_filter

        breakpoint_window = {'intrachromosomal_percent': 0.1, 'window_size': 500, 'min_separation': 1000}
        self.breakpoint_window = breakpoint_window

        self.normal_any_MAPQ = True
        self.MAPQ = None
        self.TR = None
        self.PE = None
        self.SR = None

    def configure(self, config, class_name='sFilter'):

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


    def clean_chromosomes(self, chromosomes):

        return [str(int(x)) if x in list(range(0, 23)) else str(x) for x in chromosomes ]

    def run_sFilter(self, tumor_tab, normal_bam, tumor_bam, output_tab=None, passed_tumor_df=False):

        if passed_tumor_df is True:
            df_tumor = tumor_tab.copy()
        else:
            df_tumor = pd.read_csv(tumor_tab, sep="\t", index_col=False)
        
        df_tumor.Chromosome1 = self.clean_chromosomes(df_tumor.Chromosome1)
        df_tumor.Chromosome2 = self.clean_chromosomes(df_tumor.Chromosome2)

        nb_obj = pysam.AlignmentFile(normal_bam, 'rb')
        tb_obj = pysam.AlignmentFile(tumor_bam, 'rb')        


        temp_df = pd.DataFrame(columns=[ 'passed_sFilter', 'passing_orientations', 'forward_forward_normal', 'forward_reverse_normal', 'reverse_forward_normal', 'reverse_reverse_normal', 'forward_forward_tumor', 'forward_reverse_tumor', 'reverse_forward_tumor', 'reverse_reverse_tumor'])

        temp = []
        for idx, row in df_tumor.iterrows():
            result = self.compare_tumor_to_normal(row=row, nb_obj=nb_obj, tb_obj=tb_obj)
            if result is not None: temp.append(result)
        
        if len(temp) > 0:
            sFilterResults = pd.concat(temp)
            sFilterResults = sFilterResults.reset_index(drop=True)
            results = pd.merge(df_tumor, sFilterResults, how='left', on=['Chromosome1', 'Position1', 'Chromosome2', 'Position2'])
        else:
            results = df_tumor

        if output_tab is not None:
            results.to_csv(output_tab, sep="\t", index=False)

        return results

    def prepare_output(self, result):

        row = result[['Chromosome1', 'Position1', 'Chromosome2', 'Position2']].drop_duplicates().reset_index(drop=True)

        r = pd.DataFrame(columns=['forward_forward_normal', 'forward_reverse_normal', 'reverse_forward_normal', 'reverse_reverse_normal', 'forward_forward_tumor', 'forward_reverse_tumor', 'reverse_forward_tumor', 'reverse_reverse_tumor'])
        r = r.transpose()

        result['event_label'] = result.apply( lambda x: '_'.join([x.read_1_orientation, x.read_2_orientation]), axis=1)

        tm = result[['event_label', 'number_of_reads_tumor']].set_index('event_label').transpose()
        tm.columns = [x + '_tumor' for x in tm.columns.values]
        tm = tm.reset_index(drop=True)
        nm = result[['event_label', 'number_of_reads_normal']].set_index('event_label').transpose()
        nm.columns = [x + '_normal' for x in nm.columns.values]
        nm = nm.reset_index(drop=True)
        tnm = pd.concat([tm, nm], axis=1).transpose()

        q = pd.merge(r, tnm, how='left', right_index=True, left_index=True ).transpose().fillna(0)
        q['passed_sFilter'] = np.nan
        q['passing_orientations'] = 'None'
        q['passed_sFilter'] = any(result.passed_sFilter)
        q['delly_passed_sFilter'] = result.delly_passed_sFilter.unique()[0]
        q['delly_tier'] = result.delly_tier.unique()[0]

        q.ix[(q.passed_sFilter == True), 'passing_orientations'] = ';'.join( list(result[result.passed_sFilter == True].event_label) )

        q['sFilter_intrachromosomal'] = result.intrachromosomal.unique()[0]

        result = pd.concat([row.transpose(), q.transpose()], axis=0).transpose()

        column_order = list(row.columns.values)+['passed_sFilter', 'delly_passed_sFilter', 'delly_tier', 'sFilter_intrachromosomal', 'passing_orientations' ]+list(r.index.values)

        result = result[column_order]
        return result

    def compare_tumor_to_normal(self, row, nb_obj, tb_obj):

        CT, chromosome1, position1, chromosome2, position2 = ( str(row.CT), str(row.Chromosome1), int(row.Position1), str(row.Chromosome2), int(row.Position2))

        window_size, nearby_intrachromosomal = self.check_if_nearby_intrachromomal(CT, chromosome1, position1, chromosome2, position2)

        both_directions = True
        resolve_duplicates = False
        if nearby_intrachromosomal is True:
            both_directions = True
            resolve_duplicates = True

        tumor_reads = self.fetch_from_bam(CT, chromosome1, position1, chromosome2, position2, tb_obj, window_size, both_directions=both_directions, resolve_duplicates=resolve_duplicates)
        normal_reads = self.fetch_from_bam(CT, chromosome1, position1, chromosome2, position2, nb_obj, window_size, both_directions=both_directions, resolve_duplicates=resolve_duplicates, is_normal=True)

        tumor_read_counts = self.count_reads_per_event_type(tumor_reads)
        normal_read_counts = self.count_reads_per_event_type(normal_reads)
        
        comparison = pd.merge(tumor_read_counts, normal_read_counts, how='outer', on=['is_reverse_1', 'is_reverse_2'], suffixes=['_tumor', '_normal'])
        comparison['Chromosome1'], comparison['Position1'], comparison['Chromosome2'], comparison['Position2'], comparison['intrachromosomal'] = (chromosome1, position1, chromosome2, position2, nearby_intrachromosomal)

        comparison = self.create_filter(comparison)
        comparison = self.apply_tier1_filter(comparison)
        comparison = self.apply_tier2_filter(comparison)
        comparison = self.check_delly_orientation(CT, comparison)
        comparison = self.fill_orientations(comparison)
        comparison = self.fill_total_read_counts(comparison)

        result = self.prepare_output(comparison)
        return result

    def check_delly_orientation(self, CT, comparison):

        prime1, prime2 = CT.split('to')

        delly_is_reverse_1, delly_is_reverse_2 = self._get_delly_orientation(prime1, prime2)

        delly_call = comparison[(comparison.is_reverse_1 == delly_is_reverse_1) & (comparison.is_reverse_2 == delly_is_reverse_2)]
        if len(delly_call) != 1:
            delly_passed_sFilter = np.nan
            delly_tier = np.nan
        else:
            delly_passed_sFilter = delly_call.passed_sFilter.item()
            delly_tier = delly_call.tier_sFilter.item()
        
        comparison['delly_passed_sFilter'] = delly_passed_sFilter
        comparison['delly_tier'] = delly_tier

        return comparison

    def fill_total_read_counts(self, comparison):

        totals = comparison.groupby(['Chromosome1', 'Position1', 'Chromosome2', 'Position2']).agg({'number_of_reads_tumor':  'sum', 'number_of_reads_normal': 'sum'})
        totals = totals.reset_index().rename(columns={'number_of_reads_normal': 'total_reads_normal', 'number_of_reads_tumor': 'total_reads_tumor'})
        comparison = pd.merge(comparison, totals, how='left', on=['Chromosome1', 'Position1', 'Chromosome2', 'Position2'])
        return comparison

    def fill_orientations(self, comparison):

        comparison['read_1_orientation'] = np.nan
        comparison['read_2_orientation'] = np.nan

        comparison.ix[(comparison.is_reverse_1) == True,'read_1_orientation'] = 'reverse'
        comparison.ix[(comparison.is_reverse_1) == False,'read_1_orientation'] = 'forward'

        comparison.ix[(comparison.is_reverse_2) == True,'read_2_orientation'] = 'reverse'
        comparison.ix[(comparison.is_reverse_2) == False,'read_2_orientation'] = 'forward'
        comparison = comparison.drop(['is_reverse_1', 'is_reverse_2'], 1)
        return comparison


    def create_filter(self, comparison):

        comparison['passed_sFilter'] = False
        comparison['tier_sFilter'] = np.nan
        return comparison

    def apply_tier1_filter(self, comparison):

        comparison.ix[
            (comparison.number_of_reads_tumor >= self.tier1_filter['tumor_reads_min'])
            & (comparison.number_of_reads_normal <= self.tier1_filter['normal_reads_max']), 'passed_sFilter'
        ] = True
        comparison.ix[
            (comparison.number_of_reads_tumor >= self.tier1_filter['tumor_reads_min'])
            & (comparison.number_of_reads_normal <= self.tier1_filter['normal_reads_max']), 'tier_sFilter'
        ] = 'TIER1'
        return comparison

    def apply_tier2_filter(self, comparison):

        comparison.ix[
            (comparison.number_of_reads_tumor >= self.tier2_filter['tumor_reads_min'])
            & (comparison.number_of_reads_normal <= self.tier2_filter['normal_reads_max']), 'passed_sFilter'
        ] = True
        comparison.ix[
            (comparison.number_of_reads_tumor >= self.tier2_filter['tumor_reads_min'])
            & (comparison.number_of_reads_normal <= self.tier2_filter['normal_reads_max'])
            & (comparison.tier_sFilter.isnull()), 'tier_sFilter'
        ] = 'TIER2'
        return comparison

    def count_reads_per_event_type(self, reads):

        reads_out = pd.DataFrame({'is_reverse_1': [True, True, False, False], 'is_reverse_2': [True, False, True, False]})

        if reads is not None and len(reads) > 0:
            event_counts = reads.groupby(['is_reverse_1', 'is_reverse_2']).agg({'read_name': 'nunique'})
            event_counts = event_counts.reset_index()
            event_counts = event_counts.rename(columns={'read_name': 'number_of_reads'})
        else:
            event_counts = pd.DataFrame({'is_reverse_1': [True, True, False, False], 'is_reverse_2': [True, False, True, False], 'number_of_reads' : [0,0,0,0]})
        
        reads_out = pd.merge(reads_out, event_counts, how='left', on=['is_reverse_1', 'is_reverse_2'])
        reads_out['number_of_reads'] = reads_out.number_of_reads.fillna(0)

        return reads_out

    def check_if_nearby_intrachromomal(self, CT, chr1, pos1, chr2, pos2):

        nearby_intrachromosomal = False
        if (chr1 == chr2) and (np.abs(pos1 - pos2) < self.breakpoint_window['min_separation']) and self.breakpoint_window['intrachromosomal_percent'] is not None:
            bases = int(np.abs(pos1 - pos2) * self.breakpoint_window['intrachromosomal_percent'])
            nearby_intrachromosomal = True
        else:
            bases = self.breakpoint_window['window_size']

        return (bases, nearby_intrachromosomal)

    def fetch_from_bam(self, CT, chr1, pos1, chr2, pos2, bam_obj, window_size, both_directions=False, is_normal=False, resolve_duplicates=True):


        prime1, prime2 = CT.split('to')

        pos1a, pos1b = self._create_genomic_window(prime1, pos1, window_size, both_directions)
        pos2a, pos2b = self._create_genomic_window(prime2, pos2, window_size, both_directions)

        j1 = self._fetch_from_bam(chr1, pos1, pos1a, pos1b, bam_obj, is_normal)
        j2 = self._fetch_from_bam(chr2, pos2, pos2a, pos2b, bam_obj, is_normal)

        if j1 is None or j2 is None:
            return None

        j1 = self._set_pair_of_interest(j1)
        reads = pd.merge(j1, j2, how='inner', left_on=['read_name', 'pair_number_of_interest'], right_on=['read_name', 'pair_number'], suffixes=['_1', '_2'])
        
        reads['junction_distances'] = reads.distance_to_junction_1 + reads.distance_to_junction_2
        if resolve_duplicates is True:
            # Keep any orientation being mapped to one junction, only keep the delly orientation mapped to more than one
            delly_is_reverse_1, delly_is_reverse_2 = self._get_delly_orientation(prime1, prime2)
            singletons = reads.groupby(['read_name']).size()
            singletons = singletons[singletons == 1].index.values
            reads = reads[(reads.read_name.isin(singletons) ) | ((reads.is_reverse_1 == delly_is_reverse_1) & (reads.is_reverse_2 == delly_is_reverse_2))]
            reads = reads.sort_values(['junction_distances'], ascending=True).groupby(['read_name'], as_index=True).first().reset_index()

        reads = reads.drop(['pair_number_of_interest', 'pair_number_1', 'pair_number_2', 'read_start_1', 'read_start_2', 'distance_to_junction_1', 'distance_to_junction_2', 'junction_distances'], 1)
        return reads

    def _get_delly_orientation(self, prime1, prime2):

        if prime1 == '3': delly_is_reverse_1 = False
        if prime2 == '3': delly_is_reverse_2 = False
        if prime1 == '5': delly_is_reverse_1 = True
        if prime2 == '5': delly_is_reverse_2 = True
        return (delly_is_reverse_1, delly_is_reverse_2)


    def _set_pair_of_interest(self, reads_from_first_junction):

        reads_from_first_junction['pair_number_of_interest'] = np.nan
        reads_from_first_junction.ix[(reads_from_first_junction.pair_number == 1), 'pair_number_of_interest'] = 2
        reads_from_first_junction.ix[(reads_from_first_junction.pair_number == 2), 'pair_number_of_interest'] = 1
        return reads_from_first_junction

    def _set_pair_number(self, reads):

        reads['pair_number'] = 2
        reads.ix[(reads.is_read1 == True), 'pair_number'] = 1
        reads = reads.drop('is_read1', 1)
        return reads

    def _create_genomic_window(self, junction_prime, position, bases, both_directions=False):

        if both_directions is True:
            pos1 = position - bases
            pos2 =  position + bases
        elif junction_prime == '3':
            pos1 = position - bases
            pos2 = position
        elif junction_prime == '5':
            pos1 = position
            pos2 = position + bases
        return (pos1, pos2)

    def _determine_direction(self, pos1, pos2, junction_position):

        if pos1 < junction_position and pos2 > junction_position:
            direction = 'both'
        elif pos1 == junction_position and pos2 > junction_position:
            direction = 'right'
        elif pos1 < junction_position and pos2 == junction_position:
            direction = 'left'
        else:
            direction = np.nan
        return direction

    def _fetch_from_bam(self, chromosome, junction_position, pos1, pos2, bam_obj, is_normal=False):

        junction = []

        for read in bam_obj.fetch(chromosome, pos1 + self.PYSAM_ZERO_BASED, pos2 + self.PYSAM_ZERO_BASED):
            if (1024 <= read.flag < 2048): continue
            if not is_normal:
                # Tumor reads should be MAPQ > 0
                if (read.mapping_quality == 0): continue
            elif not self.normal_any_MAPQ:
                # Check the config for whether normal reads must be MAPQ > 0
                if (read.mapping_quality == 0): continue
            if (read.is_proper_pair): continue
            if read.mapping_quality > 0:
                if 'H' in read.cigarstring: continue
            junction.append([chromosome, junction_position, read.query_name, read.is_read1, read.is_reverse, read.reference_start + self.PYSAM_ZERO_BASED, np.abs(int(junction_position) - int(read.reference_start + self.PYSAM_ZERO_BASED) )])
        
        if len(junction) > 0:
            junction_df = pd.DataFrame(junction, columns=['chromosome', 'position', 'read_name', 'is_read1', 'is_reverse', 'read_start', 'distance_to_junction'])
            junction_df = self._set_pair_number(junction_df)
            return junction_df

    def filter_chromosomes(self, df, chromosomes=None):

        # Remove anything not on the correct chromosome
        if chromosomes is None:
            chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']
        df = df[(df.Chromosome1.isin(chromosomes)) & (df.Chromosome2.isin(chromosomes))]
        return df

    def filter_quality(self, df):

        # Remove anything that does not meet passed mapping quality

        df['TR'] = int(df.PE) + 2*(int(df.SR))
        if self.PE is not None:
            df = df[df.PE > self.PE]
        if self.SR is not None:
            df = df[df.SR > self.SR]
        if self.TR is not None:
            df = df[df.TR > self.TR]
        if self.MAPQ is not None:
            df = df[df.MAPQ > self.MAPQ]
        return df
