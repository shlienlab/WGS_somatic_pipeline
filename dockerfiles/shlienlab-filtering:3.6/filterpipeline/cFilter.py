# Version              Date            Developer                   Comments
#--------------------------------------------------------------------------
#     0.1                 ?                    ?        initial development
#   0.1.1        2019-07-31        Drew Thompson    remove hardcoded values
#   0.1.2        2019-07-27        Drew Thompson                    bug fix
#   0.3.2        2020-05-12        Lisa-Monique Edward  Allow gatk depth of 
#                                                           coverage sample 
#                                                          interval summary
#                                                                file to be 
#                                                          found for cancer 
#                                                                     panel
#   0.3.2        2020-05-12        Scott Davidson  Force pysam fetch to use
#                                                   zero for start position
#                                                    when window size equal
#                                                            negative value


import pysam
import sys
import pandas as pd
import numpy as np
import os
import yaml
import warnings
import re
import itertools

class cFilter():

    def __init__(self, centromere_list, dustmaker_file):

        self.centromeres = pd.read_csv(centromere_list, sep="\t", index_col=False)
        
        self.PYSAM_ZERO_BASED = 1

        self.normal_coverage_threshold = {
            'window': 10000,
            'threshold_min': 10
        }
        self.high_depth_filter = {
            'window': 200,
        }
        self.unique_mapping_tumor = {
            'uniquely_mapped_ratio': 0.8
        }
        self.low_complexity = {
            'distance_to_low_complexity': 400,
            'low_complexity_min_length': 50,
            'dustmaker_file': dustmaker_file
        }
        self.multi_mapping = {
            'chromosomes': 4,
            'reads': 2,
            'window': 200
        }

    def configure(self, config, class_name='cFilter'):

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

    def run_cFilter(self, tumor_tab, tumor_bam, normal_bam, gatk_path=None, read_length=None, passed_tumor_df=False, centromeres=True, coverage_ratio=True, unique_mapping=True, high_depth=True, low_complexity=True, multi_mapping=True, output_tab=None):
        
        if high_depth is True:
            if gatk_path is not None:
                self.gatk_path = gatk_path
            elif gatk_path is None:
                candidate_gatk_path = tumor_bam.strip('.bam') + '.gatk.depthofcoverage.sample_interval_summary'
                if os.path.exists(candidate_gatk_path):
                    warnings.warn("Found gatk statistics file: {g}".format(g=candidate_gatk_path))
                    self.gatk_path = candidate_gatk_path
                else:
                    candidate_gatk_path = os.path.abspath(os.path.dirname(tumor_bam) + "/../gatkCovCalGP/" + os.path.basename(tumor_bam).split(".")[0] + '.' + os.path.basename(tumor_bam).split(".")[1] + '.genepanel.dp.sample_interval_summary')
                    if os.path.exists(candidate_gatk_path):
                        warnings.warn("Found gatk statistics file: {g}".format(g=candidate_gatk_path))
                        self.gatk_path = candidate_gatk_path
                    else:
                        self.gatk_path = None
            if self.gatk_path is None and read_length is None:
                raise ValueError('Either pass gatk_path or read_length to run the high depth filter.')
            else:
                self.read_length = read_length

        tb_obj = pysam.AlignmentFile(tumor_bam, 'rb')
        nb_obj = pysam.AlignmentFile(normal_bam, 'rb')
        
        if passed_tumor_df is True:
            df_tumor = tumor_tab
        else:
            df_tumor = pd.read_csv(tumor_tab, sep="\t", index_col=False)

        df_tumor = df_tumor.reset_index(drop=True)
        df_tumor.Chromosome1 = self.clean_chromosomes(df_tumor.Chromosome1)
        df_tumor.Chromosome2 = self.clean_chromosomes(df_tumor.Chromosome2)

        if centromeres is True:
            df_tumor['in_centromere'] = self.flag_centromeres(df_tumor, self.centromeres)
        if coverage_ratio is True:
            df_tumor['normal_coverage_threshold'] = self.flag_normal_coverage_threshold(df_tumor, nb_obj)
        if unique_mapping is True:
            df_tumor['unique_mapping'] = self.flag_unique_mapping_tumor(df_tumor, tb_obj)
        if high_depth is True:
            df_tumor['high_depth'] = self.flag_high_read_depth(df_tumor, nb_obj)
        if low_complexity is True and self.low_complexity['dustmaker_file'] is not None:
            df_tumor = self.flag_low_complexity(df_tumor)
        if multi_mapping is True:
            df_tumor['multi_mapping'] = self.flag_multi_mapping(df_tumor, nb_obj)

        if output_tab is not None:
            df_tumor.to_csv(output_tab, self="\t", index=False)

        return df_tumor
    
    def flag_multi_mapping(self, df_tumor, nb_obj, return_full=False):

        df_tumor['multi_mapping'] = '0'

        df_tumor['chromosome1_mate_count'] = df_tumor.apply(lambda x: self._count_chromosomes_of_mates(x.Chromosome1, x.Position1, nb_obj), axis=1)
        df_tumor['chromosome2_mate_count'] = df_tumor.apply(lambda x: self._count_chromosomes_of_mates(x.Chromosome2, x.Position2, nb_obj), axis=1)

        df_tumor.ix[((df_tumor.chromosome1_mate_count) > self.multi_mapping['chromosomes']) | ((df_tumor.chromosome2_mate_count) > self.multi_mapping['chromosomes']), 'multi_mapping'] = 'FLAGGED'

        if return_full is True: return df_tumor

        return df_tumor['multi_mapping']

    def _count_chromosomes_of_mates(self, chromosome, position, bam_obj):

        chromosomes_of_mates = []
        eligible_chromosomes = [str(x) for x in list(range(1, 25))]
        my_start_pos = position - self.PYSAM_ZERO_BASED - self.multi_mapping['window']
        if my_start_pos < self.PYSAM_ZERO_BASED:
            my_start_pos = self.PYSAM_ZERO_BASED
        my_end_pos = position - self.PYSAM_ZERO_BASED + self.multi_mapping['window']
        # todo have to check if end position is beyond the end of chromosome, 
        #  and then replace with chromosome end position
        for read in bam_obj.fetch(chromosome, my_start_pos, my_end_pos):
            chromosomes_of_mates.append(int(read.next_reference_id) + 1)
        mate_count = pd.DataFrame(chromosomes_of_mates, columns=['chromosome']).chromosome.value_counts().to_frame().reset_index()
        mate_count.columns = ['chromosome', 'mates']
        mate_count['pass_limit'] = mate_count.mates > self.multi_mapping['reads']
        mate_count['junction_chromosome'] = chromosome
        mate_count['junction_position'] = position
        mate_count.chromosome = mate_count.chromosome.astype(str)
        mate_count = mate_count[mate_count.chromosome.isin(eligible_chromosomes)]
        return len(mate_count[mate_count.pass_limit == True])

    def _prepare_dustmaker(self):

        dustmaker = pd.read_csv(self.low_complexity['dustmaker_file'], sep="\t", header=None, index_col=False)
        dustmaker.columns = ['chromosome', 'start', 'end']

        dustmaker.chromosome = [str(int(x)) if x in list(range(0, 23)) else str(x) for x in dustmaker.chromosome ]
        dustmaker = dustmaker[dustmaker.chromosome.isin( [str(x) for x in list(range(0, 23))] + ['X', 'Y']  ) ]
        
        return dustmaker


    def flag_low_complexity(self, df_tumor):


        dustmaker = self._prepare_dustmaker()
        df_tumor['distance_to_low_complexity_1'] = df_tumor.apply(lambda x: self._find_low_complexity(x.Chromosome1, x.Position1, dustmaker), axis=1)
        df_tumor['distance_to_low_complexity_2'] = df_tumor.apply(lambda x: self._find_low_complexity(x.Chromosome2, x.Position2, dustmaker), axis=1)
        
        return df_tumor

    def _find_low_complexity(self, chromosome, position, dustmaker):

        dusted = dustmaker[dustmaker.chromosome == chromosome]

        within_low_complexity = dusted[((dusted.start - position) < 0) & ((dusted.end - position) > 0)]
        
        if len(within_low_complexity) > 0:
            return 0

        candidates = dusted[(np.abs(dusted.start - position) < self.low_complexity['distance_to_low_complexity']) | (np.abs(dusted.end - position) < self.low_complexity['distance_to_low_complexity']) ]

        if len(candidates) > 0:
            candidates['distance_start'] = candidates.start - position
            candidates['distance_end'] = candidates.end - position

            candidates['abs_distance_start'] = np.abs(candidates.start - position)
            candidates['abs_distance_end'] = np.abs(candidates.end - position)

            candidates['nearest_end'] = candidates[['abs_distance_start', 'abs_distance_end']].min(axis=1)
            closest_segment = candidates[candidates.nearest_end == candidates.nearest_end.min()]

            if all( np.abs(closest_segment.distance_start) == closest_segment.nearest_end):
                return closest_segment.distance_start.item()
            elif all( np.abs(closest_segment.distance_end) == closest_segment.nearest_end):
                return closest_segment.distance_end.item()
            else:
                return np.nan

    def flag_centromeres(self, df_tumor, centromeres, return_full=False):

        df_tumor['centromere'] = '0'

        df_tumor['in_centromere_1'] = False
        df_tumor = pd.merge(df_tumor, centromeres, how='left', left_on=['Chromosome1'], right_on=['chrom'])
        df_tumor.ix[ (df_tumor.Position1 > df_tumor.chromStart) & (df_tumor.Position1 < df_tumor.chromEnd) , 'in_centromere_1'] = True 
        df_tumor = df_tumor.drop(['chromStart', 'chromEnd'], 1)

        df_tumor['in_centromere_2'] = False
        df_tumor = pd.merge(df_tumor, centromeres, how='left', left_on=['Chromosome2'], right_on=['chrom'])
        df_tumor.ix[ (df_tumor.Position2 > df_tumor.chromStart) & (df_tumor.Position2 < df_tumor.chromEnd) , 'in_centromere_2'] = True 
        df_tumor = df_tumor.drop(['chromStart', 'chromEnd'], 1)

        df_tumor.ix[(df_tumor.in_centromere_1 == True) | (df_tumor.in_centromere_2 == True), 'centromere'] = 'FLAGGED'

        if return_full is False:
            return df_tumor['centromere']

        return df_tumor

    def flag_normal_coverage_threshold(self, df_tumor, nb_obj, return_full=False):

        df_tumor['normal_coverage_min'] = '0'
        df_tumor['normal_coverage_1'] = df_tumor.apply(
            lambda x: self._fetch_coverage(
                chromosome=x.Chromosome1,
                start=x.Position1 - self.normal_coverage_threshold['window']/2,
                end=x.Position1 + self.normal_coverage_threshold['window']/2,
                bam_obj=nb_obj,
                remove_unmapped=True
            ), axis=1
        )

        df_tumor['normal_coverage_2'] = df_tumor.apply(
            lambda x: self._fetch_coverage(
                chromosome=x.Chromosome2,
                start=x.Position2 - self.normal_coverage_threshold['window']/2,
                end=x.Position2 + self.normal_coverage_threshold['window']/2,
                bam_obj=nb_obj,
                remove_unmapped=True
            ), axis=1
        )
        
        df_tumor.ix[(df_tumor.normal_coverage_1 < self.normal_coverage_threshold['threshold_min']) | (df_tumor.normal_coverage_2 < self.normal_coverage_threshold['threshold_min']), 'normal_coverage_min'] = 'FLAGGED'
        
        if return_full is False:
            df_tumor = df_tumor.drop(['normal_coverage_1', 'normal_coverage_2'], 1)
            return df_tumor['normal_coverage_min']
        else:
            return df_tumor

    def flag_unique_mapping_tumor(self, df_tumor, tb_obj, return_full=False):

        df_tumor['unique_mapping'] = '0'
        df_tumor['coverage_ratio_1'] = df_tumor.apply(lambda x: self._fetch_coverage(chromosome=x.Chromosome1, start=x.Position1, bam_obj=tb_obj, remove_unmapped=False), axis=1)
        df_tumor['coverage_ratio_2'] = df_tumor.apply(lambda x: self._fetch_coverage(chromosome=x.Chromosome2, start=x.Position2, bam_obj=tb_obj, remove_unmapped=False), axis=1)
        df_tumor['ucoverage_ratio_1'] = df_tumor.apply(lambda x: self._fetch_coverage(chromosome=x.Chromosome1, start=x.Position1, bam_obj=tb_obj, remove_unmapped=True), axis=1)
        df_tumor['ucoverage_ratio_2'] = df_tumor.apply(lambda x: self._fetch_coverage(chromosome=x.Chromosome2, start=x.Position2, bam_obj=tb_obj, remove_unmapped=True), axis=1)
        df_tumor.ix[((df_tumor.ucoverage_ratio_1 / df_tumor.coverage_ratio_1) < self.unique_mapping_tumor['uniquely_mapped_ratio']) | ((df_tumor.ucoverage_ratio_2 / df_tumor.coverage_ratio_2) < self.unique_mapping_tumor['uniquely_mapped_ratio']), 'unique_mapping' ] = 'FLAGGED'
        
        if return_full is False:
            df_tumor = df_tumor.drop(['coverage_ratio_1', 'coverage_ratio_2', 'ucoverage_ratio_1', 'ucoverage_ratio_2'], 1)
            return df_tumor['unique_mapping']
        else:
            return df_tumor


    def _fetch_coverage(self, chromosome, start, bam_obj, end=None, include_end=True, remove_unmapped=False):

        if end is None:
            end = start + 1
            include_end = False
        
        start = int(start)
        end = int(end)

        if start < 1:
            start = 1

        coverage = []
        for pileup in bam_obj.pileup(chromosome, start - self.PYSAM_ZERO_BASED, end - self.PYSAM_ZERO_BASED):
            if (pileup.pos >= int(start - self.PYSAM_ZERO_BASED)) and ( (pileup.pos <= int(end - self.PYSAM_ZERO_BASED) and (include_end)) or ((pileup.pos < int(end - self.PYSAM_ZERO_BASED))) ):
                if remove_unmapped is True:
                    counter = 0
                    for pileup_read in pileup.pileups:
                        if pileup_read.alignment.mapq == 0: continue
                        counter += 1
                    coverage.append(counter)
                    continue
                coverage.append(pileup.n)
        return np.mean(coverage)

    def flag_high_read_depth(self, df_tumor, nb_obj, return_full=False):

        gatk_path = getattr(self, 'gatk_path', None)
        read_length = getattr(self, 'read_length', None)
        if gatk_path is None and read_length is None:
            raise ValueError('Must either specify gatk_path or read_length.')

        df_tumor.Chromosome1 = [str(int(x)) if x in list(range(0, 23)) else x for x in df_tumor.Chromosome1]
        df_tumor.Chromosome2 = [str(int(x)) if x in list(range(0, 23)) else x for x in df_tumor.Chromosome2]

        df_tumor['high_depth'] = 0
        df_tumor['coverage_1a'] = df_tumor.apply(lambda x: self._fetch_coverage(chromosome = x.Chromosome1, start=x.Position1 - self.high_depth_filter['window'], end=x.Position1, include_end=True, bam_obj=nb_obj), axis=1 )
        df_tumor['coverage_1b'] = df_tumor.apply(lambda x: self._fetch_coverage(chromosome = x.Chromosome1, start=x.Position1, end=x.Position1 + self.high_depth_filter['window'], include_end=True, bam_obj=nb_obj), axis=1 )

        df_tumor['coverage_2a'] = df_tumor.apply(lambda x: self._fetch_coverage(chromosome = x.Chromosome2, start=x.Position2 - self.high_depth_filter['window'], end=x.Position2, include_end=True, bam_obj=nb_obj), axis=1 )
        df_tumor['coverage_2b'] = df_tumor.apply(lambda x: self._fetch_coverage(chromosome = x.Chromosome2, start=x.Position2, end=x.Position2 + self.high_depth_filter['window'], include_end=True, bam_obj=nb_obj), axis=1 )

        if gatk_path is None:
            average_depth = self._samtools_idxstats(nb_obj.filename.decode('utf-8'), read_length=self.read_length)
        else:
            average_depth = self._gatk_depth_of_coverage(self.gatk_path)

        # Filter if the depth is greater than 4 standard deviations above the mean
        average_depth['upper_limit_1'] = average_depth.average_coverage + (4 * average_depth.std_coverage)
        average_depth['upper_limit_2'] = average_depth['upper_limit_1']
        df_tumor = pd.merge(df_tumor, average_depth[['chromosome', 'upper_limit_1']], how='left', left_on=['Chromosome1'], right_on=['chromosome'])
        df_tumor['high_read_depth_1'] = (df_tumor.coverage_1a > df_tumor.upper_limit_1) | (df_tumor.coverage_1b > df_tumor.upper_limit_1)
        df_tumor = df_tumor.drop('chromosome', 1)
        df_tumor = pd.merge(df_tumor, average_depth[['chromosome', 'upper_limit_2']], how='left', left_on=['Chromosome2'], right_on=['chromosome'])
        df_tumor['high_read_depth_2'] = (df_tumor.coverage_2a > df_tumor.upper_limit_2) | (df_tumor.coverage_2b > df_tumor.upper_limit_2)

        df_tumor.ix[(df_tumor.high_read_depth_1 == True) | (df_tumor.high_read_depth_2 == True), 'high_depth'] = 'FLAGGED'

        if return_full is False:
            df_tumor = df_tumor.drop(['high_read_depth_1', 'high_read_depth_2', 'upper_limit_1', 'upper_limit_2'], 1)
            return df_tumor['high_depth']
        else:
            return df_tumor

    @staticmethod
    def _gatk_depth_of_coverage(gatk_path):
        """ Return a dataframe with columns with columns:
                - chromosome (one row per chromsome)
                - average_coverage (the average coverage on that chromosome)
                - std_coverage (the standard-deviation of coverage on that chromosome.
                    For genomes this is estimated according the poisson distribution: sqrt(average_coverage)
                    For exomes this is estimated using the empirical standard deviation)
            
            This Dataframe is used in flag_high_read_depth to filter out variants that occur in regions of exceedingly high depth.
        """
        df = pd.read_csv(gatk_path, sep="\t", index_col=False)
        chromosome = list()
        start = list()
        stop = list()
        for interval in df.Target:
            m = re.match(r'(\d{1,2}|X|Y|MT):(\d+)-(\d+)', interval)
            if m is None:
                raise ValueError("Interval {} in file {} could not be parsed".format(interval, gatk_path))
            else:
                chromosome.append(m.groups()[0])
                start.append(int(m.groups()[1]))
                stop.append(int(m.groups()[2]))
        df['chromosome'] = chromosome
        df['start'] = start
        df['stop'] = stop 
        df['lengths'] = df['stop'] - df['start']
        # Remove regions what we're captured at all
        # These regions where either deleted, or were not part of the capture kit,
        # either way they should be excluded from the average coverage
        df = df[df['average_coverage'] >= 5] 
        df['weighted_coverage'] =  df['lengths'] * df['average_coverage']
        df['weighted_coverage_squared'] = df['weighted_coverage'] * df['average_coverage']
        chr_group = df.groupby('chromosome')
        result = pd.DataFrame({'average_coverage': 
                               chr_group['weighted_coverage'].sum() / chr_group['lengths'].sum()})
        result['std_coverage'] = np.where(
                                  chr_group.size() >= 50,
                                  # Emperical standard deviation: used for exomes
                                  np.sqrt(
                                      (chr_group['weighted_coverage_squared'].sum() 
                                      / chr_group['lengths'].sum() 
                                      - result['average_coverage']**2)
                                      ),
                                  # Poisson distribution standard deviation
                                  # used for genomes exomes with limited numbers of baits
                                  np.sqrt(result['average_coverage']))
            
        for chrom in itertools.chain((str(i) for i in range(1, 23)), ('X', 'Y', 'MT')):
            if chrom not in result.index:
                # If no capture kits succeeded obtain more than five reads coverage on a chromsome,
                # then we give it mean coverage 0 and standard deviation coverage 1.
                # This arises most often with MT, or Y for females
                result.loc[chrom] = [0, 1]
        return result.reset_index()

    def _samtools_idxstats(self, bamfile, read_length, idxstats_path=None):

        if idxstats_path is None:
            idxstats_path = bamfile.strip('.bam') + '.idxstats.tab'

        os.system('samtools idxstats {bamfile} > {idxstats_path}'.format(bamfile=bamfile, idxstats_path=idxstats_path))

        df = pd.read_csv(idxstats_path, sep="\t", index_col=False, header=None)
        df.columns = ['chromosome', 'reference_length', 'mapped', 'unmapped']
        allowed_chromosomes = [str(int(x)) for x in list(range(0, 23))] + ['X', 'Y']
        df = df[df.chromosome.isin(allowed_chromosomes)]
        df['average_coverage'] = np.round(df.mapped * 1.0 * read_length / df.reference_length, 2)
        return df[['chromosome', 'average_coverage']]
