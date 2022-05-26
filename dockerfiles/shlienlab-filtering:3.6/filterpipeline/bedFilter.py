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

class bedFilter():

    def __init__(self):

        self.position1_uncertainty = 1000
        self.position2_uncertainty = 1000

    def configure(self, config, class_name='bedFilter'):

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

    def run_BedFilter(self, tumor_tab, normal_tab, output=None, passed_tumor_df=False):

        if passed_tumor_df is False:
            df_tumor = pd.read_csv(tumor_tab, sep="\t", index_col=False)
        else:
            df_tumor = tumor_tab.copy()
        df_normal = pd.read_csv(normal_tab, sep="\t", index_col=False)
        results = self.label_germline_events(df_tumor, df_normal)

        if output is not None:
            results.to_csv(output, sep="\t", index=False)
        return results

    def label_germline_events(self, df_tumor, df_normal):

        df_tumor = df_tumor.copy()
        df_tumor['BedFilter'] = df_tumor.apply(
            lambda x: self.calculate_closest_normal_event(x, df_normal), axis=1)
        return df_tumor 

    def calculate_closest_normal_event(self, x, df_normal):
        
        event_distance = df_normal[(df_normal.Chromosome1 == x.Chromosome1) & (df_normal.Chromosome2 == x.Chromosome2)].copy()
        event_distance['distance_Position1'] = np.abs(event_distance.Position1 - x.Position1)
        event_distance['distance_Position2'] = np.abs(event_distance.Position2 - x.Position2)
        germline = event_distance[((event_distance.distance_Position1 < self.position1_uncertainty) & (event_distance.distance_Position2 < self.position2_uncertainty))]
        
        if len(germline) == 0:
            return 0
        else:
            return 'FLAGGED'
