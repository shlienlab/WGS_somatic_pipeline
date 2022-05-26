# Version              Date            Developer                   Comments
#--------------------------------------------------------------------------
#     0.1                 ?                    ?        initial development
#   0.1.1        2019-07-31        Drew Thompson    remove hardcoded values
#   0.1.2        2019-12-09  Lisa-Monique Edward   remove dataframe copying

import argparse
import cFilter
import sFilter
import bedFilter
import relffab
import pandas as pd
import numpy as np
import os

def main(tumor_tab, tumor_bam, normal_tab, normal_bam, config, output_tab, reference, centromeres, dustmaker, gatk_path=None, run_cFilter=False, run_sFilter=False, run_bedFilter=False, run_baffler=False):
    df_tumor = pd.read_csv(tumor_tab, sep="\t", index_col=False)

    BAFFLER_LOG = output_tab.strip('.tab') + '.baffler.log'

    if len(df_tumor.index) == 0:
        if output_tab is not None:
            df_tumor.to_csv(output_tab, sep="\t", index=False)
        return df_tumor

    if run_cFilter is True:
        print("Running cFilter...")
        cf = cFilter.cFilter(centromere_list=centromeres, dustmaker_file=dustmaker)
        cf.configure(config)
        cFilter_results = cf.run_cFilter(
            tumor_tab=df_tumor,
            passed_tumor_df=True,
            tumor_bam=tumor_bam,
            normal_bam=normal_bam,
            gatk_path=gatk_path
        )
    else:
        cFilter_results = df_tumor

    if run_bedFilter is True:
        print("Running bedFilter...")
        bdf = bedFilter.bedFilter()
        bdf.configure(config)
        bedFilter_results = bdf.run_BedFilter(
            tumor_tab=cFilter_results,
            passed_tumor_df=True,
            normal_tab=normal_tab
        )
    else:
        bedFilter_results = cFilter_results

    if run_sFilter is True:
        print("Running sFilter...")
        sf = sFilter.sFilter()
        sf.configure(config)
        sFilter_results = sf.run_sFilter(
            tumor_tab=bedFilter_results,
            passed_tumor_df=True,
            tumor_bam=tumor_bam,
            normal_bam=normal_bam
        )
    else:
        sFilter_results = bedFilter_results

    if run_baffler is True:
        print("Running baffler...")
        baffler = relffab.baffler(reference_fasta=reference, log=BAFFLER_LOG)
        baffler.configure(config)
        baffler_results = baffler.run_Baffler(
            tumor_tab=sFilter_results,
            passed_tumor_df=True,
            tumor_bam=tumor_bam,
            merge=True,
            limited=True,
        )
    else:
        baffler_results = sFilter_results

    if output_tab is not None:
        baffler_results.to_csv(output_tab, sep="\t", index=False)

    return baffler_results


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tumor_tab', required=True)
    parser.add_argument('--tumor_bam', required=True)
    parser.add_argument('--normal_tab', required=True)
    parser.add_argument('--normal_bam', required=True)
    parser.add_argument('--config', required=True)
    parser.add_argument('--output_tab', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--centromeres', required=True)
    parser.add_argument('--dustmaker', required=True)

    parser.add_argument('--sFilter', required=False, action='store_true', dest='sFilter', default=False)
    parser.add_argument('--cFilter', required=False, action='store_true', dest='cFilter', default=False)
    parser.add_argument('--bedFilter', required=False, action='store_true', dest='bedFilter', default=False)
    parser.add_argument('--baffler', required=False, action='store_true', dest='baffler', default=False)
    parser.add_argument('--gatk_path', required=False, default=None)

    args = parser.parse_args()
    main(
        tumor_tab=args.tumor_tab,
        normal_tab=args.normal_tab,
        normal_bam=args.normal_bam,
        tumor_bam=args.tumor_bam,
        config=args.config,
        output_tab=args.output_tab,
        reference=args.reference,
        centromeres=args.centromeres,
        dustmaker=args.dustmaker,
        run_cFilter=args.cFilter,
        run_sFilter=args.sFilter,
        run_bedFilter=args.bedFilter,
        run_baffler=args.baffler,
        gatk_path=args.gatk_path
    )
