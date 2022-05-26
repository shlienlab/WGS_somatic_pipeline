version 1.0
import "common_tasks.wdl" as CommonTasks

workflow SSM {
    input {
        String output_dir
        String patient
        String normal_sample
        String sample
        String tumor_bam
        String tumor_bai
        String normal_bam
        String normal_bai
        File tumor_gatk_depth_of_coverage
        Array[String] intervals = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
        File reference
        File reference_index
        File reference_dict
        String build_version
        String cosmic
        String cosmic_index
        String dbsnp
        String dbsnp_index
        File annovar_humandb_tar
        File annovar_bedfile
        String cosmic_rda
        String centromeres
        String dustmaker
        String source
        String pon
        Int pon_max
        Int min_dist_complex
        Int max_number_of_flagged
        Int normal_alt_count_max

        File common_biallelic_germline_variants
        File common_biallelic_germline_variants_index

        String repo_dir
        File sa_key
        String environment
    }

    String tumor_sample = sample
    Array[String] pileup_intervals = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

    # Directories
    String log_dir = "~{output_dir}/wdl_logs+run_info"
    String progress_dir = "~{log_dir}/progress"

    # Task names
    # String GetNormalSampleNames_task_name = "~{patient}.~{sample}.normal.names"
    String Mutect2_task_name_base = "~{patient}.~{sample}.mutect"
    scatter (mutect_interval in intervals) {
        String Mutect2_task_names = "~{Mutect2_task_name_base}.~{mutect_interval}"
    }
    String MergeVCFs_Mutect2_task_name = "~{patient}.~{sample}.mergeVCFs_Mutect2"

    String GetNormalPileupSummaries_task_name_base = "~{patient}.~{sample}.normal_pileups"
    String GetTumorPileupSummaries_task_name_base = "~{patient}.~{sample}.tumor_pileups"
    scatter (pileup_interval in pileup_intervals) {
        String GetNormalPileupSummaries_task_names = "~{GetNormalPileupSummaries_task_name_base}.~{pileup_interval}"
        String GetTumorPileupSummaries_task_names = "~{GetTumorPileupSummaries_task_name_base}.~{pileup_interval}"
    }

    # GATK filtering tasks
    String GatherNormalPileupSummaries_task_name = "~{patient}.~{sample}.gather_normal_pileups"
    String GatherTumorPileupSummaries_task_name = "~{patient}.~{sample}.gather_tumor_pileups"
    String CalculateContamination_task_name = "~{patient}.~{sample}.contamination"
    String LearnReadOrientationModel_task_name = "~{patient}.~{sample}.learn_orientation"
    String MergeStats_task_name = "~{patient}.~{sample}.mergeStats"

    # String FilterMutectCalls_basic_task_name_base = "~{patient}.~{sample}.FilterMutectCalls_basic"
    # String FilterMutectCalls_OB_task_name_base = "~{patient}.~{sample}.FilterMutectCalls_OB"
    # String FilterMutectCalls_CONT_task_name_base = "~{patient}.~{sample}.FilterMutectCalls_CONT"
    String FilterMutectCalls_OB_CONT_task_name_base = "~{patient}.~{sample}.FilterMutectCalls_OB_CONT"
    String NormalizeVCF_task_name_base = "~{patient}.~{sample}.NormalizeVCF"
    String VCF2Annovar_task_name_base = "~{patient}.~{sample}.vcf2a"
    String Annovar_task_name_base = "~{patient}.~{sample}.annovar"
    String SSM_Filtering_task_name_base = "~{patient}.~{sample}.standard.filter"
    # String RDA2Tab_task_name_base = "~{patient}.~{sample}.rda2tab"
    String SNVFiltering_task_name_base = "~{patient}.~{sample}.filter"
    String SNVPostFiltering_task_name_base = "~{patient}.~{sample}.postfilter"
    String ClipFiltering_task_name_base = "~{patient}.~{sample}.clipfilter"
    scatter (filtering_interval in intervals) {
        String VCF2Annovar_task_names = "~{VCF2Annovar_task_name_base}.~{filtering_interval}"
        String Annovar_task_names = "~{Annovar_task_name_base}.~{filtering_interval}"
        String SSM_Filtering_task_names = "~{SSM_Filtering_task_name_base}.~{filtering_interval}"
        # String RDA2Tab_task_names = "~{RDA2Tab_task_name_base}.~{filtering_interval}"
        String SNVFiltering_task_names = "~{SNVFiltering_task_name_base}.~{filtering_interval}"
        String SNVPostFiltering_task_names = "~{SNVPostFiltering_task_name_base}.~{filtering_interval}"
        String ClipFiltering_task_names = "~{ClipFiltering_task_name_base}.~{filtering_interval}"
        # String FilterMutectCalls_basic_task_names = "~{FilterMutectCalls_basic_task_name_base}.~{filtering_interval}"
        # String FilterMutectCalls_OB_task_names = "~{FilterMutectCalls_OB_task_name_base}.~{filtering_interval}"
        # String FilterMutectCalls_CONT_task_names = "~{FilterMutectCalls_CONT_task_name_base}.~{filtering_interval}"
        String FilterMutectCalls_OB_CONT_task_names = "~{FilterMutectCalls_OB_CONT_task_name_base}.~{filtering_interval}"
        String NormalizeVCF_task_names = "~{NormalizeVCF_task_name_base}.~{filtering_interval}"
    }

    # String MergeVCFs_FilterMutectCalls_basic_task_name = "~{patient}.~{sample}.MergeVCFs_FilterMutectCalls_basic"
    # String MergeVCFs_FilterMutectCalls_OB_task_name = "~{patient}.~{sample}.MergeVCFs_FilterMutectCalls_OB"
    # String MergeVCFs_FilterMutectCalls_CONT_task_name = "~{patient}.~{sample}.MergeVCFs_FilterMutectCalls_CONT"
    String MergeVCFs_FilterMutectCalls_OB_CONT_task_name = "~{patient}.~{sample}.MergeVCFs_FilterMutectCalls_OB_CONT"
    String MergeVCFs_NormalizeVCF_task_name = "~{patient}.~{sample}.MergeVCFs_Normalize"

    String MergeAnnotation_task_name = "~{patient}.~{sample}.merge"

    String FilterVCF_task_name = "~{patient}.~{sample}.FilterVCF"
    String AnnotVCF_task_name = "~{patient}.~{sample}.AnnotVCF"
    String FilterGraph_task_name = "~{patient}.~{sample}.FilterGraph"
    String Cleanup_task_name = "~{patient}.~{sample}.cleanup"

    Array[String] task_names = flatten([Mutect2_task_names, [MergeVCFs_Mutect2_task_name], FilterMutectCalls_OB_CONT_task_names, NormalizeVCF_task_names, VCF2Annovar_task_names, Annovar_task_names, SSM_Filtering_task_names, SNVFiltering_task_names, SNVPostFiltering_task_names, ClipFiltering_task_names, GetNormalPileupSummaries_task_names, GetTumorPileupSummaries_task_names, [GatherNormalPileupSummaries_task_name, GatherTumorPileupSummaries_task_name, CalculateContamination_task_name, LearnReadOrientationModel_task_name, MergeStats_task_name, MergeVCFs_FilterMutectCalls_OB_CONT_task_name, MergeVCFs_NormalizeVCF_task_name, FilterVCF_task_name, AnnotVCF_task_name, FilterGraph_task_name, MergeAnnotation_task_name, Cleanup_task_name]])

    if (environment=="LOCAL") {
        call CommonTasks.FindWorkDir as SSMFindWorkDir {
            input:
                sa_key = sa_key,
                patient=patient,
                sample=sample,
                my_workflow='SSM',
                environment=environment,
                log_dir=log_dir,
                output_dir = output_dir
        }
    }

    # Task Calls and Output Filenames
    call CommonTasks.CheckProgress {
        input:
            sa_key = sa_key,
            sample = sample,
            task_names = task_names,
            progress_dir = progress_dir,
            log_dir = log_dir,
            environment = environment
    }

    # TODO add upstream task that breaks the tumor bam into 1 bam per interval (for preemptible mutect)
    scatter (index in range(length(intervals))) {
        String mutect2_interval = intervals[index]
        # captures chr1-7 and X; these will run >24h and therefore cannot be run on preemptible VMs
        Int mutect_preemptible = if (index <= 6 || index == 22) then 0 else 2

        String Mutect2_task_name = "~{Mutect2_task_name_base}.~{mutect2_interval}"
        String Mutect2_output_vcf = "~{output_dir}/~{sample}.~{mutect2_interval}.unfiltered.vcf.gz"
        String Mutect2_output_tbi = "~{Mutect2_output_vcf}.tbi"
        String Mutect2_output_stats = "~{Mutect2_output_vcf}.stats"
        String Mutect2_output_f1r2 = "~{output_dir}/~{sample}.~{mutect2_interval}.f1r2.tar.gz"
        if (CheckProgress.task_completion[Mutect2_task_name] == "false") {
            call Mutect2 {
                input:
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    tumor_name = tumor_sample,
                    normal_name = normal_sample,
                    interval = mutect2_interval,
                    log_dir = log_dir,
                    task_name = Mutect2_task_name,
                    out_vcf = Mutect2_output_vcf,
                    out_tbi = Mutect2_output_tbi,
                    out_stats = Mutect2_output_stats,
                    out_f1r2 = Mutect2_output_f1r2,
                    preemptible = mutect_preemptible,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }
    }

    String MergeVCFs_Mutect2_output_vcf = "~{output_dir}/~{sample}.mutect2.vcf.gz"
    String MergeVCFs_Mutect2_output_idx = "~{MergeVCFs_Mutect2_output_vcf}.tbi"
    if (CheckProgress.task_completion[MergeVCFs_Mutect2_task_name] == "false") {
        Array [String] MergeVCFs_Mutect2_inputs = if (length(select_all(Mutect2.output_vcf)) == length(intervals)) then select_all(Mutect2.output_vcf) else Mutect2_output_vcf
        call MergeVCFs as MergeVCFs_Mutect2 {
            input:
                input_vcfs = MergeVCFs_Mutect2_inputs,
                log_dir = log_dir,
                task_name = MergeVCFs_Mutect2_task_name,
                out_vcf = MergeVCFs_Mutect2_output_vcf,
                out_idx = MergeVCFs_Mutect2_output_idx,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String MergeStats_output_stats = "~{output_dir}/~{sample}.merged.stats"
    if (CheckProgress.task_completion[MergeStats_task_name] == "false") {
        Array[String] MergeStats_inputs = if (length(select_all(Mutect2.output_stats)) == length(intervals)) then select_all(Mutect2.output_stats) else Mutect2_output_stats
        call MergeStats {
            input:
                input_stats = MergeStats_inputs,
                log_dir = log_dir,
                task_name = MergeStats_task_name,
                out_stats = MergeStats_output_stats,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    scatter (pileup_interval in pileup_intervals) {
        String GetNormalPileupSummaries_task_name = "~{GetNormalPileupSummaries_task_name_base}.~{pileup_interval}"
        String GetNormalPileupSummaries_output_table = "~{output_dir}/~{sample}.~{pileup_interval}.normal.pileups.table"
        if (CheckProgress.task_completion[GetNormalPileupSummaries_task_name] == "false") {
            call GetPileupSummaries as GetNormalPileupSummaries {
                input:
                    input_bam = normal_bam,
                    input_bai = normal_bai,
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    common_biallelic_germline_variants = common_biallelic_germline_variants,
                    common_biallelic_germline_variants_index = common_biallelic_germline_variants_index,
                    interval = pileup_interval,
                    log_dir = log_dir,
                    task_name = GetNormalPileupSummaries_task_name,
                    out_table = GetNormalPileupSummaries_output_table,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String GetTumorPileupSummaries_task_name = "~{GetTumorPileupSummaries_task_name_base}.~{pileup_interval}"
        String GetTumorPileupSummaries_output_table = "~{output_dir}/~{sample}.~{pileup_interval}.tumor.pileups.table"
        if (CheckProgress.task_completion[GetTumorPileupSummaries_task_name] == "false") {
            call GetPileupSummaries as GetTumorPileupSummaries {
                input:
                    input_bam = tumor_bam,
                    input_bai = tumor_bai,
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    common_biallelic_germline_variants = common_biallelic_germline_variants,
                    common_biallelic_germline_variants_index = common_biallelic_germline_variants_index,
                    interval = pileup_interval,
                    log_dir = log_dir,
                    task_name = GetTumorPileupSummaries_task_name,
                    out_table = GetTumorPileupSummaries_output_table,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }
    }

    String GatherNormalPileupSummaries_output_table = "~{output_dir}/~{sample}.normal.pileups.table"
    if (CheckProgress.task_completion[GatherNormalPileupSummaries_task_name] == "false") {
        Array[String] GatherNormalPileupSummaries_inputs = if (length(select_all(GetNormalPileupSummaries.output_table)) == length(pileup_intervals)) then select_all(GetNormalPileupSummaries.output_table) else GetNormalPileupSummaries_output_table
        call GatherPileupSummaries as GatherNormalPileupSummaries {
            input:
                input_tables = GatherNormalPileupSummaries_inputs,
                reference_dict = reference_dict,
                log_dir = log_dir,
                task_name = GatherNormalPileupSummaries_task_name,
                out_table = GatherNormalPileupSummaries_output_table,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String GatherTumorPileupSummaries_output_table = "~{output_dir}/~{sample}.tumor.pileups.table"
    if (CheckProgress.task_completion[GatherTumorPileupSummaries_task_name] == "false") {
        Array[String] GatherTumorPileupSummaries_inputs = if (length(select_all(GetTumorPileupSummaries.output_table)) == length(pileup_intervals)) then select_all(GetTumorPileupSummaries.output_table) else GetTumorPileupSummaries_output_table
        call GatherPileupSummaries as GatherTumorPileupSummaries {
            input:
                input_tables = GatherTumorPileupSummaries_inputs,
                reference_dict = reference_dict,
                log_dir = log_dir,
                task_name = GatherTumorPileupSummaries_task_name,
                out_table = GatherTumorPileupSummaries_output_table,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String CalculateContamination_output_contamination = "~{output_dir}/~{sample}.contamination.table"
    String CalculateContamination_output_segmentation = "~{output_dir}/~{sample}.tumor-segments.table"
    if (CheckProgress.task_completion[CalculateContamination_task_name] == "false") {
        call CalculateContamination {
            input:
                normal_pileups = select_first([GatherNormalPileupSummaries.output_table, GatherNormalPileupSummaries_output_table]),
                tumor_pileups = select_first([GatherTumorPileupSummaries.output_table, GatherTumorPileupSummaries_output_table]),
                log_dir = log_dir,
                task_name = CalculateContamination_task_name,
                out_contamination = CalculateContamination_output_contamination,
                out_segmentation = CalculateContamination_output_segmentation,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String LearnReadOrientationModel_output_artifact_priors = "~{output_dir}/~{sample}.artifact_priors.tar.gz"
    if (CheckProgress.task_completion[LearnReadOrientationModel_task_name] == "false") {
        Array[String] LearnReadOrientationModel_inputs = if (length(select_all(Mutect2.output_f1r2)) == length(intervals)) then select_all(Mutect2.output_f1r2) else Mutect2_output_f1r2
        call LearnReadOrientationModel {
            input:
                input_f1r2s = LearnReadOrientationModel_inputs,
                log_dir = log_dir,
                task_name = LearnReadOrientationModel_task_name,
                out_artifact_priors = LearnReadOrientationModel_output_artifact_priors,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    scatter (index in range(length(intervals))) {
        String interval = intervals[index]

        # String FilterMutectCalls_basic_task_name = "~{FilterMutectCalls_basic_task_name_base}.~{interval}"
        # String FilterMutectCalls_basic_output_vcf = "~{output_dir}/~{sample}.basic_filters_applied.~{interval}.vcf.gz"
        # String FilterMutectCalls_basic_output_idx = "~{FilterMutectCalls_basic_output_vcf}.tbi"
        # String FilterMutectCalls_basic_output_filtering_stats = "~{output_dir}/~{sample}.basic_filtering.~{interval}.stats"
        # if (CheckProgress.task_completion[FilterMutectCalls_basic_task_name] == "false") {
        #     call FilterMutectCalls as FilterMutectCalls_basic {
        #         input:
        #             reference = reference,
        #             reference_index = reference_index,
        #             reference_dict = reference_dict,
        #             input_vcf = select_first([MergeVCFs_Mutect2.output_vcf, MergeVCFs_Mutect2_output_vcf]),
        #             input_vcf_index = select_first([MergeVCFs_Mutect2.output_idx, MergeVCFs_Mutect2_output_idx]),
        #             mutect_stats = select_first([MergeStats.output_stats, MergeStats_output_stats]),
        #             filter_ob = false,
        #             filter_cont = false,
        #             interval = interval,
        #             log_dir = log_dir,
        #             task_name = FilterMutectCalls_basic_task_name,
        #             out_vcf = FilterMutectCalls_basic_output_vcf,
        #             out_idx =  FilterMutectCalls_basic_output_idx,
        #             out_stats = FilterMutectCalls_basic_output_filtering_stats,
        #             environment = environment,
        #             sa_key = sa_key,
        #             progress_dir = progress_dir,
        #             output_dir = output_dir
        #     }
        # }

        # String FilterMutectCalls_OB_task_name = "~{FilterMutectCalls_OB_task_name_base}.~{interval}"
        # String FilterMutectCalls_OB_output_vcf = "~{output_dir}/~{sample}.OB_filters_applied.~{interval}.vcf.gz"
        # String FilterMutectCalls_OB_output_idx = "~{FilterMutectCalls_OB_output_vcf}.tbi"
        # String FilterMutectCalls_OB_output_filtering_stats = "~{output_dir}/~{sample}.OB_filtering.~{interval}.stats"
        # if (CheckProgress.task_completion[FilterMutectCalls_OB_task_name] == "false") {
        #     call FilterMutectCalls as FilterMutectCalls_OB {
        #         input:
        #             reference = reference,
        #             reference_index = reference_index,
        #             reference_dict = reference_dict,
        #             input_vcf = select_first([MergeVCFs_Mutect2.output_vcf, MergeVCFs_Mutect2_output_vcf]),
        #             input_vcf_index = select_first([MergeVCFs_Mutect2.output_idx, MergeVCFs_Mutect2_output_idx]),
        #             mutect_stats = select_first([MergeStats.output_stats, MergeStats_output_stats]),
        #             filter_ob = true,
        #             filter_cont = false,
        #             artifact_priors = select_first([LearnReadOrientationModel.output_artifact_priors, LearnReadOrientationModel_output_artifact_priors]),
        #             interval = interval,
        #             log_dir = log_dir,
        #             task_name = FilterMutectCalls_OB_task_name,
        #             out_vcf = FilterMutectCalls_OB_output_vcf,
        #             out_idx = FilterMutectCalls_OB_output_idx,
        #             out_stats = FilterMutectCalls_OB_output_filtering_stats,
        #             environment = environment,
        #             sa_key = sa_key,
        #             progress_dir = progress_dir,
        #             output_dir = output_dir
        #     }
        # }

        # String FilterMutectCalls_CONT_task_name = "~{FilterMutectCalls_CONT_task_name_base}.~{interval}"
        # String FilterMutectCalls_CONT_output_vcf = "~{output_dir}/~{sample}.CONT_filters_applied.~{interval}.vcf.gz"
        # String FilterMutectCalls_CONT_output_idx = "~{FilterMutectCalls_CONT_output_vcf}.tbi"
        # String FilterMutectCalls_CONT_output_filtering_stats = "~{output_dir}/~{sample}.CONT_filtering.~{interval}.stats"
        # if (CheckProgress.task_completion[FilterMutectCalls_CONT_task_name] == "false") {
        #     call FilterMutectCalls as FilterMutectCalls_CONT {
        #         input:
        #             reference = reference,
        #             reference_index = reference_index,
        #             reference_dict = reference_dict,
        #             input_vcf = select_first([MergeVCFs_Mutect2.output_vcf, MergeVCFs_Mutect2_output_vcf]),
        #             input_vcf_index = select_first([MergeVCFs_Mutect2.output_idx, MergeVCFs_Mutect2_output_idx]),
        #             mutect_stats = select_first([MergeStats.output_stats, MergeStats_output_stats]),
        #             filter_ob = false,
        #             filter_cont = true,
        #             contamination_table = select_first([CalculateContamination.output_contamination, CalculateContamination_output_contamination]),
        #             tumor_segmentation_table = select_first([CalculateContamination.output_segmentation, CalculateContamination_output_segmentation]),
        #             interval = interval,
        #             log_dir = log_dir,
        #             task_name = FilterMutectCalls_CONT_task_name,
        #             out_vcf = FilterMutectCalls_CONT_output_vcf,
        #             out_idx = FilterMutectCalls_CONT_output_idx,
        #             out_stats = FilterMutectCalls_CONT_output_filtering_stats,
        #             environment = environment,
        #             sa_key = sa_key,
        #             progress_dir = progress_dir,
        #             output_dir = output_dir
        #     }
        # }

        String FilterMutectCalls_OB_CONT_task_name = "~{FilterMutectCalls_OB_CONT_task_name_base}.~{interval}"
        String FilterMutectCalls_OB_CONT_output_vcf = "~{output_dir}/~{sample}.OB_CONT_filters_applied.~{interval}.vcf.gz"
        String FilterMutectCalls_OB_CONT_output_idx = "~{FilterMutectCalls_OB_CONT_output_vcf}.tbi"
        String FilterMutectCalls_OB_CONT_output_filtering_stats = "~{output_dir}/~{sample}.OB_CONT_filtering.~{interval}.stats"
        if (CheckProgress.task_completion[FilterMutectCalls_OB_CONT_task_name] == "false") {
            call FilterMutectCalls as FilterMutectCalls_OB_CONT {
                input:
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    input_vcf = select_first([MergeVCFs_Mutect2.output_vcf, MergeVCFs_Mutect2_output_vcf]),
                    input_vcf_index = select_first([MergeVCFs_Mutect2.output_idx, MergeVCFs_Mutect2_output_idx]),
                    mutect_stats = select_first([MergeStats.output_stats, MergeStats_output_stats]),
                    filter_ob = true,
                    filter_cont = true,
                    artifact_priors = select_first([LearnReadOrientationModel.output_artifact_priors, LearnReadOrientationModel_output_artifact_priors]),
                    contamination_table = select_first([CalculateContamination.output_contamination, CalculateContamination_output_contamination]),
                    tumor_segmentation_table = select_first([CalculateContamination.output_segmentation, CalculateContamination_output_segmentation]),
                    interval = interval,
                    log_dir = log_dir,
                    task_name = FilterMutectCalls_OB_CONT_task_name,
                    out_vcf = FilterMutectCalls_OB_CONT_output_vcf,
                    out_idx = FilterMutectCalls_OB_CONT_output_idx,
                    out_stats = FilterMutectCalls_OB_CONT_output_filtering_stats,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String NormalizeVCF_task_name = "~{NormalizeVCF_task_name_base}.~{interval}"
        String NormalizeVCF_output_vcf = "~{output_dir}/~{sample}.OB_CONT_filters_applied.normalized.~{interval}.vcf.gz"
        String NormalizeVCF_output_idx = "~{NormalizeVCF_output_vcf}.tbi"
        if (CheckProgress.task_completion[NormalizeVCF_task_name] == "false") {
            call NormalizeVCF {
                input:
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    input_vcf = select_first([FilterMutectCalls_OB_CONT.output_vcf,FilterMutectCalls_OB_CONT_output_vcf]),
                    log_dir = log_dir,
                    task_name = NormalizeVCF_task_name,
                    out_vcf = NormalizeVCF_output_vcf,
                    out_idx = NormalizeVCF_output_idx,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String VCF2Annovar_task_name = "~{VCF2Annovar_task_name_base}.~{interval}"
        String VCF2Annovar_output_file = "~{output_dir}/~{sample}.~{interval}.annovar"
        if (CheckProgress.task_completion[VCF2Annovar_task_name] == "false") {
            call VCF2Annovar {
                input:
                    vcf = select_first([NormalizeVCF.output_vcf, NormalizeVCF_output_vcf]),
                    idx = select_first([NormalizeVCF.output_idx, NormalizeVCF_output_idx]),
                    tumor_name = tumor_sample,
                    normal_name = normal_sample,
                    log_dir = log_dir,
                    task_name = VCF2Annovar_task_name,
                    out_file = VCF2Annovar_output_file,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String Annovar_task_name = "~{Annovar_task_name_base}.~{interval}"
        String Annovar_output_file = "~{VCF2Annovar_output_file}.hg19_multianno.txt"
        if (CheckProgress.task_completion[Annovar_task_name] == "false") {
            call Annovar {
                input:
                    input_file = select_first([VCF2Annovar.output_file, VCF2Annovar_output_file]),
                    humandb_tar = annovar_humandb_tar,
                    bedfile = annovar_bedfile,
                    build_version = build_version,
                    log_dir = log_dir,
                    task_name = Annovar_task_name,
                    out_prefix = VCF2Annovar_output_file,
                    out_file = Annovar_output_file,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String SSM_Filtering_task_name = "~{SSM_Filtering_task_name_base}.~{interval}"
        String SSM_Filtering_output_snvindel_file = "~{output_dir}/~{sample}_annotated.~{interval}.txt"
        String SSM_Filtering_output_snv_file = "~{output_dir}/~{sample}_annotated_filtered_snv.~{interval}.txt"
        String SSM_Filtering_output_pype_file = "~{output_dir}/~{sample}_annotated_filtered_snv.pype.~{interval}.txt"
        String SSM_Filtering_output_indel_file = sub(SSM_Filtering_output_snv_file, "_snv.~{interval}.txt$", "_indel.~{interval}.txt")
        String SSM_Filtering_output_filter_info = "~{output_dir}/~{sample}_filter_info.~{interval}.txt"
        String SSM_Filtering_output_annot_file = "~{output_dir}/~{sample}_annots.~{interval}.txt"
        if (CheckProgress.task_completion[SSM_Filtering_task_name] == "false") {
            call SSM_Filtering {
                input:
                    annotations = select_first([Annovar.output_file, Annovar_output_file]),
                    vcf = select_first([NormalizeVCF.output_vcf, NormalizeVCF_output_vcf]),
                    idx = select_first([NormalizeVCF.output_idx, NormalizeVCF_output_idx]),
                    tumor_name = tumor_sample,
                    normal_name = normal_sample,
                    source = source,
                    cosmic_rda = cosmic_rda,
                    log_dir = log_dir,
                    task_name = SSM_Filtering_task_name,
                    out_snvindel_file = SSM_Filtering_output_snvindel_file,
                    out_snv_file = SSM_Filtering_output_snv_file,
                    out_pype_file = SSM_Filtering_output_pype_file,
                    out_indel_file = SSM_Filtering_output_indel_file,
                    filter_info = SSM_Filtering_output_filter_info,
                    annot_file = SSM_Filtering_output_annot_file,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        # String RDA2Tab_task_name = "~{RDA2Tab_task_name_base}.~{interval}"
        # String RDA2Tab_output_tab = "~{output_dir}/~{sample}_annotated_filtered_snv.pype.~{interval}.txt"
        # if (CheckProgress.task_completion[RDA2Tab_task_name] == "false") {
        #     call RDA2Tab {
        #         input:
        #             input_file = select_first([SSM_Filtering.output_snv_rda, SSM_Filtering_output_snv_rda]),
        #             log_dir = log_dir,
        #             task_name = RDA2Tab_task_name,
        #             out_tab = RDA2Tab_output_tab,
        #             environment = environment,
        #             sa_key = sa_key,
        #             progress_dir = progress_dir,
        #             output_dir = output_dir
        #     }
        # }

        String SNVFiltering_task_name = "~{SNVFiltering_task_name_base}.~{interval}"
        String SNVFiltering_output_tab = "~{output_dir}/~{sample}_cfilter.~{interval}.txt"
        if (CheckProgress.task_completion[SNVFiltering_task_name] == "false") {
            call SNVFiltering {
                input:
                    input_file = select_first([SSM_Filtering.output_pype_file, SSM_Filtering_output_pype_file]),
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    centromeres = centromeres,
                    dustmaker = dustmaker,
                    gatk_path = tumor_gatk_depth_of_coverage,
                    log_dir = log_dir,
                    task_name = SNVFiltering_task_name,
                    out_tab = SNVFiltering_output_tab,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String SNVPostFiltering_task_name = "~{SNVPostFiltering_task_name_base}.~{interval}"
        String SNVPostFiltering_output_txt = "~{output_dir}/~{sample}_cfilter_pon.~{interval}.txt"
        if (CheckProgress.task_completion[SNVPostFiltering_task_name] == "false") {
            call SNVPostFiltering {
                input:
                    cfilter_input = select_first([SNVFiltering.output_tab, SNVFiltering_output_tab]),
                    snv_input = select_first([SSM_Filtering.output_snv_file, SSM_Filtering_output_snv_file]),
                    pon = pon,
                    pon_max = pon_max,
                    min_dist_complex = min_dist_complex,
                    max_number_of_flagged = max_number_of_flagged,
                    normal_alt_count_max = normal_alt_count_max,
                    log_dir = log_dir,
                    task_name = SNVPostFiltering_task_name,
                    out_file = SNVPostFiltering_output_txt,
                    filter_info = select_first([SSM_Filtering.output_filter_info, SSM_Filtering_output_filter_info]),
                    annot_file = select_first([SSM_Filtering.output_annot_file, SSM_Filtering_output_annot_file]),
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String ClipFiltering_task_name = "~{ClipFiltering_task_name_base}.~{interval}"
        String ClipFiltering_output_txt = "~{output_dir}/~{sample}_annotated_filtered_clipped.~{interval}.txt"
        if (CheckProgress.task_completion[ClipFiltering_task_name] == "false") {
            call ClipFiltering {
                input:
                    input_file = select_first([SNVPostFiltering.output_txt, SNVPostFiltering_output_txt]),
                    bam_file = tumor_bam,
                    bam_index = tumor_bai,
                    log_dir = log_dir,
                    task_name = ClipFiltering_task_name,
                    out_file = ClipFiltering_output_txt,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }
    }

    # String MergeVCFs_FilterMutectCalls_basic_output_vcf = "~{output_dir}/~{sample}.basic_filters_applied.vcf.gz"
    # String MergeVCFs_FilterMutectCalls_basic_output_idx = "~{MergeVCFs_FilterMutectCalls_basic_output_vcf}.tbi"
    # if (CheckProgress.task_completion[MergeVCFs_FilterMutectCalls_basic_task_name] == "false") {
    #     Array [String] MergeVCFs_FilterMutectCalls_basic_inputs = if (length(select_all(FilterMutectCalls_basic.output_vcf)) == length(intervals)) then select_all(FilterMutectCalls_basic.output_vcf) else FilterMutectCalls_basic_output_vcf
    #     call MergeVCFs as MergeVCFs_FilterMutectCalls_basic {
    #         input:
    #             input_vcfs = MergeVCFs_FilterMutectCalls_basic_inputs,
    #             log_dir = log_dir,
    #             task_name = MergeVCFs_FilterMutectCalls_basic_task_name,
    #             out_vcf = MergeVCFs_FilterMutectCalls_basic_output_vcf,
    #             out_idx = MergeVCFs_FilterMutectCalls_basic_output_idx,
    #             environment = environment,
    #             sa_key = sa_key,
    #             progress_dir = progress_dir,
    #             output_dir = output_dir
    #     }
    # }

    # String MergeVCFs_FilterMutectCalls_OB_output_vcf = "~{output_dir}/~{sample}.OB_filtering_applied.vcf.gz"
    # String MergeVCFs_FilterMutectCalls_OB_output_idx = "~{MergeVCFs_FilterMutectCalls_OB_output_vcf}.tbi"
    # if (CheckProgress.task_completion[MergeVCFs_FilterMutectCalls_OB_task_name] == "false") {
    #     Array [String] MergeVCFs_FilterMutectCalls_OB_inputs = if (length(select_all(FilterMutectCalls_OB.output_vcf)) == length(intervals)) then select_all(FilterMutectCalls_OB.output_vcf) else FilterMutectCalls_OB_output_vcf
    #     call MergeVCFs as MergeVCFs_FilterMutectCalls_OB {
    #         input:
    #             input_vcfs = MergeVCFs_FilterMutectCalls_OB_inputs,
    #             log_dir = log_dir,
    #             task_name = MergeVCFs_FilterMutectCalls_OB_task_name,
    #             out_vcf = MergeVCFs_FilterMutectCalls_OB_output_vcf,
    #             out_idx = MergeVCFs_FilterMutectCalls_OB_output_idx,
    #             environment = environment,
    #             sa_key = sa_key,
    #             progress_dir = progress_dir,
    #             output_dir = output_dir
    #     }
    # }

    # String MergeVCFs_FilterMutectCalls_CONT_output_vcf = "~{output_dir}/~{sample}.CONT_filters_applied.vcf.gz"
    # String MergeVCFs_FilterMutectCalls_CONT_output_idx = "~{MergeVCFs_FilterMutectCalls_CONT_output_vcf}.tbi"
    # if (CheckProgress.task_completion[MergeVCFs_FilterMutectCalls_CONT_task_name] == "false") {
    #     Array [String] MergeVCFs_FilterMutectCalls_CONT_inputs = if (length(select_all(FilterMutectCalls_CONT.output_vcf)) == length(intervals)) then select_all(FilterMutectCalls_CONT.output_vcf) else FilterMutectCalls_CONT_output_vcf
    #     call MergeVCFs as MergeVCFs_FilterMutectCalls_CONT {
    #         input:
    #             input_vcfs = MergeVCFs_FilterMutectCalls_CONT_inputs,
    #             log_dir = log_dir,
    #             task_name = MergeVCFs_FilterMutectCalls_CONT_task_name,
    #             out_vcf = MergeVCFs_FilterMutectCalls_CONT_output_vcf,
    #             out_idx = MergeVCFs_FilterMutectCalls_CONT_output_idx,
    #             environment = environment,
    #             sa_key = sa_key,
    #             progress_dir = progress_dir,
    #             output_dir = output_dir
    #     }
    # }

    String MergeVCFs_FilterMutectCalls_OB_CONT_output_vcf = "~{output_dir}/~{sample}.OB_CONT_filters_applied.vcf.gz"
    String MergeVCFs_FilterMutectCalls_OB_CONT_output_idx = "~{MergeVCFs_FilterMutectCalls_OB_CONT_output_vcf}.tbi"
    if (CheckProgress.task_completion[MergeVCFs_FilterMutectCalls_OB_CONT_task_name] == "false") {
        Array [String] MergeVCFs_FilterMutectCalls_OB_CONT_inputs = if (length(select_all(FilterMutectCalls_OB_CONT.output_vcf)) == length(intervals)) then select_all(FilterMutectCalls_OB_CONT.output_vcf) else FilterMutectCalls_OB_CONT_output_vcf
        call MergeVCFs as MergeVCFs_FilterMutectCalls_OB_CONT {
            input:
                input_vcfs = MergeVCFs_FilterMutectCalls_OB_CONT_inputs,
                log_dir = log_dir,
                task_name = MergeVCFs_FilterMutectCalls_OB_CONT_task_name,
                out_vcf = MergeVCFs_FilterMutectCalls_OB_CONT_output_vcf,
                out_idx = MergeVCFs_FilterMutectCalls_OB_CONT_output_idx,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String MergeVCFs_NormalizeVCF_output_vcf = "~{output_dir}/~{sample}.OB_CONT_filters_applied.normalized.vcf.gz"
    String MergeVCFs_NormalizeVCF_output_idx = "~{MergeVCFs_NormalizeVCF_output_vcf}.tbi"
    if (CheckProgress.task_completion[MergeVCFs_NormalizeVCF_task_name] == "false") {
        Array [String] MergeVCFs_NormalizeVCF_inputs = if (length(select_all(NormalizeVCF.output_vcf)) == length(intervals)) then select_all(NormalizeVCF.output_vcf) else NormalizeVCF_output_vcf
        call MergeVCFs as MergeVCFs_NormalizeVCF {
            input:
                input_vcfs = MergeVCFs_NormalizeVCF_inputs,
                log_dir = log_dir,
                task_name = MergeVCFs_NormalizeVCF_task_name,
                out_vcf = MergeVCFs_NormalizeVCF_output_vcf,
                out_idx = MergeVCFs_NormalizeVCF_output_idx,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String MergeAnnotation_output_snv_file = "~{output_dir}/~{sample}_annotated_filtered_snv.txt"
    String MergeAnnotation_output_indel_file = sub(MergeAnnotation_output_snv_file, "_snv.txt$", "_indel.txt")
    String MergeAnnotation_output_target_file = "~{output_dir}/~{sample}.tmp_filtered_variants.tsv"
    String FilterInfoFile = "~{output_dir}/~{sample}_filter_info.txt"
    String AnnotationFile = "~{output_dir}/~{sample}_annots.txt"
    if (CheckProgress.task_completion[MergeAnnotation_task_name] == "false") {
        Array[String] SNVFiles = if (length(select_all(ClipFiltering.output_txt)) == length(intervals)) then select_all(ClipFiltering.output_txt) else ClipFiltering_output_txt
        Array[String] IndelFiles = if (length(select_all(SSM_Filtering.output_indel_file)) == length(intervals)) then select_all(SSM_Filtering.output_indel_file) else SSM_Filtering_output_indel_file
        Array[String] FilterInfoFiles = if (length(select_all(SNVPostFiltering.output_filter_info)) == length(intervals)) then select_all(SNVPostFiltering.output_filter_info) else SSM_Filtering_output_filter_info
        Array[String] AnnotationFiles = if (length(select_all(SNVPostFiltering.output_annot_file)) == length(intervals)) then select_all(SNVPostFiltering.output_annot_file) else SSM_Filtering_output_annot_file
        call MergeAnnotation {
            input:
                snvs = SNVFiles,
                indels = IndelFiles,
                filter_info_files = FilterInfoFiles,
                annot_files = AnnotationFiles,
                log_dir = log_dir,
                task_name = MergeAnnotation_task_name,
                out_snv_file = MergeAnnotation_output_snv_file,
                out_indel_file = MergeAnnotation_output_indel_file,
                out_target_file = MergeAnnotation_output_target_file,
                filter_info = FilterInfoFile,
                annot_file = AnnotationFile,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String FilterVCF_output_vcf = "~{output_dir}/~{sample}.filtered.vcf.gz"
    String FilterVCF_output_idx = "~{FilterVCF_output_vcf}.tbi"
    if (CheckProgress.task_completion[FilterVCF_task_name] == "false") {
        call FilterVCF {
            input:
                merged_vcf = select_first([MergeVCFs_NormalizeVCF.output_vcf, MergeVCFs_NormalizeVCF_output_vcf]),
                reference = reference,
                reference_index = reference_index,
                reference_dict = reference_dict,
                log_dir = log_dir,
                task_name = FilterVCF_task_name,
                targets = select_first([MergeAnnotation.output_target_file, MergeAnnotation_output_target_file]),
                out_vcf = FilterVCF_output_vcf,
                out_idx = FilterVCF_output_idx,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String AnnotVCF_output_vcf = "~{output_dir}/~{sample}.annotated_filtered.vcf.gz"
    String AnnotVCF_output_idx = "~{AnnotVCF_output_vcf}.tbi"
    if (CheckProgress.task_completion[AnnotVCF_task_name] == "false") {
        call AnnotVCF {
            input:
                merged_vcf = select_first([MergeVCFs_NormalizeVCF.output_vcf, MergeVCFs_NormalizeVCF_output_vcf]),
                reference = reference,
                reference_index = reference_index,
                reference_dict = reference_dict,
                annot_file = select_first([MergeAnnotation.output_annot_file, AnnotationFile]),
                log_dir = log_dir,
                task_name = AnnotVCF_task_name,
                out_vcf = AnnotVCF_output_vcf,
                out_idx = AnnotVCF_output_idx,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String FilterGraph_output_file = "~{output_dir}/~{sample}.filtering.png"
    if (CheckProgress.task_completion[FilterGraph_task_name] == "false") {
        call FilterGraph {
            input:
                sample = sample,
                filter_info = select_first([MergeAnnotation.output_filter_info, FilterInfoFile]),
                log_dir = log_dir,
                task_name = FilterGraph_task_name,
                out_file = FilterGraph_output_file,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    if (CheckProgress.task_completion[Cleanup_task_name] == "false") {
        Array[String] Mutect2_vcfs = if (length(select_all(Mutect2.output_vcf)) == length(intervals)) then select_all(Mutect2.output_vcf) else Mutect2_output_vcf
        Array[String] Mutect2_idxs = if (length(select_all(Mutect2.output_tbi)) == length(intervals)) then select_all(Mutect2.output_tbi) else Mutect2_output_tbi
        Array[String] Mutect2_stats = if (length(select_all(Mutect2.output_stats)) == length(intervals)) then select_all(Mutect2.output_stats) else Mutect2_output_stats
        Array[String] Mutect2_f1r2s = if (length(select_all(Mutect2.output_f1r2)) == length(intervals)) then select_all(Mutect2.output_f1r2) else Mutect2_output_f1r2
        Array[String] Normal_pileups = if (length(select_all(GetNormalPileupSummaries.output_table)) == length(pileup_intervals)) then select_all(GetNormalPileupSummaries.output_table) else GetNormalPileupSummaries_output_table
        Array[String] Tumor_pileups = if (length(select_all(GetTumorPileupSummaries.output_table)) == length(pileup_intervals)) then select_all(GetTumorPileupSummaries.output_table) else GetTumorPileupSummaries_output_table
        # Array[String] basic_vcfs = if (length(select_all(FilterMutectCalls_basic.output_vcf)) == length(intervals)) then select_all(FilterMutectCalls_basic.output_vcf) else FilterMutectCalls_basic_output_vcf
        # Array[String] basic_idxs = if (length(select_all(FilterMutectCalls_basic.output_idx)) == length(intervals)) then select_all(FilterMutectCalls_basic.output_idx) else FilterMutectCalls_basic_output_idx
        # Array[String] basic_stats = if (length(select_all(FilterMutectCalls_basic.output_filtering_stats)) == length(intervals)) then select_all(FilterMutectCalls_basic.output_filtering_stats) else FilterMutectCalls_basic_output_filtering_stats
        # Array[String] OB_vcfs = if (length(select_all(FilterMutectCalls_OB.output_vcf)) == length(intervals)) then select_all(FilterMutectCalls_OB.output_vcf) else FilterMutectCalls_OB_output_vcf
        # Array[String] OB_idxs = if (length(select_all(FilterMutectCalls_OB.output_idx)) == length(intervals)) then select_all(FilterMutectCalls_OB.output_idx) else FilterMutectCalls_OB_output_idx
        # Array[String] OB_stats = if (length(select_all(FilterMutectCalls_OB.output_filtering_stats)) == length(intervals)) then select_all(FilterMutectCalls_OB.output_filtering_stats) else FilterMutectCalls_OB_output_filtering_stats
        # Array[String] CONT_vcfs = if (length(select_all(FilterMutectCalls_CONT.output_vcf)) == length(intervals)) then select_all(FilterMutectCalls_CONT.output_vcf) else FilterMutectCalls_CONT_output_vcf
        # Array[String] CONT_idxs = if (length(select_all(FilterMutectCalls_CONT.output_idx)) == length(intervals)) then select_all(FilterMutectCalls_CONT.output_idx) else FilterMutectCalls_CONT_output_idx
        # Array[String] CONT_stats = if (length(select_all(FilterMutectCalls_CONT.output_filtering_stats)) == length(intervals)) then select_all(FilterMutectCalls_CONT.output_filtering_stats) else FilterMutectCalls_CONT_output_filtering_stats
        Array[String] OB_CONT_vcfs = if (length(select_all(FilterMutectCalls_OB_CONT.output_vcf)) == length(intervals)) then select_all(FilterMutectCalls_OB_CONT.output_vcf) else FilterMutectCalls_OB_CONT_output_vcf
        Array[String] OB_CONT_idxs = if (length(select_all(FilterMutectCalls_OB_CONT.output_idx)) == length(intervals)) then select_all(FilterMutectCalls_OB_CONT.output_idx) else FilterMutectCalls_OB_CONT_output_idx
        Array[String] OB_CONT_stats = if (length(select_all(FilterMutectCalls_OB_CONT.output_filtering_stats)) == length(intervals)) then select_all(FilterMutectCalls_OB_CONT.output_filtering_stats) else FilterMutectCalls_OB_CONT_output_filtering_stats
        Array[String] Norm_vcfs = if (length(select_all(NormalizeVCF.output_vcf)) == length(intervals)) then select_all(NormalizeVCF.output_vcf) else NormalizeVCF_output_vcf
        Array[String] Norm_idxs = if (length(select_all(NormalizeVCF.output_idx)) == length(intervals)) then select_all(NormalizeVCF.output_idx) else NormalizeVCF_output_idx
        Array[String] VCF2Annovar_files = if (length(select_all(VCF2Annovar.output_file)) == length(intervals)) then select_all(VCF2Annovar.output_file) else VCF2Annovar_output_file
        Array[String] Annovar_files = if (length(select_all(Annovar.output_file)) == length(intervals)) then select_all(Annovar.output_file) else Annovar_output_file
        Array[String] SSM_Filtering_snvs_indels = if (length(select_all(SSM_Filtering.output_snvindel_file)) == length(intervals)) then select_all(SSM_Filtering.output_snvindel_file) else SSM_Filtering_output_snvindel_file
        Array[String] SSM_Filtering_snvs = if (length(select_all(SSM_Filtering.output_snv_file)) == length(intervals)) then select_all(SSM_Filtering.output_snv_file) else SSM_Filtering_output_snv_file
        Array[String] SSM_Filtering_indels = if (length(select_all(SSM_Filtering.output_indel_file)) == length(intervals)) then select_all(SSM_Filtering.output_indel_file) else SSM_Filtering_output_indel_file
        Array[String] SSM_Filtering_pypes = if (length(select_all(SSM_Filtering.output_pype_file)) == length(intervals)) then select_all(SSM_Filtering.output_pype_file) else SSM_Filtering_output_pype_file
        Array[String] SNVFiltering_files = if (length(select_all(SNVFiltering.output_tab)) == length(intervals)) then select_all(SNVFiltering.output_tab) else SNVFiltering_output_tab
        Array[String] SNVPostFiltering_files = if (length(select_all(SNVPostFiltering.output_txt)) == length(intervals)) then select_all(SNVPostFiltering.output_txt) else SNVPostFiltering_output_txt
        Array[String] ClipFiltering_files = if (length(select_all(ClipFiltering.output_txt)) == length(intervals)) then select_all(ClipFiltering.output_txt) else ClipFiltering_output_txt
        

        Array[Array[String]?] Cleanup_inputs = [Mutect2_vcfs, Mutect2_idxs, Mutect2_stats, Mutect2_f1r2s, Normal_pileups, Tumor_pileups, OB_CONT_vcfs, OB_CONT_idxs, OB_CONT_stats, Norm_vcfs, Norm_idxs, VCF2Annovar_files, Annovar_files, SSM_Filtering_snvs_indels, SSM_Filtering_snvs, SSM_Filtering_indels, SSM_Filtering_pypes, SNVFiltering_files, SNVPostFiltering_files, ClipFiltering_files, FilterInfoFiles, AnnotationFiles, select_all([GatherNormalPileupSummaries.output_table, GatherTumorPileupSummaries.output_table, MergeVCFs_FilterMutectCalls_OB_CONT.output_vcf, MergeVCFs_FilterMutectCalls_OB_CONT.output_idx])]
        call Cleanup {
            input:
                filtered_vcf = select_first([FilterVCF.output_vcf, FilterVCF_output_vcf]),
                annotated_vcf = select_first([AnnotVCF.output_vcf, AnnotVCF_output_vcf]),
                inputs = flatten(select_all(Cleanup_inputs)),
                log_dir = log_dir,
                task_name = Cleanup_task_name,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }


    output {
        String annotated_filtered_snv_file = select_first([MergeAnnotation.output_snv_file, MergeAnnotation_output_snv_file])
        String annotated_filtered_indel_file = select_first([MergeAnnotation.output_indel_file, MergeAnnotation_output_indel_file])
        String unfiltered_vcf = select_first([MergeVCFs_Mutect2.output_vcf, MergeVCFs_Mutect2_output_vcf])
        String filtered_vcf = select_first([Cleanup.output_vcf, FilterVCF_output_vcf])
        String filtered_vcf_idx = select_first([FilterVCF.output_idx, FilterVCF_output_idx])
        String annotated_filtered_vcf = select_first([Cleanup.output_annotated_vcf, AnnotVCF_output_vcf])
        String annotated_filtered_idx = select_first([AnnotVCF.output_idx, AnnotVCF_output_idx])
        String filtering_graph = select_first([FilterGraph.output_file, FilterGraph_output_file])

        # String basic_filtered_vcf = select_first([MergeVCFs_FilterMutectCalls_basic.output_vcf, MergeVCFs_FilterMutectCalls_basic_output_vcf])
        # String OB_filtered_vcf = select_first([MergeVCFs_FilterMutectCalls_OB.output_vcf, MergeVCFs_FilterMutectCalls_OB_output_vcf])
        # String CONT_filtered_vcf = select_first([MergeVCFs_FilterMutectCalls_CONT.output_vcf, MergeVCFs_FilterMutectCalls_CONT_output_vcf])
        String OB_CONT_filtered_vcf = select_first([MergeVCFs_FilterMutectCalls_OB_CONT.output_vcf, MergeVCFs_FilterMutectCalls_OB_CONT_output_vcf])
    }
}

task Mutect2 {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        File reference
        File reference_index
        File reference_dict
        String tumor_name
        String normal_name
        String interval
        String log_dir
        String task_name
        String out_vcf
        String out_tbi
        String out_stats
        String out_f1r2
        Int preemptible

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String sa_key_base = basename(sa_key)

    String out_vcf_base = basename(out_vcf)
    String out_tbi_base = basename(out_tbi)
    String out_stats_base = basename(out_stats)
    String out_f1r2_base = basename(out_f1r2)

    Int disk_size = ceil((size(reference, "GB") + size(tumor_bam, "GB") + size(normal_bam, "GB")) * 2 + 50)

    command <<<
        set -eo pipefail

        echo "In Mutect2"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        my_out_vcf=~{out_vcf_base}
        my_out_tbif=~{out_vcf_base}
        my_out_stats=~{out_vcf_base}
        my_out_f1r2=~{out_vcf_base}

        env

        echo "Going into our logic statments"
        if [ "~{environment}" = "LOCAL" ]; then
            echo "We are on LOCAL"
            # Check the size of the scratch space, df -k gives size in 1-K blocks
            df -k
            df -k /scratch
            MYSIZE=$(df -k /scratch | grep scratch | awk '{print $2}')
            echo "$MYSIZE"
            if [ "$MYSIZE" -gt 1073741824 ]; then
                # We are on one of the larger nodes
                # Check the free space available
                FREE=$(df -k /scratch | grep scratch | awk '{print $4}')
                echo "$FREE"
                if [ "$FREE" -gt 4294967296 ]; then
                    echo "Skipping moving bams for this task, it seems to slow down this task."
                    #echo "Copying bams to scratch space..."
                    #cp ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai} /scratch/
                    #rm ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai}
                    #my_tumor_bam=/scratch/~{tumor_bam_base}
                    #my_tumor_bai=/scratch/~{tumor_bai_base}
                    #my_normal_bam=/scratch/~{normal_bam_base}
                    #my_normal_bai=/scratch/~{normal_bai_base}
                fi
            fi

            echo "Copying reference to scratch space..."
            cp ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key} /scratch/
            rm ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key}
            my_reference=/scratch/~{reference_base}
            my_reference_index=/scratch/~{reference_index_base}
            my_reference_dict=/scratch/~{reference_dict_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_vcf=/scratch/~{out_vcf_base}
            my_out_tbif=/scratch/~{out_tbi_base}
            my_out_stats=/scratch/~{out_stats_base}
            my_out_f1r2=/scratch/~{out_f1r2_base}
        fi


        time \
        gatk --java-options "-Xmx24g" \
        Mutect2 \
        -R $my_reference \
        -I $my_tumor_bam \
        -tumor ~{tumor_name} \
        -I $my_normal_bam \
        -normal ~{normal_name} \
        -L ~{interval} \
        -O $my_out_vcf \
        --max-reads-per-alignment-start 0 \
        --max-mnp-distance 0 \
        --f1r2-tar-gz $my_out_f1r2

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_vcf ~{output_dir}
                cp $my_out_tbif ~{output_dir}
                cp $my_out_stats ~{output_dir}
                cp $my_out_f1r2 ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_vcf_base} ~{output_dir}
                gsutil cp ~{out_tbi_base} ~{output_dir}
                gsutil cp ~{out_stats_base} ~{output_dir}
                gsutil cp ~{out_f1r2_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_vcf = "~{out_vcf}"
        String output_tbi = "~{out_tbi}"
        String output_stats = "~{out_stats}"
        String output_f1r2 = "~{out_f1r2}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 4
        memory: "48 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible
        walltime: "96:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task VCF2Annovar {
    input {
        File vcf
        File idx
        String tumor_name
        String normal_name
        String log_dir
        String task_name
        String out_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String vcf_base = basename(vcf)
    String idx_base = basename(idx)
    String sa_key_base = basename(sa_key)

    String out_file_base = basename(out_file)
    Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

    command {
        set -eo pipefail

        echo "In VCF2Annovar"
        env

        my_vcf=~{vcf}
        my_idx=~{idx}
        my_sa_key=~{sa_key}

        my_out_file=~{out_file_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{vcf} ~{idx} ~{sa_key} /scratch/
            rm ~{vcf} ~{idx} ~{sa_key}
            my_vcf=/scratch/~{vcf_base}
            my_idx=/scratch/~{idx_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_file=/scratch/~{out_file_base}
        fi

        time \
            mutect2annovar.pl \
            --vcf $my_vcf \
            --filter false \
            --header false \
            --output $my_out_file \
            --tumour ~{tumor_name} \
            --normal ~{normal_name}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_file ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_file_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_file = "~{out_file}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/vcftools-perl:5.22.1"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task Annovar {
    input {
        File input_file
        File humandb_tar
        String bedfile
        String build_version
        String log_dir
        String task_name
        String out_prefix
        String out_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String humandb_tar_base = basename(humandb_tar)
    String sa_key_base = basename(sa_key)

    String out_prefix_base = basename(out_prefix)
    String out_file_base = basename(out_file)
    Int disk_size = ceil(size(input_file, "GB") * 2 + 100)

    command {
        set -eo pipefail

        echo "In Mutect2"
        env

        my_humandb_tar=~{humandb_tar}
        HUMANDB=$(pwd)/$(basename ~{humandb_tar} ".tar.gz")
        echo "$HUMANDB"
        my_sa_key=~{sa_key}

        my_out_prefix=~{out_prefix_base}
        my_out_file=~{out_file_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{humandb_tar} ~{sa_key} /scratch/
            rm ~{humandb_tar} ~{sa_key}
            my_humandb_tar=/scratch/~{humandb_tar_base}
            HUMANDB=/scratch/$(basename ~{humandb_tar} ".tar.gz")
            echo "$HUMANDB"
            my_sa_key=/scratch/~{sa_key_base}

            my_out_prefix=/scratch/~{out_prefix_base}
            my_out_file=/scratch/~{out_file_base}

            cd /scratch

            #ls "$HUMANDB"
        fi

        tar -zxvf $my_humandb_tar

        if [ "~{environment}" = "CLOUD" ]; then
            mv ~{input_file} .
            INFILE_LOC="$(pwd)/~{basename(input_file)}"
        else
            INFILE_LOC="~{input_file}"
        fi

        if [ ! -e "$HUMANDB/~{basename(bedfile)}" ]; then
            cp ~{bedfile} $HUMANDB
        fi

        time \
        table_annovar.pl \
            "$INFILE_LOC" \
            $HUMANDB \
            --protocol refGene,ensGene,snp132,1000g2012feb_all,esp6500si_all,cg69,cosmic70,clinvar_20150330,exac03,bed \
            --operation g,g,f,f,f,f,f,f,f,r \
            --buildver ~{build_version} \
            --remove \
            --otherinfo \
            --bedfile ~{basename(bedfile)} \
            --outfile $my_out_prefix

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_file ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_file_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_file = "~{out_file}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/annovar:2013.08.23"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task SSM_Filtering {
    input {
        File annotations
        File vcf
        File idx

        String tumor_name
        String normal_name
        String source
        File cosmic_rda
        String log_dir
        String task_name

        String out_snvindel_file
        String out_snv_file
        String out_pype_file
        String out_indel_file
        String filter_info
        String annot_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String annotations_base = basename(annotations)
    String vcf_base = basename(vcf)
    String idx_base = basename(idx)
    String cosmic_rda_base = basename(cosmic_rda)
    String sa_key_base = basename(sa_key)

    String out_snvindel_file_base = basename(out_snvindel_file)
    String out_snv_file_base = basename(out_snv_file)
    String out_pype_file_base = basename(out_pype_file)
    String out_indel_file_base = basename(out_indel_file)
    String filter_info_base = basename(filter_info)
    String annot_file_base = basename(annot_file)

    command {
        set -eo pipefail

        echo "In SSM_Filtering"
        env

        my_annotations=~{annotations}
        my_vcf=~{vcf}
        my_idx=~{idx}
        my_cosmic_rda=~{cosmic_rda}
        my_sa_key=~{sa_key}

        my_out_snvindel_file=~{out_snvindel_file_base}
        my_out_snv_file=~{out_snv_file_base}
        my_out_pype_file=~{out_pype_file_base}
        my_out_indel_file=~{out_indel_file_base}
        my_filter_info=~{filter_info_base}
        my_annot_file=~{annot_file_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{annotations} ~{vcf} ~{idx} ~{cosmic_rda} ~{sa_key} /scratch/
            rm ~{annotations} ~{vcf} ~{idx} ~{cosmic_rda} ~{sa_key}
            my_annotations=/scratch/~{annotations_base}
            my_vcf=/scratch/~{vcf_base}
            my_idx=/scratch/~{idx_base}
            my_cosmic_rda=/scratch/~{cosmic_rda_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_snvindel_file=/scratch/~{out_snvindel_file_base}
            my_out_snv_file=/scratch/~{out_snv_file_base}
            my_out_pype_file=/scratch/~{out_pype_file_base}
            my_out_indel_file=/scratch/~{out_indel_file_base}
            my_filter_info=/scratch/~{filter_info_base}
            my_annot_file=/scratch/~{annot_file_base}
        fi

        time \
        Rscript /opt/scripts/cromwell_run_ssm_standard_filters.R \
            $my_cosmic_rda \
            ~{source} \
            $my_out_snvindel_file \
            $my_out_snv_file \
            $my_out_pype_file \
            $my_out_indel_file \
            $my_filter_info \
            $my_annot_file \
            $my_annotations \
            $my_vcf \
            ~{tumor_name} \
            ~{normal_name}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_snvindel_file ~{output_dir}
                cp $my_out_snv_file ~{output_dir}
                cp $my_out_pype_file ~{output_dir}
                cp $my_out_indel_file ~{output_dir}
                cp $my_filter_info ~{output_dir}
                cp $my_annot_file ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_snvindel_file_base} ~{output_dir}
                gsutil cp ~{out_snv_file_base} ~{output_dir}
                gsutil cp ~{out_pype_file_base} ~{output_dir}
                gsutil cp ~{out_indel_file_base} ~{output_dir}
                gsutil cp ~{filter_info_base} ~{output_dir}
                gsutil cp ~{annot_file_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_snvindel_file = "~{out_snvindel_file}"
        String output_snv_file = "~{out_snv_file}"
        String output_pype_file = "~{out_pype_file}"
        String output_indel_file = "~{out_indel_file}"
        String output_filter_info = "~{filter_info}"
        String output_annot_file = "~{annot_file}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk 50 HDD"
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task SNVFiltering {
    input {
        File input_file
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        File gatk_path
        File reference
        File reference_index
        File reference_dict
        File centromeres
        File dustmaker
        String log_dir
        String task_name
        String out_tab

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String input_file_base = basename(input_file)
    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String gatk_path_base = basename(gatk_path)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String centromeres_base = basename(centromeres)
    String dustmaker_base = basename(dustmaker)
    String sa_key_base = basename(sa_key)

    String out_tab_base = basename(out_tab)

    Int disk_size = ceil((size(input_file, "GB") + size(tumor_bam, "GB") + size(normal_bam, "GB")) * 2 + 20)

    command <<<
        set -eo pipefail

        echo "In SNVFiltering"
        env

        my_input_file=~{input_file}
        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_gatk_path=~{gatk_path}
        my_centromeres=~{centromeres}
        my_dustmaker=~{dustmaker}
        my_sa_key=~{sa_key}

        my_out_tab=~{out_tab_base}

        if [ "~{environment}" = "LOCAL" ]; then
            echo "We are on LOCAL"
            # Check the size of the scratch space, df -k gives size in 1-K blocks
            df -k
            df -k /scratch
            MYSIZE=$(df -k /scratch | grep scratch | awk '{print $2}')
            echo "$MYSIZE"
            if [ "$MYSIZE" -gt 1073741824 ]; then
                # We are on one of the larger nodes
                # Check the free space available
                FREE=$(df -k /scratch | grep scratch | awk '{print $4}')
                echo "$FREE"
                if [ "$FREE" -gt 4294967296 ]; then
                    echo "Skipping moving bams for this task, it seems to slow down this task."
                    #cp ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai} /scratch/
                    #rm ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai}
                    #my_tumor_bam=/scratch/~{tumor_bam_base}
                    #my_tumor_bai=/scratch/~{tumor_bai_base}
                    #my_normal_bam=/scratch/~{normal_bam_base}
                    #my_normal_bai=/scratch/~{normal_bai_base}
                fi
            fi

            cp ~{input_file} ~{gatk_path} /scratch/
            rm ~{input_file} ~{gatk_path}
            my_input_file=/scratch/~{input_file_base}
            my_gatk_path=/scratch/~{gatk_path_base}

            cp ~{reference} ~{reference_dict} ~{centromeres} ~{dustmaker} ~{sa_key} /scratch/
            rm ~{reference} ~{reference_dict} ~{centromeres} ~{dustmaker} ~{sa_key}
            my_reference=/scratch/~{reference_base}
            my_reference_dict=/scratch/~{reference_dict_base}
            my_centromeres=/scratch/~{centromeres_base}
            my_dustmaker=/scratch/~{dustmaker_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_tab=/scratch/~{out_tab_base}
        fi

        time \
        python3 -B /opt/filterpipeline/runFilters.py \
            --tumor_tab $my_input_file \
            --tumor_bam $my_tumor_bam \
            --normal_tab none \
            --normal_bam $my_normal_bam \
            --config /opt/filterpipeline/config.yaml \
            --reference $my_reference \
            --centromeres $my_centromeres \
            --dustmaker $my_dustmaker \
            --cFilter \
            --output_tab $my_out_tab \
            --gatk_path $my_gatk_path

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_tab ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_tab_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_tab = "~{out_tab}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-filtering:3.6"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task SNVPostFiltering {
    input {
        File cfilter_input
        File snv_input
        File pon
        Int pon_max
        Int min_dist_complex
        Int max_number_of_flagged
        Int normal_alt_count_max
        String log_dir
        String task_name
        String out_file
        File filter_info
        File annot_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String cfilter_input_base = basename(cfilter_input)
    String snv_input_base = basename(snv_input)
    String pon_base = basename(pon)
    String sa_key_base = basename(sa_key)

    String out_file_base = basename(out_file)
    String filter_info_base = if (environment == "LOCAL") then filter_info else basename(filter_info)
    String annot_file_base = if (environment == "LOCAL") then annot_file else basename(annot_file)

    Int disk_size = ceil((size(cfilter_input, "GB") + size(snv_input, "GB")) * 3 + 20)

    command {
        set -eo pipefail

        echo "In SNVPostFiltering"
        env

        my_cfilter_input=~{cfilter_input}
        my_snv_input=~{snv_input}
        my_pon=~{pon}
        my_sa_key=~{sa_key}

        my_out_file=~{out_file_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{cfilter_input} ~{snv_input} ~{pon} ~{sa_key} /scratch/
            rm ~{cfilter_input} ~{snv_input} ~{pon} ~{sa_key}
            my_cfilter_input=/scratch/~{cfilter_input_base}
            my_snv_input=/scratch/~{snv_input_base}
            my_pon=/scratch/~{pon_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_file=/scratch/~{out_file_base}
        fi

        cp ~{filter_info} .
        cp ~{annot_file} .

        time \
        Rscript /opt/scripts/cromwell_run_snv_postfiltering.R \
            --snv $my_snv_input \
            --cfilter $my_cfilter_input \
            --output_file $my_out_file \
            --pon $my_pon \
            --pon_max ~{pon_max} \
            --min_dist_complex ~{min_dist_complex} \
            --max_number_of_flagged ~{max_number_of_flagged} \
            --normal_alt_count_max ~{normal_alt_count_max} \
            --filter_info ~{filter_info_base} \
            --annot_file ~{annot_file_base}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp -uf $my_out_file ~{output_dir}
                # overwrite the input filter file with the updated one
                #cp -uf ~{filter_info_base} ~{output_dir}
                #cp -uf ~{annot_file_base} ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_file_base} ~{output_dir}
                gsutil cp ~{filter_info_base} ~{output_dir}
                gsutil cp ~{annot_file_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_txt = "~{out_file}"
        String output_filter_info = "~{filter_info}"
        String output_annot_file = "~{annot_file}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task ClipFiltering {
    input {
        File input_file
        File bam_file
        File bam_index
        String log_dir
        String task_name
        String out_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String input_file_base = basename(input_file)
    String bam_file_base = basename(bam_file)
    String bam_index_base = basename(bam_index)
    String sa_key_base = basename(sa_key)

    String out_file_base = basename(out_file)
    Int disk_size = ceil((size(input_file, "GB") + size(bam_file, "GB")) * 2 + 20)

    command <<<
        set -eo pipefail

        echo "In ClipFiltering"
        env

        my_input_file=~{input_file}
        my_bam_file=~{bam_file}
        my_bam_index=~{bam_index}
        my_sa_key=~{sa_key}

        my_out_file=~{out_file_base}

        if [ "~{environment}" = "LOCAL" ]; then
            echo "We are on LOCAL"
            # Check the size of the scratch space, df -k gives size in 1-K blocks
            df -k
            df -k /scratch
            MYSIZE=$(df -k /scratch | grep scratch | awk '{print $2}')
            echo "$MYSIZE"
            if [ "$MYSIZE" -gt 1073741824 ]; then
                # We are on one of the larger nodes
                # Check the free space available
                FREE=$(df -k /scratch | grep scratch | awk '{print $4}')
                echo "$FREE"
                if [ "$FREE" -gt 4294967296 ]; then
                    echo "Skipping moving bams for this task, it seems to slow down this task."
                    #cp ~{bam_file} ~{bam_index} /scratch/
                    #rm ~{bam_file} ~{bam_index}
                    #my_bam_file=/scratch/~{bam_file_base}
                    #my_bam_index=/scratch/~{bam_index_base}
                fi
            fi

            cp ~{input_file} ~{sa_key} /scratch/
            rm ~{input_file} ~{sa_key}
            my_input_file=/scratch/~{input_file_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_file=/scratch/~{out_file_base}
        fi

        time \
        python3 -B /opt/scripts/filter_clips.py \
            --tumor_bam $my_bam_file \
            --input_file $my_input_file \
            --output_file $my_out_file

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_file ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_file_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_txt = "~{out_file}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-filtering:3.6"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task MergeAnnotation {
    input {
        Array[File] snvs
        Array[File] indels
        Array[File] filter_info_files
        Array[File] annot_files
        String log_dir
        String task_name
        String out_snv_file
        String out_indel_file
        String out_target_file
        String filter_info
        String annot_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String sa_key_base = basename(sa_key)

    String out_snv_file_base = if (environment == "LOCAL") then out_snv_file else basename(out_snv_file)
    String out_indel_file_base = if (environment == "LOCAL") then out_indel_file else basename(out_indel_file)
    String out_target_file_base = if (environment == "LOCAL") then out_target_file else basename(out_target_file)
    String filter_info_base = if (environment == "LOCAL") then filter_info else basename(filter_info)
    String annot_file_base = if (environment == "LOCAL") then annot_file else basename(annot_file)

    Int disk_size = ceil(size(snvs[0], "GB") * length(snvs) + 50)

    command {
        set -eo pipefail

        echo "In MergeAnnotation"
        env

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            rm ~{sa_key}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        Rscript /opt/scripts/merge_annotations.R \
            ~{out_snv_file_base} \
            ~{out_indel_file_base} \
            ~{out_target_file} \
            ~{filter_info_base} \
            ~{annot_file_base} \
            ~{sep="," snvs} \
            ~{sep="," indels} \
            ~{sep="," filter_info_files} \
            ~{sep="," annot_files}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_snv_file_base} ~{output_dir}
                gsutil cp ~{out_indel_file_base} ~{output_dir}
                gsutil cp ~{out_target_file_base} ~{output_dir}
                gsutil cp ~{filter_info_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_snv_file = "~{out_snv_file}"
        String output_indel_file = "~{out_indel_file}"
        String output_filter_info = "~{filter_info}"
        String output_annot_file = "~{annot_file}"
        String output_target_file = "~{out_target_file}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task MergeVCFs {
    input {
        Array[File] input_vcfs
        String log_dir
        String task_name
        String out_vcf
        String out_idx

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String sa_key_base = basename(sa_key)

    String out_vcf_base = if (environment == "LOCAL") then out_vcf else basename(out_vcf)
    String out_idx_base = if (environment == "LOCAL") then out_idx else basename(out_idx)

    Int disk_size = ceil(size(input_vcfs[0], "GB")*length(input_vcfs) * 2 + 20)

    command {
        set -eo pipefail

        echo "In MergeVCFs"
        env

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            rm ~{sa_key}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        gatk --java-options "-Xmx24g" \
        MergeVcfs \
            -I ~{sep=' -I ' input_vcfs} \
            -O ~{out_vcf_base}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_vcf_base} ~{output_dir}
                gsutil cp ~{out_idx_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_vcf = "~{out_vcf}"
        String output_idx = "~{out_idx}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 2
        memory: "48 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "4:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task FilterVCF {
    input {
        File merged_vcf
        File reference
        File reference_index
        File reference_dict
        String log_dir
        String task_name
        File targets
        String out_vcf
        String out_idx

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String sa_key_base = basename(sa_key)

    String out_vcf_base = if (environment == "LOCAL") then out_vcf else basename(out_vcf)
    String out_idx_base = if (environment == "LOCAL") then out_idx else basename(out_idx)
    String targets_base = if (environment == "LOCAL") then targets else basename(targets)
    Int disk_size = ceil((size(targets_base, "GB") + size(merged_vcf, "GB")) * 2 + 50)

    command <<<
        set -eo pipefail

        echo "In FilterVCF"
        env

        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key} /scratch/
            rm ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key}
            my_reference=/scratch/~{reference_base}
            my_reference_index=/scratch/~{reference_index_base}
            my_reference_dict=/scratch/~{reference_dict_base}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        bcftools view \
            -T ~{targets_base} \
            -Oz \
            -o ~{out_vcf} \
            ~{merged_vcf} \
        && \
        tabix ~{out_vcf} \
        && \
        java -Xmx2g \
            -jar $GATK \
            -T ValidateVariants \
            -V ~{out_vcf_base} \
            -R $my_reference

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_vcf_base} ~{output_dir}
                gsutil cp ~{out_idx_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_vcf = "~{out_vcf}"
        String output_idx = "~{out_vcf}.idx"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk-r:3.5_3.4.0"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "4:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task AnnotVCF {
    input {
        File merged_vcf
        File reference
        File reference_index
        File reference_dict
        File annot_file
        String log_dir
        String task_name
        String out_vcf
        String out_idx

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String sa_key_base = basename(sa_key)

    String out_vcf_base = if (environment == "LOCAL") then out_vcf else basename(out_vcf)
    String out_idx_base = if (environment == "LOCAL") then out_idx else basename(out_idx)
    String annot_file_base = if (environment == "LOCAL") then annot_file else basename(annot_file)
    Int disk_size = ceil((size(annot_file_base, "GB") + size(merged_vcf, "GB")) * 2 + 50)

    command <<<
        set -eo pipefail

        echo "In AnnotVCF"
        env

        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key} /scratch/
            rm ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key}
            my_reference=/scratch/~{reference_base}
            my_reference_index=/scratch/~{reference_index_base}
            my_reference_dict=/scratch/~{reference_dict_base}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        FILTER_NAMES=( `tail -n +2 ~{annot_file} | cut -f4 | sort | uniq` )
        if [[ -f "~{out_vcf_base}.annots.hdr" ]]
        then
            rm "~{out_vcf_base}.annots.hdr"
        fi
        for f in ${FILTER_NAMES[@]}
        do
            printf "%sFILTER=<ID=%s,Description=\"Custom filter: %s\">\n" '##' "$f" "$f" >> "~{out_vcf_base}.annots.hdr"
        done \
        && \
        bgzip -f ~{annot_file_base} \
        && \
        tabix -S 1 -f -s1 -b2 -e3 ~{annot_file_base}.gz \
        && \
        bcftools annotate \
            -a ~{annot_file_base}.gz \
            -h "~{out_vcf_base}.annots.hdr" \
            -c CHROM,FROM,TO,FILTER \
            ~{merged_vcf} \
            -o ~{out_vcf_base} \
            -Oz \
        && \
        tabix ~{out_vcf_base} \
        && \
        rm "~{out_vcf_base}.annots.hdr" \
        && \
        java -Xmx2g \
            -jar $GATK \
            -T ValidateVariants \
            -V ~{out_vcf_base} \
            -R $my_reference

        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_vcf_base} ~{output_dir}
                gsutil cp ~{out_idx_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_vcf = "~{out_vcf}"
        String output_idx = "~{out_idx}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk-r:3.5_3.4.0"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "4:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task FilterGraph {
    input {
        String sample
        File filter_info
        String log_dir
        String task_name
        String out_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String filter_info_base = basename(filter_info)
    String sa_key_base = basename(sa_key)

    String out_file_base = basename(out_file)

    Int disk_size = ceil(size(filter_info, "GB") + 20)

    command {
        set -eo pipefail

        echo "In FilterGraph"
        env

        my_filter_info=~{filter_info}
        my_sa_key=~{sa_key}

        my_out_file=~{out_file_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{filter_info} ~{sa_key} /scratch/
            rm ~{filter_info} ~{sa_key}
            my_filter_info=/scratch/~{filter_info_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_file=/scratch/~{out_file_base}
        fi

        time \
        Rscript /opt/scripts/ssm_filtering_graph.R \
            ~{sample} \
            $my_filter_info \
            $my_out_file

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_file ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_file_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_file = "~{out_file}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "0:30:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task MergeStats {
    input {
        Array[File] input_stats
        String log_dir
        String task_name
        String out_stats

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String sa_key_base = basename(sa_key)

    String out_stats_base = if (environment == "LOCAL") then out_stats else basename(out_stats)

    Int disk_size = ceil(size(input_stats[0], "GB") * length(input_stats) * 4 + 20)

    command {
        set -eo pipefail

        echo "In MergeStats"
        env

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            rm ~{sa_key}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        gatk --java-options "-Xmx24g" \
        MergeMutectStats \
        --stats ~{sep=' --stats ' input_stats} \
        -O ~{out_stats_base}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_stats_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_stats = "~{out_stats}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 1
        memory: "48 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "4:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task GetPileupSummaries {
    input {
        File input_bam
        File input_bai
        File reference
        File reference_index
        File reference_dict
        File common_biallelic_germline_variants
        File common_biallelic_germline_variants_index
        String interval
        String log_dir
        String task_name
        String out_table

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String input_bam_base = basename(input_bam)
    String input_bai_base = basename(input_bai)
    String common_biallelic_germline_variants_base = basename(common_biallelic_germline_variants)
    String common_biallelic_germline_variants_index_base = basename(common_biallelic_germline_variants_index)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String sa_key_base = basename(sa_key)

    String out_table_base = basename(out_table)

    Int disk_size = ceil((size(input_bam, "GB") + size(reference, "GB")) * 1.2 + 20)

    command <<<
        set -eo pipefail

        echo "In GetPileupSummaries"
        env

        my_input_bam=~{input_bam}
        my_input_bai=~{input_bai}
        my_common_biallelic_germline_variants=~{common_biallelic_germline_variants}
        my_common_biallelic_germline_variants_index=~{common_biallelic_germline_variants_index}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        my_out_table=~{out_table_base}

        if [ "~{environment}" = "LOCAL" ]; then
            echo "We are on LOCAL"
            # Check the size of the scratch space, df -k gives size in 1-K blocks
            df -k
            df -k /scratch
            MYSIZE=$(df -k /scratch | grep scratch | awk '{print $2}')
            echo "$MYSIZE"
            if [ "$MYSIZE" -gt 1073741824 ]; then
                # We are on one of the larger nodes
                # Check the free space available
                FREE=$(df -k /scratch | grep scratch | awk '{print $4}')
                echo "$FREE"
                if [ "$FREE" -gt 4294967296 ]; then
                    echo "Skipping moving bams for this task, it seems to slow down this task."
                    #echo "Copying bams to scratch space..."
                    #cp ~{input_bam} ~{input_bai} /scratch/
                    #rm ~{input_bam} ~{input_bai}
                    #my_input_bam=/scratch/~{input_bam_base}
                    #my_input_bai=/scratch/~{input_bai_base}
                fi
            fi

            cp ~{common_biallelic_germline_variants} ~{common_biallelic_germline_variants_index} ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key} /scratch/
            #rm ~{common_biallelic_germline_variants} ~{common_biallelic_germline_variants_index} ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key}
            my_common_biallelic_germline_variants=/scratch/~{common_biallelic_germline_variants_base}
            my_common_biallelic_germline_variants_index=/scratch/~{common_biallelic_germline_variants_index_base}
            my_reference=/scratch/~{reference_base}
            my_reference_index=/scratch/~{reference_index_base}
            my_reference_dict=/scratch/~{reference_dict_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_table=/scratch/~{out_table_base}

            ls -l /scratch
            biallelic_dir=$(dirname $my_common_biallelic_germline_variants)
        elif [ "~{environment}" = "CLOUD" ]; then
            biallelic_dir=$(dirname $my_common_biallelic_germline_variants)
            mv -vn $my_common_biallelic_germline_variants_index $biallelic_dir
        fi

        #biallelic_dir=$(dirname $my_common_biallelic_germline_variants)
        #if [ ! -f biallelic_dir/$my_common_biallelic_germline_variants_index ]; then
        #    echo $biallelic_dir
        #    echo $my_common_biallelic_germline_variants
        #    echo $my_common_biallelic_germline_variants_index
        #    ls -l $biallelic_dir
        #    mv -vn $my_common_biallelic_germline_variants_index $biallelic_dir
        #fi

        time \
        gatk --java-options "-Xmx24g" \
        GetPileupSummaries \
        -R $my_reference \
        -I $my_input_bam \
        --interval-set-rule INTERSECTION \
        -L ~{interval} \
        -V $my_common_biallelic_germline_variants \
        -L $my_common_biallelic_germline_variants \
        -O $my_out_table

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_table ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_table_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_table = "~{out_table}"
    }

    # n1-highmem-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 4
        memory: "52 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptibe: 2
        walltime: "96:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task GatherPileupSummaries {
    input {
        Array[File] input_tables
        File reference_dict
        String log_dir
        String task_name
        String out_table

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String reference_dict_base = basename(reference_dict)
    String sa_key_base = basename(sa_key)

    String out_table_base = basename(out_table)

    command {
        set -eo pipefail

        echo "In GatherPileupSummaries"
        env

        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        my_out_table=~{out_table_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{reference_dict} ~{sa_key} /scratch/
            rm ~{reference_dict} ~{sa_key}
            my_reference_dict=/scratch/~{reference_dict_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_table=/scratch/~{out_table_base}
        fi

        time \
        gatk --java-options "-Xmx6g" \
        GatherPileupSummaries \
        --sequence-dictionary $my_reference_dict \
        -I ~{sep=' -I ' input_tables} \
        -O $my_out_table

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_table ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_table_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_table = "~{out_table}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 2
        memory: "12 GB"
        disks: "local-disk 20 HDD"
        preemptible: 2
        walltime: "24:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task CalculateContamination {
    input {
        File normal_pileups
        File tumor_pileups
        String log_dir
        String task_name
        String out_contamination
        String out_segmentation

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String normal_pileups_base = basename(normal_pileups)
    String tumor_pileups_base = basename(tumor_pileups)
    String sa_key_base = basename(sa_key)

    String out_contamination_base = basename(out_contamination)
    String out_segmentation_base = basename(out_segmentation)

    command {
        set -eo pipefail

        echo "In CalculateContamination"
        env

        my_normal_pileups=~{normal_pileups}
        my_tumor_pileups=~{tumor_pileups}
        my_sa_key=~{sa_key}

        my_out_contamination=~{out_contamination_base}
        my_out_segmentation=~{out_segmentation_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{normal_pileups} ~{tumor_pileups} ~{sa_key} /scratch/
            rm ~{normal_pileups} ~{tumor_pileups} ~{sa_key}
            my_normal_pileups=/scratch/~{normal_pileups_base}
            my_tumor_pileups=/scratch/~{tumor_pileups_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_contamination=/scratch/~{out_contamination_base}
            my_out_segmentation=/scratch/~{out_segmentation_base}
        fi

        time \
        gatk --java-options "-Xmx12g" \
        CalculateContamination \
        -I $my_tumor_pileups \
        -matched $my_normal_pileups \
        -O $my_out_contamination \
        --tumor-segmentation $my_out_segmentation

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_contamination ~{output_dir}
                cp $my_out_segmentation ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_contamination_base} ~{output_dir}
                gsutil cp ~{out_segmentation_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_contamination = "~{out_contamination}"
        String output_segmentation = "~{out_segmentation}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 4
        memory: "24 GB"
        disks: "local-disk 20 HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task LearnReadOrientationModel {
    input {
        Array[File] input_f1r2s
        String log_dir
        String task_name
        String out_artifact_priors

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String sa_key_base = basename(sa_key)

    String out_artifact_priors_base = if (environment == "LOCAL") then out_artifact_priors else basename(out_artifact_priors)

    Int disk_size = ceil(size(input_f1r2s[0], "GB") * length(input_f1r2s) * 1.5 + 20)

    command {
        set -eo pipefail

        echo "In LearnReadOrientationModel"
        env

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            rm ~{sa_key}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        gatk --java-options "-Xmx24g" \
        LearnReadOrientationModel \
        -I ~{sep=" -I " input_f1r2s} \
        -O ~{out_artifact_priors_base}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_artifact_priors_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_artifact_priors = "~{out_artifact_priors}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 4
        memory: "48 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task FilterMutectCalls {
    input {
        File reference
        File reference_index
        File reference_dict
        File input_vcf
        File input_vcf_index
        File mutect_stats
        Boolean filter_ob
        Boolean filter_cont
        File? artifact_priors
        File? contamination_table
        File? tumor_segmentation_table
        String interval
        String log_dir
        String task_name
        String out_vcf
        String out_idx
        String out_stats

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String input_vcf_base = basename(input_vcf)
    String input_vcf_index_base = basename(input_vcf_index)
    String mutect_stats_base = basename(mutect_stats)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String sa_key_base = basename(sa_key)

    String out_vcf_base = basename(out_vcf)
    String out_idx_base = basename(out_idx)
    String out_stats_base = basename(out_stats)

    Int disk_size = ceil((size(reference, "GB") + size(input_vcf, "GB")) * 1.5 + 20)

    command {
        set -eo pipefail

        echo "In FilterMutectCalls"
        env

        my_input_vcf=~{input_vcf}
        my_input_vcf_index=~{input_vcf_index}
        my_mutect_stats=~{mutect_stats}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        my_out_vcf=~{out_vcf_base}
        my_out_idx=~{out_idx_base}
        my_out_stats=~{out_stats_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{input_vcf} ~{input_vcf_index} ~{mutect_stats} ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key} /scratch/
            rm ~{input_vcf} ~{input_vcf_index} ~{mutect_stats} ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key}
            my_input_vcf=/scratch/~{input_vcf_base}
            my_input_vcf_index=/scratch/~{input_vcf_index_base}
            my_mutect_stats=/scratch/~{mutect_stats_base}
            my_reference=/scratch/~{reference_base}
            my_reference_index=/scratch/~{reference_index_base}
            my_reference_dict=/scratch/~{reference_dict_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_vcf=/scratch/~{out_vcf_base}
            my_out_idx=/scratch/~{out_idx_base}
            my_out_stats=/scratch/~{out_stats_base}
        fi

        time \
        gatk --java-options "-Xmx6g" \
            FilterMutectCalls \
            -R $my_reference \
            -V $my_input_vcf \
            --stats $my_mutect_stats \
            ~{ if filter_ob then '-ob-priors ' + select_first([artifact_priors]) else "" } \
            ~{ if filter_cont then '--contamination-table ' + select_first([contamination_table]) else "" } \
            ~{ if filter_cont then '--tumor-segmentation ' + select_first([tumor_segmentation_table]) else "" } \
            -O $my_out_vcf \
            -L ~{interval} \
            --filtering-stats $my_out_stats

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_vcf ~{output_dir}
                cp $my_out_stats ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_vcf_base} ~{output_dir}
                gsutil cp ~{out_stats_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_vcf = "~{out_vcf}"
        String output_idx = "~{out_idx}"
        String output_filtering_stats = "~{out_stats}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 2
        memory: "12 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task NormalizeVCF {
    input {
        File reference
        File reference_index
        File reference_dict
        File input_vcf
        String log_dir
        String task_name
        String out_vcf
        String out_idx

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String input_vcf_base = basename(input_vcf)
    String sa_key_base = basename(sa_key)

    String out_vcf_base = basename(out_vcf)
    String out_idx_base = basename(out_idx)

    Int disk_size = ceil((size(reference, "GB") + size(input_vcf, "GB")) * 1.5 + 20)

    command {
        set -eo pipefail

        echo "In NormalizeVCF"
        env

        my_input_vcf=~{input_vcf}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        my_out_vcf=~{out_vcf_base}
        my_out_idx=~{out_idx_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{input_vcf} ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key} /scratch/
            rm ~{input_vcf} ~{reference} ~{reference_index} ~{reference_dict} ~{sa_key}
            my_input_vcf=/scratch/~{input_vcf_base}
            my_reference=/scratch/~{reference_base}
            my_reference_index=/scratch/~{reference_index_base}
            my_reference_dict=/scratch/~{reference_dict_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_vcf=/scratch/~{out_vcf_base}
            my_out_idx=/scratch/~{out_idx_base}
        fi

        time \
        bcftools norm \
            -m - \
            -f $my_reference \
            -o $my_out_vcf \
            -Oz \
            $my_input_vcf

        tabix $my_out_vcf

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_vcf ~{output_dir}
                cp $my_out_idx ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_vcf_base} ~{output_dir}
                gsutil cp ~{out_idx_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_vcf = "~{out_vcf}"
        String output_idx = "~{out_idx}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 2
        memory: "12 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task Cleanup {
    input {
        String filtered_vcf
        String annotated_vcf
        Array[String] inputs
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String sa_key_base = basename(sa_key)

    command {
        set -eo pipefail

        echo "In NormalizeVCF"
        env

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            rm ~{sa_key}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        for f in ~{sep=" " inputs}; do [ ! -e "$f" ] || rm -v "$f" ; done

        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_vcf = "~{filtered_vcf}"
        String output_annotated_vcf = "~{annotated_vcf}"
    }

    runtime {
        memory: 1
        walltime: "0:30:00"
        task_name: task_name
        log_dir: log_dir
    }
}

task SelectUnfiltered {
    input {
        File reference
        File reference_index
        File reference_dict
        File input_vcf
        File input_vcf_index
        String log_dir
        String task_name
        String out_vcf
    }

    Int disk_size = ceil((size(reference, "GB") + size(input_vcf, "GB")) * 1.5 + 20)

    command {
        time \
        gatk --java-options "-Xmx24g" \
        SelectVariants \
        -R ~{reference} \
        --exclude-filtered \
        -V ~{input_vcf} \
        -O ~{out_vcf}
    }

    output {
        File output_vcf = "~{out_vcf}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.9.0"
        cpu: 4
        memory: "48 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "4:00:00"
        task_name: task_name
        log_dir: log_dir
    }
}

task AnnotateVcf {
    input {
        String input_vcf
        String normal_sample
        String tumor_sample
        String build_version
        String log_dir
        String task_name
        String out_vcf
    }

    String cancer_samples_file = "snpEff_cancer_samples.txt"

    command {
        set -e

        module purge
        . /hpf/largeprojects/adam/local/etc/shlienlab.bashrc
        module load java/1.8.0_161
        module load snpEff/4.3

        echo -e "~{normal_sample}\t~{tumor_sample}" > ~{cancer_samples_file}

        time \
        java -Xmx8g \
        -jar /hpf/tools/centos6/snpEff/4.3/snpEff.jar \
        -c /hpf/largeprojects/adam/local/wdl_pipelines/snpeff_data/snpEff.config \
        -v \
        -cancer \
        -cancerSamples ~{cancer_samples_file} \
        ~{build_version} \
        ~{input_vcf} \
        > ~{out_vcf}
    }

    output {
        String output_vcf = "~{out_vcf}"
    }

    runtime {
        memory: 20
        walltime: "8:00:00"
        task_name: task_name
        log_dir: log_dir
    }
}
