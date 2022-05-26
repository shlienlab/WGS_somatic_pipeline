version 1.0
import "common_tasks.wdl" as CommonTasks

workflow CNV {
    input {
        String output_dir
        String patient
        String sample
        String tumor_bam
        String tumor_bai
        String normal_bam
        String normal_bai
        String sex
        File reference
        File reference_index
        File reference_dict
        File impute_info
        File impute_files_tar
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar

        String repo_dir
        File sa_key
        String environment
    }

    Array[Int] allelecount_chrs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46]
    Array[Int] chrs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

    # Directories
    String log_dir = "~{output_dir}/wdl_logs+run_info"
    String progress_dir = "~{log_dir}/progress"

    # Task Names
    String BICseq_task_name = "~{patient}.~{sample}.BICseq"
    String CheckInputs_task_name = "~{patient}.~{sample}.CheckInputs"
    String AlleleCount_task_name_base = "~{patient}.~{sample}.AlleleCount"
    scatter (chr in allelecount_chrs) {
        String AlleleCount_task_names = "~{AlleleCount_task_name_base}.~{chr}"
    }
    String Baflog_task_name = "~{patient}.~{sample}.Baflog"
    String GCCorrect_task_name = "~{patient}.~{sample}.GCCorrect"
    String ImputeFromAF_task_name_base = "~{patient}.~{sample}.ImputeFromAF"
    String Impute_task_name_base = "~{patient}.~{sample}.Impute"
    String CombineImpute_task_name_base = "~{patient}.~{sample}.CombineImpute"
    String HaplotypeBAFs_task_name_base = "~{patient}.~{sample}.HaplotypeBAFs"
    String CleanupPostBAF_task_name_base = "~{patient}.~{sample}.CleanupPostBAF"
    String PlotHaplotypes_task_name_base = "~{patient}.~{sample}.PlotHaplotypes"
    scatter (chr in chrs) {
        String ImputeFromAF_task_names = "~{ImputeFromAF_task_name_base}.~{chr}"
        String Impute_task_names = "~{Impute_task_name_base}.~{chr}"
        String CombineImpute_task_names = "~{CombineImpute_task_name_base}.~{chr}"
        String HaplotypeBAFs_task_names = "~{HaplotypeBAFs_task_name_base}.~{chr}"
        String CleanupPostBAF_task_names = "~{CleanupPostBAF_task_name_base}.~{chr}"
        String PlotHaplotypes_task_names = "~{PlotHaplotypes_task_name_base}.~{chr}"
    }
    String CombineBAFs_task_name = "~{patient}.~{sample}.CombineBAFs"
    String SegmentPhased_task_name = "~{patient}.~{sample}.SegmentPhased"
    String FitCN_task_name = "~{patient}.~{sample}.FitCN"
    String Subclones_task_name = "~{patient}.~{sample}.Subclones"
    String Finalise_task_name = "~{patient}.~{sample}.Finalise"
    String VcfAscat_task_name = "~{patient}.~{sample}.VcfAscat"
    String GZSubclones_task_name = "~{patient}.~{sample}.GZSubclones"
    Array[String] task_names = flatten([[BICseq_task_name, CheckInputs_task_name], AlleleCount_task_names, [Baflog_task_name, GCCorrect_task_name], ImputeFromAF_task_names, Impute_task_names, CombineImpute_task_names, HaplotypeBAFs_task_names, CleanupPostBAF_task_names, PlotHaplotypes_task_names, [CombineBAFs_task_name, SegmentPhased_task_name, FitCN_task_name, Subclones_task_name, Finalise_task_name, VcfAscat_task_name, GZSubclones_task_name]])

    if (environment=="LOCAL") {
        call CommonTasks.FindWorkDir as CNVFindWorkDir {
            input:
                sa_key = sa_key,
                patient=patient,
                sample=sample,
                my_workflow='CNV',
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

    String BICseq_output_rda = "~{output_dir}/bicseq_data.rda"
    if (CheckProgress.task_completion[BICseq_task_name] == "false") {
        call BICseq {
            input:
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                log_dir = log_dir,
                task_name = BICseq_task_name,
                out_rda = BICseq_output_rda,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String battenberg_dir = "battenberg"
    scatter (chr in allelecount_chrs) {
        String AlleleCount_task_name = "~{AlleleCount_task_name_base}.~{chr}"
        String AlleleCount_out_dir = "~{output_dir}/~{battenberg_dir}.AlleleCount.~{chr}.tar.gz"
        if (CheckProgress.task_completion[AlleleCount_task_name] == "false") {
            call AlleleCount {
                input:
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    sex = sex,
                    reference_index = reference_index,
                    impute_info = impute_info,
                    thougen_loci_tar = thougen_loci_tar,
                    ignore_contig = ignore_contig,
                    prob_loci = prob_loci,
                    gc_correct_tar = gc_correct_tar,
                    chr = chr,
                    log_dir = log_dir,
                    battenberg_dir = battenberg_dir,
                    AlleleCount_output_dir = AlleleCount_out_dir,
                    task_name = AlleleCount_task_name,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }
    }

    String Baflog_out_dir = "~{output_dir}/~{battenberg_dir}.Baflog.tar.gz"
    Array [String] Baflog_files_to_delete = AlleleCount_out_dir
    if (CheckProgress.task_completion[Baflog_task_name] == "false") {
        call Baflog {
            input:
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                sex = sex,
                reference_index = reference_index,
                impute_info = impute_info,
                thougen_loci_tar = thougen_loci_tar,
                ignore_contig = ignore_contig,
                prob_loci = prob_loci,
                gc_correct_tar = gc_correct_tar,
                out_dir_tar = if (length(select_all(AlleleCount.out_dir)) == length(allelecount_chrs)) then select_all(AlleleCount.out_dir) else AlleleCount_out_dir,
                battenberg_dir = battenberg_dir,
                Baflog_output_dir = Baflog_out_dir,
                log_dir = log_dir,
                task_name = Baflog_task_name,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = Baflog_files_to_delete
        }
    }

    String GCCorrect_out_dir = "~{output_dir}/~{battenberg_dir}.GCCorrect.tar.gz"
    Array [String] GCCorrect_files_to_delete = [Baflog_out_dir]
    if (CheckProgress.task_completion[GCCorrect_task_name] == "false") {
        call GCCorrect {
            input:
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                sex = sex,
                reference_index = reference_index,
                impute_info = impute_info,
                thougen_loci_tar = thougen_loci_tar,
                ignore_contig = ignore_contig,
                prob_loci = prob_loci,
                gc_correct_tar = gc_correct_tar,
                out_dir_tar = select_first([Baflog.out_dir, Baflog_out_dir]),
                battenberg_dir = battenberg_dir,
                GCCorrect_output_dir = GCCorrect_out_dir,
                log_dir = log_dir,
                task_name = GCCorrect_task_name,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = GCCorrect_files_to_delete
        }
    }

    scatter (chr in chrs) {
        String ImputeFromAF_task_name = "~{ImputeFromAF_task_name_base}.~{chr}"
        String ImputeFromAF_out_dir = "~{output_dir}/~{battenberg_dir}.ImputeFromAF.~{chr}.tar.gz"
        if (CheckProgress.task_completion[Impute_task_name] == "false") {
            call ImputeFromAF {
                input:
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    sex = sex,
                    reference_index = reference_index,
                    impute_info = impute_info,
                    impute_files_tar = impute_files_tar,
                    thougen_loci_tar = thougen_loci_tar,
                    ignore_contig = ignore_contig,
                    prob_loci = prob_loci,
                    gc_correct_tar = gc_correct_tar,
                    chr = chr,
                    out_dir_tar = select_first([GCCorrect.out_dir, GCCorrect_out_dir]),
                    battenberg_dir = battenberg_dir,
                    ImputeFromAF_output_dir = ImputeFromAF_out_dir,
                    log_dir = log_dir,
                    task_name = ImputeFromAF_task_name,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String Impute_task_name = "~{Impute_task_name_base}.~{chr}"
        String Impute_out_dir = "~{output_dir}/~{battenberg_dir}.Impute.~{chr}.tar.gz"
        Array [String] Impute_files_to_delete = [ImputeFromAF_out_dir]
        if (CheckProgress.task_completion[Impute_task_name] == "false") {
            call Impute {
                input:
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    sex = sex,
                    reference_index = reference_index,
                    impute_info = impute_info,
                    impute_files_tar = impute_files_tar,
                    thougen_loci_tar = thougen_loci_tar,
                    ignore_contig = ignore_contig,
                    prob_loci = prob_loci,
                    gc_correct_tar = gc_correct_tar,
                    chr = chr,
                    out_dir_tar = select_first([ImputeFromAF.out_dir, ImputeFromAF_out_dir]),
                    battenberg_dir = battenberg_dir,
                    Impute_output_dir = Impute_out_dir,
                    log_dir = log_dir,
                    task_name = Impute_task_name,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = Impute_files_to_delete
            }
        }

        String CombineImpute_task_name = "~{CombineImpute_task_name_base}.~{chr}"
        String CombineImpute_out_dir = "~{output_dir}/~{battenberg_dir}.CombineImpute.~{chr}.tar.gz"
        Array [String] CombineImpute_files_to_delete = [Impute_out_dir]
        if (CheckProgress.task_completion[CombineImpute_task_name] == "false") {
            call CombineImpute {
                input:
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    sex = sex,
                    reference_index = reference_index,
                    impute_info = impute_info,
                    impute_files_tar = impute_files_tar,
                    thougen_loci_tar = thougen_loci_tar,
                    ignore_contig = ignore_contig,
                    prob_loci = prob_loci,
                    gc_correct_tar = gc_correct_tar,
                    chr = chr,
                    out_dir_tar = select_first([Impute.out_dir, Impute_out_dir]),
                    battenberg_dir = battenberg_dir,
                    CombineImpute_output_dir = CombineImpute_out_dir,
                    log_dir = log_dir,
                    task_name = CombineImpute_task_name,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = CombineImpute_files_to_delete
            }
        }

        String HaplotypeBAFs_task_name = "~{HaplotypeBAFs_task_name_base}.~{chr}"
        String HaplotypeBAFs_out_dir = "~{output_dir}/~{battenberg_dir}.HaplotypeBAFs.~{chr}.tar.gz"
        Array [String] HaplotypeBAFs_files_to_delete = [CombineImpute_out_dir]
        if (CheckProgress.task_completion[HaplotypeBAFs_task_name] == "false") {
            call HaplotypeBAFs {
                input:
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    sex = sex,
                    reference_index = reference_index,
                    impute_info = impute_info,
                    thougen_loci_tar = thougen_loci_tar,
                    ignore_contig = ignore_contig,
                    prob_loci = prob_loci,
                    gc_correct_tar = gc_correct_tar,
                    chr = chr,
                    out_dir_tar = select_first([CombineImpute.out_dir, CombineImpute_out_dir]),
                    battenberg_dir = battenberg_dir,
                    HaplotypeBAFs_output_dir = HaplotypeBAFs_out_dir,
                    log_dir = log_dir,
                    task_name = HaplotypeBAFs_task_name,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = HaplotypeBAFs_files_to_delete
            }
        }

        String CleanupPostBAF_task_name = "~{CleanupPostBAF_task_name_base}.~{chr}"
        String CleanupPostBAF_out_dir = "~{output_dir}/~{battenberg_dir}.CleanupPostBAF.~{chr}.tar.gz"
        Array [String] CleanupPostBAF_files_to_delete = [HaplotypeBAFs_out_dir]
        if (CheckProgress.task_completion[CleanupPostBAF_task_name] == "false") {
            call CleanupPostBAF {
                input:
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    sex = sex,
                    reference_index = reference_index,
                    impute_info = impute_info,
                    thougen_loci_tar = thougen_loci_tar,
                    ignore_contig = ignore_contig,
                    prob_loci = prob_loci,
                    gc_correct_tar = gc_correct_tar,
                    chr = chr,
                    out_dir_tar = select_first([HaplotypeBAFs.out_dir, HaplotypeBAFs_out_dir]),
                    battenberg_dir = battenberg_dir,
                    CleanupPostBAF_output_dir = CleanupPostBAF_out_dir,
                    log_dir = log_dir,
                    task_name = CleanupPostBAF_task_name,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = CleanupPostBAF_files_to_delete
            }
        }

        String PlotHaplotypes_task_name = "~{PlotHaplotypes_task_name_base}.~{chr}"
        String PlotHaplotypes_out_dir = "~{output_dir}/~{battenberg_dir}.PlotHaplotypes.~{chr}.tar.gz"
        Array [String] PlotHaplotypes_files_to_delete = [CleanupPostBAF_out_dir]
        if (CheckProgress.task_completion[PlotHaplotypes_task_name] == "false") {
            call PlotHaplotypes {
                input:
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    sex = sex,
                    reference_index = reference_index,
                    impute_info = impute_info,
                    thougen_loci_tar = thougen_loci_tar,
                    ignore_contig = ignore_contig,
                    prob_loci = prob_loci,
                    gc_correct_tar = gc_correct_tar,
                    chr = chr,
                    out_dir_tar = select_first([CleanupPostBAF.out_dir, CleanupPostBAF_out_dir]),
                    battenberg_dir = battenberg_dir,
                    PlotHaplotypes_output_dir = PlotHaplotypes_out_dir,
                    log_dir = log_dir,
                    task_name = PlotHaplotypes_task_name,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = PlotHaplotypes_files_to_delete
            }
        }
    }

    String CombineBAFs_out_dir = "~{output_dir}/~{battenberg_dir}.CombineBAFs.tar.gz"
    Array [String] CombineBAFs_files_to_delete = flatten([PlotHaplotypes_out_dir, [GCCorrect_out_dir]])
    if (CheckProgress.task_completion[CombineBAFs_task_name] == "false") {
        call CombineBAFs {
            input:
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                sex = sex,
                reference_index = reference_index,
                impute_info = impute_info,
                thougen_loci_tar = thougen_loci_tar,
                ignore_contig = ignore_contig,
                prob_loci = prob_loci,
                gc_correct_tar = gc_correct_tar,
                out_dir_tar = if (length(select_all(PlotHaplotypes.out_dir)) == length(chrs)) then select_all(PlotHaplotypes.out_dir) else PlotHaplotypes_out_dir,
                battenberg_dir = battenberg_dir,
                CombineBAFs_output_dir = CombineBAFs_out_dir,
                log_dir = log_dir,
                task_name = CombineBAFs_task_name,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = CombineBAFs_files_to_delete
        }
    }

    String SegmentPhased_out_dir = "~{output_dir}/~{battenberg_dir}.SegmentPhased.tar.gz"
    Array [String] SegmentPhased_files_to_delete = [CombineBAFs_out_dir]
    if (CheckProgress.task_completion[SegmentPhased_task_name] == "false") {
        call SegmentPhased {
            input:
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                sex = sex,
                reference_index = reference_index,
                impute_info = impute_info,
                thougen_loci_tar = thougen_loci_tar,
                ignore_contig = ignore_contig,
                prob_loci = prob_loci,
                gc_correct_tar = gc_correct_tar,
                out_dir_tar = select_first([CombineBAFs.out_dir, CombineBAFs_out_dir]),
                battenberg_dir = battenberg_dir,
                SegmentPhased_output_dir = SegmentPhased_out_dir,
                log_dir = log_dir,
                task_name = SegmentPhased_task_name,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = SegmentPhased_files_to_delete
        }
    }

    String FitCN_out_dir = "~{output_dir}/~{battenberg_dir}.FitCN.tar.gz"
    Array [String] FitCN_files_to_delete = [SegmentPhased_out_dir]
    if (CheckProgress.task_completion[FitCN_task_name] == "false") {
        call FitCN {
            input:
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                sex = sex,
                reference_index = reference_index,
                impute_info = impute_info,
                thougen_loci_tar = thougen_loci_tar,
                ignore_contig = ignore_contig,
                prob_loci = prob_loci,
                gc_correct_tar = gc_correct_tar,
                out_dir_tar = select_first([SegmentPhased.out_dir, SegmentPhased_out_dir]),
                battenberg_dir = battenberg_dir,
                FitCN_output_dir = FitCN_out_dir,
                log_dir = log_dir,
                task_name = FitCN_task_name,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = FitCN_files_to_delete
        }
    }

    String Subclones_out_dir = "~{output_dir}/~{battenberg_dir}.Subclones.tar.gz"
    Array [String] Subclones_files_to_delete = [FitCN_out_dir]
    if (CheckProgress.task_completion[Subclones_task_name] == "false") {
        call Subclones {
            input:
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                sex = sex,
                reference_index = reference_index,
                impute_info = impute_info,
                thougen_loci_tar = thougen_loci_tar,
                ignore_contig = ignore_contig,
                prob_loci = prob_loci,
                gc_correct_tar = gc_correct_tar,
                out_dir_tar = select_first([FitCN.out_dir, FitCN_out_dir]),
                battenberg_dir = battenberg_dir,
                Subclones_output_dir = Subclones_out_dir,
                log_dir = log_dir,
                task_name = Subclones_task_name,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = Subclones_files_to_delete
        }
    }

    String Finalise_out_dir = "~{output_dir}/~{battenberg_dir}.Finalise.tar.gz"
    Array [String] Finalise_files_to_delete = [Subclones_out_dir]
    if (CheckProgress.task_completion[Finalise_task_name] == "false") {
        call Finalise {
            input:
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                sex = sex,
                reference = reference,
                reference_dict = reference_dict,
                reference_index = reference_index,
                impute_info = impute_info,
                thougen_loci_tar = thougen_loci_tar,
                ignore_contig = ignore_contig,
                prob_loci = prob_loci,
                gc_correct_tar = gc_correct_tar,
                out_dir_tar = select_first([Subclones.out_dir, Subclones_out_dir]),
                battenberg_dir = battenberg_dir,
                Finalise_output_dir = Finalise_out_dir,
                log_dir = log_dir,
                task_name = Finalise_task_name,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = Finalise_files_to_delete
        }
    }

    String VcfAscat_output_file = "~{output_dir}/~{sample}.copynumber.caveman.csv"
    if (CheckProgress.task_completion[VcfAscat_task_name] == "false") {
        call VcfAscat {
            input:
                out_dir_tar = select_first([Finalise.out_dir, Finalise_out_dir]),
                log_dir = log_dir,
                task_name = VcfAscat_task_name,
                out_file = VcfAscat_output_file,
                VcfAscat_search_dir = battenberg_dir,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String GZSubclones_output_file = "~{output_dir}/~{sample}_subclones.txt"
    if (CheckProgress.task_completion[GZSubclones_task_name] == "false") {
        call GZSubclones {
            input:
                out_dir_tar = select_first([Finalise.out_dir, Finalise_out_dir]),
                log_dir = log_dir,
                task_name = GZSubclones_task_name,
                out_file = GZSubclones_output_file,
                GZSubclones_search_dir = battenberg_dir,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    output {
        File bicseq_data_rda = select_first([BICseq.output_rda, BICseq_output_rda])
        File copynumber_csv = select_first([VcfAscat.output_file, VcfAscat_output_file])
        File subclones_txt = select_first([GZSubclones.output_file, GZSubclones_output_file])
    }
}

task BICseq {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String log_dir
        String task_name
        String out_rda

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String sa_key_base = basename(sa_key)

    String out_rda_base = basename(out_rda)

    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 1.2 + 20)

    command <<<
        set -eo pipefail

        echo "In BICseq"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}
        my_sa_key=~{sa_key}

        my_out_rda=~{out_rda_base}

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
                    cp ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai} /scratch/
                    rm ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai}
                    my_tumor_bam=/scratch/~{tumor_bam_base}
                    my_tumor_bai=/scratch/~{tumor_bai_base}
                    my_normal_bam=/scratch/~{normal_bam_base}
                    my_normal_bai=/scratch/~{normal_bai_base}
                fi
            fi

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_rda=/scratch/~{out_rda_base}

            cd /scratch
        fi

        time \
        Rscript /opt/scripts/run_bicseq.R \
        --tumour $my_tumor_bam \
        --normal $my_normal_bam

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_rda ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_rda_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_rda = "~{out_rda}"
    }

    # n1-highmem-2
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 2
        memory: "13 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task AlleleCount {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        Int chr
        String battenberg_dir
        String AlleleCount_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)

    String AlleleCount_output_dir_base = if (environment == "LOCAL") then AlleleCount_output_dir else basename(AlleleCount_output_dir)

    Int threads = 2
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 1.5 + 50)

    command <<<
        set -eo pipefail

        echo "In AlleleCount"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}
        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}
        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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
                    #cp ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai} /scratch/
                    #rm ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai}
                    #my_tumor_bam=/scratch/~{tumor_bam_base}
                    #my_tumor_bai=/scratch/~{tumor_bai_base}
                    #my_normal_bam=/scratch/~{normal_bam_base}
                    #my_normal_bai=/scratch/~{normal_bai_base}
                fi
            fi

            cp ~{reference_index} ~{impute_info} ~{thougen_loci_tar} ~{ignore_contig} ~{prob_loci} ~{gc_correct_tar} ~{sa_key} /scratch/
            rm ~{reference_index} ~{impute_info} ~{thougen_loci_tar} ~{ignore_contig} ~{prob_loci} ~{gc_correct_tar} ~{sa_key}
            my_reference_index=/scratch/~{reference_index_base}
            my_impute_info=/scratch/~{impute_info_base}
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            my_ignore_contig=/scratch/~{ignore_contig_base}
            my_prob_loci=/scratch/~{prob_loci_base}
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            my_sa_key=/scratch/~{sa_key_base}

            THOUGEN_LOCI=/scratch/$(basename $my_thougen_loci_tar ".tar.gz")
            GC_CORRECT=/scratch/$(basename $my_gc_correct_tar ".tar.gz")

            cd /scratch/
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        mkdir -p ~{battenberg_dir}

        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p allelecount \
            -i ~{chr}

        tar -zcvf ~{AlleleCount_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{AlleleCount_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{AlleleCount_output_dir}"
    }

    # n1-highmem-2
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "13 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task Baflog {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        Array [File] out_dir_tar
        String battenberg_dir
        String Baflog_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)

    String Baflog_output_dir_base = if (environment == "LOCAL") then Baflog_output_dir else basename(Baflog_output_dir)

    Int threads = 8
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In Baflog"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}

        my_sa_key=~{sa_key}

        my_Baflog_output_dir=~{Baflog_output_dir_base}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} ~{impute_info} ~{thougen_loci_tar} ~{ignore_contig} ~{prob_loci} ~{gc_correct_tar} ~{sa_key} /scratch/
            rm ~{reference_index} ~{impute_info} ~{thougen_loci_tar} ~{ignore_contig} ~{prob_loci} ~{gc_correct_tar} ~{sa_key}
            my_reference_index=/scratch/~{reference_index_base}
            my_impute_info=/scratch/~{impute_info_base}
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            my_ignore_contig=/scratch/~{ignore_contig_base}
            my_prob_loci=/scratch/~{prob_loci_base}
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            my_sa_key=/scratch/~{sa_key_base}

            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cd /scratch/
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        for archive in ~{sep=' ' out_dir_tar}; do
            tar -zxvf $archive
        done

        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p baflog

        tar -zcvf ~{Baflog_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_Baflog_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{Baflog_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{Baflog_output_dir}"
    }

    # n1-highmem-8
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "52 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task GCCorrect {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        File out_dir_tar
        String battenberg_dir
        String GCCorrect_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String out_dir_tar_base = basename(out_dir_tar)
    String sa_key_base = basename(sa_key)

    String GCCorrect_output_dir_base = if (environment == "LOCAL") then GCCorrect_output_dir else basename(GCCorrect_output_dir)

    Int threads = 16
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 2 + 50)

    command <<<
        set -eo pipefail

        echo "In GCCorrect"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}
        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}
        my_out_dir_tar=~{out_dir_tar}
        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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
                    #cp ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai} /scratch/
                    #rm ~{tumor_bam} ~{tumor_bai} ~{normal_bam} ~{normal_bai}
                    #my_tumor_bam=/scratch/~{tumor_bam_base}
                    #my_tumor_bai=/scratch/~{tumor_bai_base}
                    #my_normal_bam=/scratch/~{normal_bam_base}
                    #my_normal_bai=/scratch/~{normal_bai_base}
                fi
            fi

            cp ~{reference_index} ~{impute_info} ~{thougen_loci_tar} ~{ignore_contig} ~{prob_loci} ~{gc_correct_tar} ~{out_dir_tar} ~{sa_key} /scratch/
            rm ~{reference_index} ~{impute_info} ~{thougen_loci_tar} ~{ignore_contig} ~{prob_loci} ~{gc_correct_tar} ~{out_dir_tar} ~{sa_key}
            my_reference_index=/scratch/~{reference_index_base}
            my_impute_info=/scratch/~{impute_info_base}
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            my_ignore_contig=/scratch/~{ignore_contig_base}
            my_prob_loci=/scratch/~{prob_loci_base}
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            my_sa_key=/scratch/~{sa_key_base}

            THOUGEN_LOCI=/scratch/$(basename $my_thougen_loci_tar ".tar.gz")
            GC_CORRECT=/scratch/$(basename $my_gc_correct_tar ".tar.gz")

            cd /scratch/
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_out_dir_tar

        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p gccorrect

        tar -zcvf ~{GCCorrect_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{GCCorrect_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{GCCorrect_output_dir}"
    }

    # n1-highmem-16
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "80 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task ImputeFromAF {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File impute_files_tar
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        Int chr
        File out_dir_tar
        String battenberg_dir
        String ImputeFromAF_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String impute_files_base = basename(impute_files_tar, ".tar.gz")

    String ImputeFromAF_output_dir_base = if (environment == "LOCAL") then ImputeFromAF_output_dir else basename(ImputeFromAF_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In ImputeFromAF"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        my_ImputeFromAF_output_dir=~{ImputeFromAF_output_dir_base}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} ~{impute_info} ~{thougen_loci_tar} ~{ignore_contig} ~{prob_loci} ~{gc_correct_tar} ~{out_dir_tar} ~{sa_key} /scratch/
            rm ~{reference_index} ~{impute_info} ~{thougen_loci_tar} ~{ignore_contig} ~{prob_loci} ~{gc_correct_tar} ~{out_dir_tar} ~{sa_key}
            my_reference_index=/scratch/~{reference_index_base}
            my_impute_info=/scratch/~{impute_info_base}
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            my_ignore_contig=/scratch/~{ignore_contig_base}
            my_prob_loci=/scratch/~{prob_loci_base}
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            my_sa_key=/scratch/~{sa_key_base}

            THOUGEN_LOCI=/scratch/$(basename $my_thougen_loci_tar ".tar.gz")
            GC_CORRECT=/scratch/$(basename $my_gc_correct_tar ".tar.gz")

            cd /scratch/
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar


        tar -zxvf ~{impute_files_tar}
        sed "s:/hpf/largeprojects/adam/local/etc/bberg_ref/impute/:$(pwd)/~{impute_files_base}/:g" $my_impute_info > impute_info.tmp.txt
        echo "Currently in $(pwd)..."
        tar -zxvf $my_out_dir_tar

        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $(pwd)/impute_info.tmp.txt \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p imputefromaf \
            -i ~{chr}

        tar -zcvf ~{ImputeFromAF_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp my_ImputeFromAF_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{ImputeFromAF_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{ImputeFromAF_output_dir}"
    }

    # n1-standard-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "15 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task Impute {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File impute_files_tar
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        Int chr
        File out_dir_tar
        String battenberg_dir
        String Impute_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String impute_files_tar_base = basename(impute_files_tar)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String impute_files_base = basename(impute_files_tar, ".tar.gz")

    String Impute_output_dir_base = if (environment == "LOCAL") then Impute_output_dir else basename(Impute_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In Impute"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_impute_files_tar=~{impute_files_tar}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        my_Impute_output_dir=~{Impute_output_dir_base}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{impute_files_tar} /scratch/
            my_impute_files_tar=/scratch/~{impute_files_tar_base}
            rm ~{impute_files_tar}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_impute_files_tar
        sed "s:/hpf/largeprojects/adam/local/etc/bberg_ref/impute/:$(pwd)/~{impute_files_base}/:g" $my_impute_info > impute_info.tmp.txt

        tar -zxvf $my_out_dir_tar

        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $(pwd)/impute_info.tmp.txt \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p impute \
            -i ~{chr}

        tar -zcvf ~{Impute_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_Impute_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{Impute_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{Impute_output_dir}"
    }

    # n1-standard-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "15 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task CombineImpute {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File impute_files_tar
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        Int chr
        File out_dir_tar
        String battenberg_dir
        String CombineImpute_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String impute_files_tar_base = basename(impute_files_tar)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String impute_files_base = basename(impute_files_tar, ".tar.gz")

    String CombineImpute_output_dir_base = if (environment == "LOCAL") then CombineImpute_output_dir else basename(CombineImpute_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In CombineImpute"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_impute_files_tar=~{impute_files_tar}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{impute_files_tar} /scratch/
            my_impute_files_tar=/scratch/~{impute_files_tar_base}
            rm ~{impute_files_tar}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_impute_files_tar
        sed "s:/hpf/largeprojects/adam/local/etc/bberg_ref/impute/:$(pwd)/~{impute_files_base}/:g" $my_impute_info > impute_info.tmp.txt

        tar -zxvf $my_out_dir_tar

        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $(pwd)/impute_info.tmp.txt \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p combineimpute \
            -i ~{chr}

        tar -zcvf ~{CombineImpute_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_CombineImpute_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{CombineImpute_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{CombineImpute_output_dir}"
    }

    # n1-standard-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "15 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task HaplotypeBAFs {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        Int chr
        File out_dir_tar
        String battenberg_dir
        String HaplotypeBAFs_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String HaplotypeBAFs_output_dir_base = if (environment == "LOCAL") then HaplotypeBAFs_output_dir else basename(HaplotypeBAFs_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In HaplotypeBAFs"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_out_dir_tar

        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p haplotypebafs \
            -i ~{chr}

        tar -zcvf ~{HaplotypeBAFs_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_HaplotypeBAFs_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{HaplotypeBAFs_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{HaplotypeBAFs_output_dir}"
    }

    # n1-standard-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "15 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task CleanupPostBAF {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        Int chr
        File out_dir_tar
        String battenberg_dir
        String CleanupPostBAF_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String CleanupPostBAF_output_dir_base = if (environment == "LOCAL") then CleanupPostBAF_output_dir else basename(CleanupPostBAF_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In CleanupPostBAF"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            #my_battenberg_dir=/scratch/~{battenberg_dir}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_out_dir_tar

        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p cleanuppostbaf \
            -i ~{chr}

        tar -zcvf ~{CleanupPostBAF_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_CleanupPostBAF_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{CleanupPostBAF_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{CleanupPostBAF_output_dir}"
    }

    # n1-standard-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "15 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task PlotHaplotypes {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        Int chr
        File out_dir_tar
        String battenberg_dir
        String PlotHaplotypes_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String PlotHaplotypes_output_dir_base = if (environment == "LOCAL") then PlotHaplotypes_output_dir else basename(PlotHaplotypes_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In PlotHaplotypes"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_out_dir_tar

        time \
        battenberg.pl \
        -o ~{battenberg_dir} \
        -tb $my_tumor_bam \
        -nb $my_normal_bam \
        -r $my_reference_index \
        -ge ~{sex} \
        -e $my_impute_info \
        -u $THOUGEN_LOCI \
        -ig $my_ignore_contig \
        -c $my_prob_loci \
        -gc $GC_CORRECT \
        -t ~{threads} \
        -p plothaplotypes \
        -i ~{chr}

        tar -zcvf ~{PlotHaplotypes_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_PlotHaplotypes_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{PlotHaplotypes_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{PlotHaplotypes_output_dir}"
    }

    # n1-standard-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "15 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task CombineBAFs {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        Array [File] out_dir_tar
        String battenberg_dir
        String CombineBAFs_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)

    String CombineBAFs_output_dir_base = if (environment == "LOCAL") then CombineBAFs_output_dir else basename(CombineBAFs_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In CombineBAFs"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}

        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        for archive in ~{sep=' ' out_dir_tar}; do
            tar -zxvf $archive
        done

        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p combinebafs

        tar -zcvf ~{CombineBAFs_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_CombineBAFs_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{CombineBAFs_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{CombineBAFs_output_dir}"
    }

    # n1-highmem-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "26 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task SegmentPhased {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        File out_dir_tar
        String battenberg_dir
        String SegmentPhased_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String SegmentPhased_output_dir_base = if (environment == "LOCAL") then SegmentPhased_output_dir else basename(SegmentPhased_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB") ) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In SegmentPhased"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_out_dir_tar


        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p segmentphased

        tar -zcvf ~{SegmentPhased_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_SegmentPhased_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{SegmentPhased_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{SegmentPhased_output_dir}"
    }

    # n1-highmem-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "26 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task FitCN {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        File out_dir_tar
        String battenberg_dir
        String FitCN_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String FitCN_output_dir_base = if (environment == "LOCAL") then FitCN_output_dir else basename(FitCN_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In FitCN"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_out_dir_tar


        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p fitcn

        tar -zcvf ~{FitCN_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_FitCN_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{FitCN_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{FitCN_output_dir}"
    }

    # n1-highmem-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "26 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task Subclones {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        File out_dir_tar
        String battenberg_dir
        String Subclones_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String Subclones_output_dir_base = if (environment == "LOCAL") then Subclones_output_dir else basename(Subclones_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In Subclones"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_out_dir_tar


        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p subclones

        tar -zcvf ~{Subclones_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_Subclones_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{Subclones_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{Subclones_output_dir}"
    }

    # n1-highmem-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "26 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task Finalise {
    input {
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        String sex
        File reference
        File reference_dict
        File reference_index
        File impute_info
        File thougen_loci_tar
        File ignore_contig
        File prob_loci
        File gc_correct_tar
        File out_dir_tar
        String battenberg_dir
        String Finalise_output_dir
        String log_dir
        String task_name

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_base = basename(reference)
    String reference_dict_base = basename(reference_dict)
    String reference_index_base = basename(reference_index)
    String impute_info_base = basename(impute_info)
    String thougen_loci_tar_base = basename(thougen_loci_tar)
    String ignore_contig_base = basename(ignore_contig)
    String prob_loci_base = basename(prob_loci)
    String gc_correct_tar_base = basename(gc_correct_tar)
    String sa_key_base = basename(sa_key)
    String out_dir_tar_base = basename(out_dir_tar)

    String Finalise_output_dir_base = if (environment == "LOCAL") then Finalise_output_dir else basename(Finalise_output_dir)

    Int threads = 4
    Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB")) * 3 + 50)

    command <<<
        set -eo pipefail

        echo "In Finalise"
        env

        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}

        my_reference=~{reference}
        my_reference_dict=~{reference_dict}
        my_reference_index=~{reference_index}
        my_impute_info=~{impute_info}
        my_thougen_loci_tar=~{thougen_loci_tar}
        my_ignore_contig=~{ignore_contig}
        my_prob_loci=~{prob_loci}
        my_gc_correct_tar=~{gc_correct_tar}

        my_battenberg_dir=~{battenberg_dir}
        my_out_dir_tar=~{out_dir_tar}

        my_sa_key=~{sa_key}

        THOUGEN_LOCI=$(pwd)/$(basename ~{thougen_loci_tar} ".tar.gz")
        GC_CORRECT=$(pwd)/$(basename ~{gc_correct_tar} ".tar.gz")

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

            cp ~{reference} /scratch/
            my_reference=/scratch/~{reference_base}
            rm ~{reference}

            cp ~{reference_dict} /scratch/
            my_reference_dict=/scratch/~{reference_dict_base}
            rm ~{reference_dict}

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{impute_info} /scratch/
            my_impute_info=/scratch/~{impute_info_base}
            rm ~{impute_info}

            cp ~{thougen_loci_tar} /scratch/
            my_thougen_loci_tar=/scratch/~{thougen_loci_tar_base}
            rm ~{thougen_loci_tar}
            THOUGEN_LOCI=/scratch/$(basename ~{thougen_loci_tar} ".tar.gz")

            cp ~{ignore_contig} /scratch/
            my_ignore_contig=/scratch/~{ignore_contig_base}
            rm ~{ignore_contig}

            cp ~{prob_loci} /scratch/
            my_prob_loci=/scratch/~{prob_loci_base}
            rm ~{prob_loci}

            cp ~{gc_correct_tar} /scratch/
            my_gc_correct_tar=/scratch/~{gc_correct_tar_base}
            rm ~{gc_correct_tar}
            GC_CORRECT=/scratch/$(basename ~{gc_correct_tar} ".tar.gz")

            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            cd /scratch
        fi

        export PERL5LIB=/usr/local/lib/perl5

        tar -zxvf $my_thougen_loci_tar
        tar -zxvf $my_gc_correct_tar

        tar -zxvf $my_out_dir_tar

        time \
        battenberg.pl \
            -o ~{battenberg_dir} \
            -tb $my_tumor_bam \
            -nb $my_normal_bam \
            -r $my_reference_index \
            -ge ~{sex} \
            -e $my_impute_info \
            -u $THOUGEN_LOCI \
            -ig $my_ignore_contig \
            -c $my_prob_loci \
            -gc $GC_CORRECT \
            -t ~{threads} \
            -p finalise \
            -ra GRCh37 \
            -rs human \
            -pr WGS

        tar -zcvf ~{Finalise_output_dir_base} ~{battenberg_dir}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_Finalise_output_dir ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
                #rm -rf ~{battenberg_dir}
                #rm -rf $THOUGEN_LOCI
                #rm -rf $GC_CORRECT
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{Finalise_output_dir_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String out_dir = "~{Finalise_output_dir}"
    }

    # n1-highmem-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: threads
        memory: "26 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task VcfAscat {
    input {
        File out_dir_tar
        String log_dir
        String task_name
        String out_file
        String VcfAscat_search_dir

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }
    String out_dir_tar_base = basename(out_dir_tar)
    String sa_key_base = basename(sa_key)

    String out_file_base = basename(out_file)

    Int disk_size = ceil(size(out_dir_tar, "GB") * 3 + 20)

    command <<<
        set -eo pipefail

        echo "In VcfAscat"
        env

        my_out_dir_tar=~{out_dir_tar}
        my_sa_key=~{sa_key}

        my_out_file=~{out_file_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_file=/scratch/~{out_file_base}

            cd /scratch
        fi

        tar -zxvf $my_out_dir_tar

        time \
        find ~{VcfAscat_search_dir} -type f -name '*_battenberg_cn.vcf.gz' \
        | xargs -I {} bb_vcf_to_ascat_cn.pl {} \
        > $my_out_file

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
        String output_file = "~{out_file}"
    }

    # n1-highmem-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: 4
        memory: "26 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}

task GZSubclones {
    input {
        File out_dir_tar
        String log_dir
        String task_name
        String out_file
        String GZSubclones_search_dir

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String out_dir_tar_base = basename(out_dir_tar)
    String sa_key_base = basename(sa_key)

    String out_file_base = basename(out_file)

    Int disk_size = ceil(size(out_dir_tar, "GB") * 3 + 20)

    command <<<
        set -eo pipefail

        echo "In GZSubclones"
        env

        my_out_dir_tar=~{out_dir_tar}
        my_sa_key=~{sa_key}

        my_out_file=~{out_file_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{out_dir_tar} /scratch/
            my_out_dir_tar=/scratch/~{out_dir_tar_base}
            rm ~{out_dir_tar}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_file=/scratch/~{out_file_base}
        fi
        
        tar -zxvf $my_out_dir_tar

        time \
        find ~{GZSubclones_search_dir} -type f -name '*_subclones.txt.gz' \
        | xargs -I {}  gunzip -c {} \
        > $my_out_file

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
        String output_file = "~{out_file}"
    }

    # n1-highmem-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/cgpbattenberg:3.2.2"
        cpu: 4
        memory: "26 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}
