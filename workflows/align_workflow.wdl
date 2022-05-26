version 1.0
import "common_tasks.wdl" as CommonTasks

workflow Align {
    input {
        String output_dir
        String patient
        String sample
        String tumor_normal
        Array[Int] lanes
        Array[Array[String]] fastqs
        Array[String] uuids
        String library
        String platform
        String platform_unit
        String center
        File reference
        File reference_index
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_sa
        File reference_pac
        File known_sites
        File known_sites_index
        File gatk_intervals

        String repo_dir
        File sa_key
        String environment
    }

    # Directories
    String log_dir = "~{output_dir}/wdl_logs+run_info"
    String progress_dir = "~{log_dir}/progress"

    # Task names
    String BWAMem_task_name_base = "~{patient}.~{sample}.bwa"
    String Sort_task_name_base = "~{patient}.~{sample}.sort"
    scatter (lane in lanes) {
        String BWAMem_task_names = "~{BWAMem_task_name_base}.~{lane}"
        String Sort_task_names = "~{Sort_task_name_base}.~{lane}"
    }

    String CreateIntervals_task_name = "~{patient}.~{sample}.create_intervals"
    String MergeAligned_task_name = "~{patient}.~{sample}.merge_aligned"
    String MarkDuplicates_task_name = "~{patient}.~{sample}.markdup"
    String MergeRecalibrated_task_name = "~{patient}.~{sample}.merge_recalibrated"
    String Idxstats_task_name = "~{patient}.~{sample}.idxstats"
    String Barchart_task_name = "~{patient}.~{sample}.barchart"
    String Flagstat_task_name = "~{patient}.~{sample}.flagstat"
    String DepthOfCoverage_task_name = "~{patient}.~{sample}.coverage"
    String CollectInsertSizeMetrics_task_name = "~{patient}.~{sample}.InsertSize"
    String FASTQC_task_name = "~{patient}.~{sample}.fastqc"

    if (environment=="LOCAL") {
        call CommonTasks.FindWorkDir as AlignFindWorkDir {
            input:
                sa_key = sa_key,
                patient=patient,
                sample=sample,
                my_workflow='Align',
                environment=environment,
                log_dir=log_dir,
                output_dir = output_dir
        }
    }

    call createIntervals {
        input:
            reference_dict = reference_dict,
            log_dir = log_dir,
            task_name = "~{patient}.~{sample}.create_intervals",
            output_dir = output_dir,
            environment = environment
    }

    scatter (index in range(length(createIntervals.sequence_grouping))) {
        String RealignerTargetCreator_task_names = "~{patient}.~{sample}.rtc.~{index}"
        String IndelRealigner_task_names = "~{patient}.~{sample}.IndelRealign.~{index}"
        String BaseRecalibrator_task_names = "~{patient}.~{sample}.BaseRecal.~{index}"
        String ApplyBQSR_task_names = "~{patient}.~{sample}.ApplyBQSR.~{index}"
    }

    Array[String] task_names = flatten([BWAMem_task_names, Sort_task_names, [MergeAligned_task_name, MarkDuplicates_task_name], RealignerTargetCreator_task_names, IndelRealigner_task_names, BaseRecalibrator_task_names, ApplyBQSR_task_names, [MergeRecalibrated_task_name, Idxstats_task_name, Barchart_task_name, Flagstat_task_name, DepthOfCoverage_task_name, CollectInsertSizeMetrics_task_name, FASTQC_task_name]])

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

    scatter (lane in lanes) {
        String BWAMem_task_name = "~{BWAMem_task_name_base}.~{lane}"
        String BWAMem_output_bam = "~{output_dir}/~{sample}.~{lane}.bam"

        scatter (fastq in fastqs[lane]) {
            if (environment=="LOCAL") {
                call CommonTasks.WriteManifest{
                    input:
                        sa_key = sa_key,
                        patient = patient,
                        sample = sample,
                        tumor_normal = tumor_normal,
                        file_path = fastq,
                        file_type = 'fastq',
                        environment = environment,
                        task_name = "~{patient}.~{sample}.~{fastq}",
                        log_dir = log_dir
                }
            }
        }

        if (CheckProgress.task_completion[BWAMem_task_name] == "false") {
            call BWAMem {
                input:
                    id = uuids[lane],
                    sample = sample,
                    library = library,
                    platform = platform,
                    platform_unit = platform_unit,
                    center = center,
                    r1_fastq = fastqs[lane][0],
                    r2_fastq = fastqs[lane][1],
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    reference_amb = reference_amb,
                    reference_ann = reference_ann,
                    reference_bwt = reference_bwt,
                    reference_sa = reference_sa,
                    reference_pac = reference_pac,
                    log_dir = log_dir,
                    task_name = BWAMem_task_name,
                    out_bam = BWAMem_output_bam,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String Sort_task_name = "~{Sort_task_name_base}.~{lane}"
        String Sort_input_bam = select_first([BWAMem.output_bam, BWAMem_output_bam])
        String Sort_output_bam = "~{output_dir}/~{sample}.~{lane}.sorted.bam"
        String Sort_output_bai = "~{Sort_output_bam}.bai"
        Array [String] Sort_files_to_delete = [BWAMem_output_bam]
        if (CheckProgress.task_completion[Sort_task_name] == "false") {
            call Sort {
                input:
                    input_bam = Sort_input_bam,
                    log_dir = log_dir,
                    task_name = Sort_task_name,
                    out_bam = Sort_output_bam,
                    out_bai = Sort_output_bai,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = Sort_files_to_delete
            }
        }
    }

    String MergeAligned_output_bam = "~{output_dir}/~{sample}.bam"
    String MergeAligned_output_bai = "~{output_dir}/~{sample}.bai"
    Array [String] MergeAligned_files_to_delete = flatten([Sort_output_bam, Sort_output_bai])
    if (CheckProgress.task_completion[MergeAligned_task_name] == "false") {
        Array[String] MergeAligned_input_bams = if (length(select_all(Sort.output_bam)) == length(lanes)) then select_all(Sort.output_bam) else Sort_output_bam
        Array[String] MergeAligned_input_bais = if (length(select_all(Sort.output_bai)) == length(lanes)) then select_all(Sort.output_bai) else Sort_output_bai
        call MergeAligned{
            input:
                input_bams = MergeAligned_input_bams,
                input_bais = MergeAligned_input_bais,
                reference = reference,
                reference_index = reference_index,
                reference_dict = reference_dict,
                log_dir = log_dir,
                task_name = MergeAligned_task_name,
                out_bam = MergeAligned_output_bam,
                out_bai = MergeAligned_output_bai,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = MergeAligned_files_to_delete
        }
    }

    String MarkDuplicates_output_bam = "~{output_dir}/~{sample}.markdup.bam"
    String MarkDuplicates_output_bai = "~{MarkDuplicates_output_bam}.bai"
    Array [String] MarkDuplicates_files_to_delete = [MergeAligned_output_bam, MergeAligned_output_bai]
    if (CheckProgress.task_completion[MarkDuplicates_task_name] == "false") {
        call MarkDuplicates {
            input:
                input_bam = select_first([MergeAligned.output_bam, MergeAligned_output_bam]),
                input_bai = select_first([MergeAligned.output_bai, MergeAligned_output_bai]),
                log_dir = log_dir,
                task_name = MarkDuplicates_task_name,
                out_bam = MarkDuplicates_output_bam,
                out_bai = MarkDuplicates_output_bai,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = MarkDuplicates_files_to_delete
        }
    }

    scatter (index in range(length(createIntervals.sequence_grouping))) {
        Array [String] interval_list = createIntervals.sequence_grouping[index]

        String RealignerTargetCreator_task_name = "~{patient}.~{sample}.rtc.~{index}"
        String RealignerTargetCreator_output_intervals = "~{output_dir}/~{sample}.markdup.rtc.~{index}.intervals"
        if (CheckProgress.task_completion[RealignerTargetCreator_task_name] == "false") {
            call RealignerTargetCreator {
                input:
                    input_bam = select_first([MarkDuplicates.output_bam, MarkDuplicates_output_bam]),
                    input_bai = select_first([MarkDuplicates.output_bai, MarkDuplicates_output_bai]),
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    interval_list = interval_list,
                    log_dir = log_dir,
                    task_name = RealignerTargetCreator_task_name,
                    out_intervals = RealignerTargetCreator_output_intervals,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String IndelRealigner_task_name = "~{patient}.~{sample}.IndelRealign.~{index}"
        String IndelRealigner_output_bam = "~{output_dir}/~{sample}.markdup.indelrealigned.~{index}.bam"
        String IndelRealigner_output_bai = "~{output_dir}/~{sample}.markdup.indelrealigned.~{index}.bai"
        if (CheckProgress.task_completion[IndelRealigner_task_name] == "false") {
            call IndelRealigner {
                input:
                    input_bam = select_first([MarkDuplicates.output_bam, MarkDuplicates_output_bam]),
                    input_bai = select_first([MarkDuplicates.output_bai, MarkDuplicates_output_bai]),
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    target_intervals = select_first([RealignerTargetCreator.output_intervals, RealignerTargetCreator_output_intervals]),
                    interval_list = interval_list,
                    log_dir = log_dir,
                    task_name = IndelRealigner_task_name,
                    out_bam = IndelRealigner_output_bam,
                    out_bai = IndelRealigner_output_bai,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String BaseRecalibrator_task_name = "~{patient}.~{sample}.BaseRecal.~{index}"
        String BaseRecalibrator_output_table = "~{output_dir}/~{sample}.markdup.indelrealigned.baserecal.~{index}.table"
        if (CheckProgress.task_completion[BaseRecalibrator_task_name] == "false") {
            call BaseRecalibrator {
                input:
                    input_bam = select_first([IndelRealigner.output_bam, IndelRealigner_output_bam]),
                    input_bai = select_first([IndelRealigner.output_bai, IndelRealigner_output_bai]),
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    known_sites = known_sites,
                    known_sites_index = known_sites_index,
                    interval_list = interval_list,
                    log_dir = log_dir,
                    task_name = BaseRecalibrator_task_name,
                    out_table = BaseRecalibrator_output_table,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String ApplyBQSR_task_name = "~{patient}.~{sample}.ApplyBQSR.~{index}"
        String ApplyBQSR_output_bam = "~{output_dir}/~{sample}.realigned-recalibrated.~{index}.bam"
        String ApplyBQSR_output_bai = "~{output_dir}/~{sample}.realigned-recalibrated.~{index}.bai"
        Array [String] ApplyBQSR_files_to_delete = [RealignerTargetCreator_output_intervals, IndelRealigner_output_bam, IndelRealigner_output_bai, BaseRecalibrator_output_table]
        if (CheckProgress.task_completion[ApplyBQSR_task_name] == "false") {
            call ApplyBQSR {
                input:
                    input_bam = select_first([IndelRealigner.output_bam, IndelRealigner_output_bam]),
                    input_bai = select_first([IndelRealigner.output_bai, IndelRealigner_output_bai]),
                    recal_table = select_first([BaseRecalibrator.output_table, BaseRecalibrator_output_table]),
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    interval_list = interval_list,
                    log_dir = log_dir,
                    task_name = ApplyBQSR_task_name,
                    out_bam = ApplyBQSR_output_bam,
                    out_bai = ApplyBQSR_output_bai,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = ApplyBQSR_files_to_delete
            }
        }
    }

    Array[String] MergeRecalibrated_input_bams = if (length(select_all(ApplyBQSR.output_bam)) == length(createIntervals.sequence_grouping)) then select_all(ApplyBQSR.output_bam) else ApplyBQSR_output_bam
    Array[String] MergeRecalibrated_input_bais = if (length(select_all(ApplyBQSR.output_bai)) == length(createIntervals.sequence_grouping)) then select_all(ApplyBQSR.output_bai) else ApplyBQSR_output_bai
    String MergeRecalibrated_output_bam = "~{output_dir}/~{sample}.realigned-recalibrated.bam"
    String MergeRecalibrated_output_bai = "~{output_dir}/~{sample}.realigned-recalibrated.bai"
    String MergeRecalibrated_output_bambai = "~{output_dir}/~{sample}.realigned-recalibrated.bam.bai"
    String MergeRecalibrated_output_md5 = "~{MergeRecalibrated_output_bam}.md5"
    Array [String] MergeRecalibrated_files_to_delete = flatten([ApplyBQSR_output_bam, ApplyBQSR_output_bai])
    if (CheckProgress.task_completion[MergeRecalibrated_task_name] == "false") {
        call MergeRecalibrated {
            input:
                input_bams = MergeRecalibrated_input_bams,
                input_bais = MergeRecalibrated_input_bais,
                log_dir = log_dir,
                task_name = MergeRecalibrated_task_name,
                out_bam = MergeRecalibrated_output_bam,
                out_bai = MergeRecalibrated_output_bai,
                out_bambai = MergeRecalibrated_output_bambai,
                out_md5 = MergeRecalibrated_output_md5,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = MergeRecalibrated_files_to_delete
        }
    }

    String Idxstats_output_file = "~{MergeRecalibrated_output_bam}.idxstats"
    if (CheckProgress.task_completion[Idxstats_task_name] == "false") {
        call Idxstats {
            input:
                bam = select_first([MergeRecalibrated.output_bam, MergeRecalibrated_output_bam]),
                bai = select_first([MergeRecalibrated.output_bai, MergeRecalibrated_output_bai]),
                log_dir = log_dir,
                task_name = Idxstats_task_name,
                out_file = Idxstats_output_file,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    #String Barchart_output_pdf = "~{output_dir}/~{MergeRecalibrated_output_bam}.alignment.barchart.pdf"
    String Barchart_output_pdf = "~{MergeRecalibrated_output_bam}.alignment.barchart.pdf"
    if (CheckProgress.task_completion[Barchart_task_name] == "false") {
        call Barchart {
            input:
                idx_stats = select_first([Idxstats.output_file, Idxstats_output_file]),
                sample = sample,
                log_dir = log_dir,
                task_name = Barchart_task_name,
                out_pdf = Barchart_output_pdf,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    # String Flagstat_output_file = "~{output_dir}/~{MergeRecalibrated_output_bam}.flagstat"
    String Flagstat_output_file = "~{MergeRecalibrated_output_bam}.flagstat"
    if (CheckProgress.task_completion[Flagstat_task_name] == "false") {
        call Flagstat {
            input:
                bam = select_first([MergeRecalibrated.output_bam, MergeRecalibrated_output_bam]),
                bai = select_first([MergeRecalibrated.output_bai, MergeRecalibrated_output_bai]),
                log_dir = log_dir,
                task_name = Flagstat_task_name,
                out_file = Flagstat_output_file,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    # TODO split this task up, it takes > 24h
    String MergeRecalibrated_output_bam_no_ext = sub(MergeRecalibrated_output_bam, ".bam$", "")
    String DepthOfCoverage_output_files_prefix = "~{MergeRecalibrated_output_bam_no_ext}.gatk.depthofcoverage"
    if (CheckProgress.task_completion[DepthOfCoverage_task_name] == "false") {
        call DepthOfCoverage {
            input:
                bam = select_first([MergeRecalibrated.output_bam, MergeRecalibrated_output_bam]),
                bai = select_first([MergeRecalibrated.output_bai, MergeRecalibrated_output_bai]),
                reference = reference,
                reference_index = reference_index,
                reference_dict = reference_dict,
                gatk_intervals = gatk_intervals,
                log_dir = log_dir,
                task_name = DepthOfCoverage_task_name,
                out_prefix = DepthOfCoverage_output_files_prefix,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String CollectInsertSizeMetrics_output_metrics = "~{MergeRecalibrated_output_bam_no_ext}.insertsizemetrics.txt"
    String CollectInsertSizeMetrics_output_histogram = "~{MergeRecalibrated_output_bam_no_ext}.insertsizemetrics.histogram.pdf"
    if (CheckProgress.task_completion[CollectInsertSizeMetrics_task_name] == "false") {
        call CollectInsertSizeMetrics {
            input:
                bam = select_first([MergeRecalibrated.output_bam, MergeRecalibrated_output_bam]),
                bai = select_first([MergeRecalibrated.output_bai, MergeRecalibrated_output_bai]),
                log_dir = log_dir,
                task_name = CollectInsertSizeMetrics_task_name,
                out_metrics = CollectInsertSizeMetrics_output_metrics,
                out_histogram = CollectInsertSizeMetrics_output_histogram,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    String FASTQC_output_zip = "~{MergeRecalibrated_output_bam_no_ext}_fastqc.zip"
    String FASTQC_output_html = "~{MergeRecalibrated_output_bam_no_ext}_fastqc.html"
    if (CheckProgress.task_completion[FASTQC_task_name] == "false") {
        call FASTQC {
            input:
                bam = select_first([MergeRecalibrated.output_bam, MergeRecalibrated_output_bam]),
                log_dir = log_dir,
                task_name = FASTQC_task_name,
                out_zip = FASTQC_output_zip,
                out_html = FASTQC_output_html,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    if (environment=="LOCAL") {
        Array [Pair [String, String]] output_files_types = [
        ("bam", select_first([MergeRecalibrated.output_bam, MergeRecalibrated_output_bam])),
        ("bai", select_first([MergeRecalibrated.output_bai, MergeRecalibrated_output_bai])),
        ("idxstats", select_first([Idxstats.output_file, Idxstats_output_file])),
        ("barchart", select_first([Barchart.output_pdf, Barchart_output_pdf])),
        ("insert_size_metrics", select_first([CollectInsertSizeMetrics.output_metrics, CollectInsertSizeMetrics_output_metrics])),
        ("insert_size_histogram", select_first([CollectInsertSizeMetrics.output_histogram, CollectInsertSizeMetrics_output_histogram])),
        ("fastqc_zip", select_first([FASTQC.output_zip, FASTQC_output_zip])),
        ("fastqc_html", select_first([FASTQC.output_html, FASTQC_output_html])),
        ]

        scatter (f in output_files_types) {
            String file_type = f.left
            String file_path = f.right

            call CommonTasks.WriteManifest as WriteOutputManifest {
                input:
                    sa_key = sa_key,
                    patient = patient,
                    sample = sample,
                    tumor_normal = tumor_normal,
                    file_path = file_path,
                    file_type = file_type,
                    environment = environment,
                    task_name = "~{patient}.~{sample}.~{file_path}",
                    log_dir = log_dir
            }
        }
    }

    output {
        String bam = select_first([MergeRecalibrated.output_bam, MergeRecalibrated_output_bam])
        String bai = select_first([MergeRecalibrated.output_bai, MergeRecalibrated_output_bai])
        String markdup_bam = select_first([MarkDuplicates.output_bam, MarkDuplicates_output_bam])
        String markdup_bai = select_first([MarkDuplicates.output_bai, MarkDuplicates_output_bai])
        String bambai = select_first([MergeRecalibrated.output_bambai, MergeRecalibrated_output_bambai])
        String bam_md5sum = select_first([MergeRecalibrated.output_md5, MergeRecalibrated_output_md5])
        String idxstats = select_first([Idxstats.output_file, Idxstats_output_file])
        String barchart = select_first([Barchart.output_pdf, Barchart_output_pdf])
        String flagstat = select_first([Flagstat.output_file, Flagstat_output_file])
        String coverage_files_prefix = select_first([DepthOfCoverage.output_files_prefix, DepthOfCoverage_output_files_prefix])
        String insert_size_metrics = select_first([CollectInsertSizeMetrics.output_metrics, CollectInsertSizeMetrics_output_metrics])
        String insert_size_histogram = select_first([CollectInsertSizeMetrics.output_histogram, CollectInsertSizeMetrics_output_histogram])
        String fastqc_zip = select_first([FASTQC.output_zip, FASTQC_output_zip])
        String fastqc_html = select_first([FASTQC.output_html, FASTQC_output_html])
        Pair[String, Pair[String, String]] sample_bam_pair = (sample, (select_first([MergeRecalibrated.output_bam, MergeRecalibrated_output_bam]), select_first([MarkDuplicates.output_bam, MarkDuplicates_output_bam])))
    }
}

task BWAMem {
    input {
        String id
        String sample
        String library
        String platform
        String platform_unit
        String center
        File r1_fastq
        File r2_fastq
        File reference
        File reference_index
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_sa
        File reference_pac
        Int? num_threads
        String log_dir
        String task_name
        String out_bam

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String r1_fastq_base = basename(r1_fastq)
    String r2_fastq_base = basename(r2_fastq)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String reference_amb_base = basename(reference_amb)
    String reference_ann_base = basename(reference_ann)
    String reference_bwt_base = basename(reference_bwt)
    String reference_sa_base = basename(reference_sa)
    String reference_pac_base = basename(reference_pac)
    String sa_key_base = basename(sa_key)

    String out_bam_base = basename(out_bam)

    Int n_threads = if (defined(num_threads)) then select_first([num_threads]) else 32
    Int disk_size = ceil((size(r1_fastq, "GB") + size(r2_fastq, "GB") + size(reference, "GB")) * 2 + 50)

    command <<<
        set -eo pipefail

        echo "In BWAMem"
        env

        cd /scratch

        echo $PWD

        # Bring the File variables into bash
        my_r1_fastq=~{r1_fastq}
        my_r2_fastq=~{r2_fastq}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_reference_amb=~{reference_amb}
        my_reference_ann=~{reference_ann}
        my_reference_bwt=~{reference_bwt}
        my_reference_sa=~{reference_sa}
        my_reference_pac=~{reference_pac}
        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then

            # On local environment copy the data files to the scratch space
            cp ~{r1_fastq} /scratch/
            # redefine path to the data to the local scratc hspace copy
            my_r1_fastq=/scratch/~{r1_fastq_base}
            # remove the file hard linked to the cromwell execution directory
            rm ~{r1_fastq}
            echo "done fastq1"

            cp ~{r2_fastq} /scratch/
            my_r2_fastq=/scratch/~{r2_fastq_base}
            rm ~{r2_fastq}
            echo "done fastq2"

            cp ~{reference} /scratch/
            my_reference=/scratch/~{reference_base}
            rm ~{reference}
            echo "done reference"

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}
            echo "done my_reference_index"

            cp ~{reference_dict} /scratch/
            my_reference_dict=/scratch/~{reference_dict_base}
            rm ~{reference_dict}
            echo "done my_reference_dict"

            cp ~{reference_amb} /scratch/
            my_reference_amb=/scratch/~{reference_amb_base}
            rm ~{reference_amb}
            echo "done my_reference_amb"

            cp ~{reference_ann} /scratch/
            my_reference_ann=/scratch/~{reference_ann_base}
            rm ~{reference_ann}
            echo "done my_reference_ann"

            cp ~{reference_bwt} /scratch/
            my_reference_bwt=/scratch/~{reference_bwt_base}
            rm ~{reference_bwt}
            echo "my_reference_bwt"

            cp ~{reference_sa} /scratch/
            my_reference_sa=/scratch/~{reference_sa_base}
            rm ~{reference_sa}
            echo "done my_reference_sa"

            cp ~{reference_pac} /scratch/
            my_reference_pac=/scratch/~{reference_pac_base}
            rm ~{reference_pac}
            echo "done my_reference_pac"

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}
            echo "done fastq2"

        fi

        ls -l /scratch

        #ref_dir=$(dirname ~{reference})
        #mv ~{reference_index} ~{reference_dict} ~{reference_amb} ~{reference_ann} ~{reference_bwt} ~{reference_sa} ~{reference_pac} $ref_dir

        time \
        bwa mem \
            -t ~{n_threads} \
            -R '@RG\tID:~{id}\tSM:~{sample}\tLB:~{library}\tPL:~{platform}\tPU:~{platform_unit}\tCN:~{center}' \
            $my_reference \
            $my_r1_fastq \
            $my_r2_fastq \
        | sambamba view \
            -S \
            -f bam \
            -o ~{out_bam_base} \
            /dev/stdin

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp ~{out_bam_base} ~{out_bam}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_bam_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_bam = "~{out_bam}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/bwa-sambamba:0.7.8"
        cpu: n_threads
        memory: "57.6 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task Sort {
    input {
        File input_bam
        Int? num_threads
        String log_dir
        String task_name
        String out_bam
        String out_bai

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String input_bam_base = basename(input_bam)
    String sa_key_base = basename(sa_key)

    String out_bam_base = basename(out_bam)
    String out_bai_base = basename(out_bai)

    Int n_threads = if (defined(num_threads)) then select_first([num_threads]) else 8
    Int disk_size = ceil(size(input_bam, "GB")*2.5 + 50)

    command {
        set -eo pipefail

        echo "In Sort"
        env

        my_input_bam=~{input_bam}
        my_sa_key=~{sa_key}
        my_out_bam=~{out_bam_base}
        my_out_bai=~{out_bai_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{input_bam} /scratch/
            my_input_bam=/scratch/~{input_bam_base}
            rm ~{input_bam}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key}
            rm ~{sa_key}

            my_out_bam=/scratch/~{out_bam_base}
            my_out_bai=/scratch/~{out_bai_base}
        fi

        time \
        sambamba sort \
            -t ~{n_threads} \
            -o $my_out_bam \
            $my_input_bam

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_bam ~{out_bam}
                cp $my_out_bai ~{out_bai}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_bam_base} ~{output_dir}
                gsutil cp ~{out_bai_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_bam = "~{out_bam}"
        String output_bai = "~{out_bai}"
    }

    # n1-standard-8; might need 32 GB memory
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/bwa-sambamba:0.7.8"
        cpu: n_threads
        memory: "30 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "80:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task MergeAligned {
    input {
        Array[File] input_bams
        Array[File] input_bais
        File reference
        File reference_index
        File reference_dict
        String log_dir
        String task_name
        String out_bam
        String out_bai

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String sa_key_base = basename(sa_key)


    String out_bam_base = basename(out_bam)
    String out_bai_base = basename(out_bai)

    Int disk_size = ceil(size(input_bams[0], "GB") * length(input_bams) * 2.5 + 50)

    command <<<
        set -eo pipefail

        echo "In MergeAligned"
        env

        my_input_bam_string="~{sep=' -I ' input_bams}"

        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        my_out_bam=~{out_bam_base}
        my_out_bai=~{out_bai_base}

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
                    #cp ~{reference} ~{reference_index} ~{reference_dict} /scratch/
                    #rm ~{reference} ~{reference_index} ~{reference_dict}

                    #my_reference=/scratch/~{reference_base}
                    #my_reference_index=/scratch/~{reference_index_base}
                    #my_reference_dict=/scratch/~{reference_dict_base}
                fi
            fi

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}
            echo "done my_sa_key"

            #echo ~{sep=' -I ' input_bams}
            #cp ~{sep=' ' input_bams} /scratch/
            #rm ~{sep=' ' input_bams}
            #rm ~{sep=' ' input_bais}
            #ls -li /scratch/
            #file_list=(/scratch/*.bam)
            #echo "${file_list[@]}"
            #echo "${file_list[*]}"
            #separator=" -I "
            #regex="$( printf "${separator}%s" "${file_list[@]}" )"
            #my_input_bam_string="${regex:${#separator}}"
            #echo my_input_bam_string is: $my_input_bam_string

            my_out_bam=~{out_bam}
            my_out_bai=~{out_bai}
        fi

        #ref_dir=$(dirname ~{reference})
        #mv ~{reference_index} ~{reference_dict} $ref_dir

        time \
        gatk \
            --java-options "-Xmx20g" \
            MergeSamFiles \
            -O /dev/stdout \
            -I $my_input_bam_string \
            --ASSUME_SORTED true \
            --SORT_ORDER coordinate \
        | gatk \
            --java-options "-Xmx20g" \
            SetNmMdAndUqTags \
            -I /dev/stdin \
            -O $my_out_bam \
            -R $my_reference \
            --CREATE_INDEX true

        ls -li /scratch/*

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_out_bam ~{out_bam}
                #cp $my_out_bai ~{out_bai}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_bam_base} ~{output_dir}
                gsutil cp ~{out_bai_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_bam = "~{out_bam}"
        String output_bai = "~{out_bai}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.3.0"
        cpu: 2
        memory: "80 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "80:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task MarkDuplicates {
    input {
        File input_bam
        File input_bai
        Int? num_threads
        String log_dir
        String task_name
        String out_bam
        String out_bai

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String input_bam_base = basename(input_bam)
    String input_bai_base = basename(input_bai)
    String sa_key_base = basename(sa_key)

    String out_bam_base = basename(out_bam)
    String out_bai_base = basename(out_bai)

    Int n_threads = if (defined(num_threads)) then select_first([num_threads]) else 8
    Int disk_size = ceil(size(input_bam, "GB") * 2 + 20)

    command <<<
        set -eo pipefail

        echo "In MarkDuplicates"
        env

        my_input_bam=~{input_bam}
        my_input_bai=~{input_bai}
        my_sa_key=~{sa_key}
        my_out_bam=~{out_bam_base}
        my_out_bai=~{out_bai_base}

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
                    #cp ~{input_bam} ~{input_bai} /scratch/
                    #rm ~{input_bam} ~{input_bai}
                    #my_input_bam=/scratch/~{input_bam_base}
                    #my_input_bai=/scratch/~{input_bai_base}
                fi
            fi

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_bam=/scratch/~{out_bam_base}
            my_out_bai=/scratch/~{out_bai_base}
        fi

        time \
        sambamba markdup \
            -t ~{n_threads} \
            $my_input_bam \
            $my_out_bam

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_bam ~{out_bam}
                cp $my_out_bai ~{out_bai}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_bam_base} ~{output_dir}
                gsutil cp ~{out_bai_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_bam = "~{out_bam}"
        String output_bai = "~{out_bai}"
    }

    # n1-standard-8
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/bwa-sambamba:0.7.8"
        cpu: n_threads
        memory: "45 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "72:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task createIntervals {
    input {
        File reference_dict
        String log_dir
        String task_name
        String output_dir
        String environment
    }

    String reference_dict_base = basename(reference_dict)

    command <<<
        set -eo pipefail

        my_reference_dict=~{reference_dict}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{reference_dict} /scratch/
            my_reference_dict=/scratch/~{reference_dict_base}
            rm ~{reference_dict}
        fi

        create_sequence_intervals.py \
            --dict $my_reference_dict
    >>>

    output {
        Array [Array [String]] sequence_grouping = read_tsv(stdout())
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/dnastack_toolkit:0.0.3"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 HDD"
        preemptible: 2
        walltime: "01:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

# TODO no -nt ~{n_threads} here?
task RealignerTargetCreator {
    input {
        File input_bam
        File input_bai
        File reference
        File reference_index
        File reference_dict
        Array [String] interval_list
        Int n_threads = 4
        String log_dir
        String task_name
        String out_intervals

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String input_bam_base = basename(input_bam)
    String input_bai_base = basename(input_bai)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String sa_key_base = basename(sa_key)

    String out_intervals_base = basename(out_intervals)

    Int disk_size = ceil((size(input_bam, "GB") + size(reference, "GB")) * 3 + 100)

    command <<<
        set -eo pipefail

        echo "In RealignerTargetCreator"
        env

        my_input_bam=~{input_bam}
        my_input_bai=~{input_bai}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        my_out_intervals=~{out_intervals_base}

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
                    #cp ~{input_bam} ~{input_bai} /scratch/
                    #rm ~{input_bam} ~{input_bai}
                    #my_input_bam=/scratch/~{input_bam_base}
                    #my_input_bai=/scratch/~{input_bai_base}
                fi
            fi

            cp ~{reference} /scratch/
            my_reference=/scratch/~{reference_base}
            rm ~{reference}

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{reference_dict} /scratch/
            my_reference_dict=/scratch/~{reference_dict_base}
            rm ~{reference_dict}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_intervals=/scratch/~{out_intervals_base}
        fi

        #ref_dir=$(dirname ~{reference})
        #mv ~{reference_index} ~{reference_dict} $ref_dir

        time \
        java \
            -Xmx20g \
            -jar "$GATK" \
            -T RealignerTargetCreator \
            -I $my_input_bam \
            -R $my_reference \
            -L ~{sep=" -L " interval_list} \
            -l INFO \
            -o $my_out_intervals

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_intervals ~{out_intervals}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_intervals_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_intervals = "~{out_intervals}"
    }

    # n1-highmem-4
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:3.8"
        cpu: n_threads
        memory: "26 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "100:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task IndelRealigner {
    input {
        File input_bam
        File input_bai
        File reference
        File reference_index
        File reference_dict
        File target_intervals
        Array [String] interval_list
        String log_dir
        String task_name
        String out_bam
        String out_bai

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String input_bam_base = basename(input_bam)
    String input_bai_base = basename(input_bai)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String target_intervals_base = basename(target_intervals)
    String sa_key_base = basename(sa_key)

    String out_bam_base = basename(out_bam)
    String out_bai_base = basename(out_bai)

    Int disk_size = ceil((size(input_bam, "GB") + size(reference, "GB")) * 4 + 200)

    command <<<
        set -eo pipefail

        echo "In IndelRealigner"
        env

        my_input_bam=~{input_bam}
        my_input_bai=~{input_bai}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_target_intervals=~{target_intervals}
        my_sa_key=~{sa_key}

        my_out_bam=~{out_bam_base}
        my_out_bai=~{out_bai_base}


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
                    #echo "Skipping moving bams for this task, it seems to slow down this task."
                    cp ~{input_bam} ~{input_bai} /scratch/
                    rm ~{input_bam} ~{input_bai}
                    my_input_bam=/scratch/~{input_bam_base}
                    my_input_bai=/scratch/~{input_bai_base}
                fi
            fi

            cp ~{reference} /scratch/
            my_reference=/scratch/~{reference_base}
            rm ~{reference}

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{reference_dict} /scratch/
            my_reference_dict=/scratch/~{reference_dict_base}
            rm ~{reference_dict}

            cp ~{target_intervals} /scratch/
            my_target_intervals=/scratch/~{target_intervals_base}
            rm ~{target_intervals}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_bam=/scratch/~{out_bam_base}
            my_out_bai=/scratch/~{out_bai_base}
        fi

        #ref_dir=$(dirname ~{reference})
        #mv ~{reference_index} ~{reference_dict} $ref_dir

        time \
        java \
            -Xmx24g \
            -jar "$GATK" \
            -T IndelRealigner \
            -I $my_input_bam \
            -R $my_reference \
            -L ~{sep=" -L " interval_list} \
            -targetIntervals $my_target_intervals \
            -o $my_out_bam

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_bam ~{out_bam}
                cp $my_out_bai ~{out_bai}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_bam_base} ~{output_dir}
                gsutil cp ~{out_bai_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_bam = "~{out_bam}"
        String output_bai = "~{out_bai}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:3.8"
        cpu: 4
        memory: "50 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        walltime: "100:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

# TODO no more -nct ~{n_threads} here?
task BaseRecalibrator {
    input {
        File input_bam
        File input_bai
        File reference
        File reference_index
        File reference_dict
        Array [String] interval_list
        File known_sites
        File known_sites_index
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
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String known_sites_base = basename(known_sites)
    String known_sites_index_base = basename(known_sites_index)
    String sa_key_base = basename(sa_key)

    String out_table_base = basename(out_table)

    Int disk_size = ceil((size(input_bam, "GB") + size(reference, "GB")) * 3 + 100)

    command <<<
        set -eo pipefail

        echo "In BaseRecalibrator"
        env

        my_input_bam=~{input_bam}
        my_input_bai=~{input_bai}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        known_sites=~{known_sites}
        known_sites_index=~{known_sites_index}
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
                    #echo "Skipping moving bams for this task, it seems to slow down this task."
                    cp ~{input_bam} ~{input_bai} /scratch/
                    rm ~{input_bam} ~{input_bai}
                    my_input_bam=/scratch/~{input_bam_base}
                    my_input_bai=/scratch/~{input_bai_base}
                fi
            fi

            cp ~{reference} /scratch/
            my_reference=/scratch/~{reference_base}
            rm ~{reference}

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{reference_dict} /scratch/
            my_reference_dict=/scratch/~{reference_dict_base}
            rm ~{reference_dict}

            cp ~{known_sites} /scratch/
            my_known_sites=/scratch/~{known_sites_base}
            rm ~{known_sites}

            cp ~{known_sites_index} /scratch/
            my_known_sites_index=/scratch/~{known_sites_index_base}
            rm ~{known_sites_index}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_table=/scratch/~{out_table_base}
        fi

        #known_sites_dir=$(dirname ~{known_sites})
        #mv ~{known_sites_index} $known_sites_dir

        #ref_dir=$(dirname ~{reference})
        #mv ~{reference_index} ~{reference_dict} $ref_dir

        time \
        gatk \
            --java-options "-Xmx24g" \
            BaseRecalibrator \
            -R $my_reference \
            -I $my_input_bam \
            -L ~{sep=" -L " interval_list} \
            -O $my_out_table \
            --verbosity INFO \
            --known-sites $my_known_sites

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_table ~{out_table}
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

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.3.0"
        cpu: 8
        memory: "30 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        walltime: "100:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task ApplyBQSR {
    input {
        File input_bam
        File input_bai
        File recal_table
        File reference
        File reference_index
        File reference_dict
        Array [String] interval_list
        String log_dir
        String task_name
        String out_bam
        String out_bai

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String input_bam_base = basename(input_bam)
    String input_bai_base = basename(input_bai)
    String recal_table_base = basename(recal_table)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String sa_key_base = basename(sa_key)

    String out_bam_base = basename(out_bam)
    String out_bai_base = basename(out_bai)

    Int disk_size = ceil((size(input_bam, "GB") + size(reference, "GB")) * 2 + 100)

    command <<<
        set -eo pipefail

        echo "In ApplyBQSR"
        env

        my_input_bam=~{input_bam}
        my_input_bai=~{input_bai}
        my_recal_table=~{recal_table}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_sa_key=~{sa_key}

        my_out_bam=~{out_bam_base}
        my_out_bai=~{out_bai_base}


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
                    #echo "Skipping moving bams for this task, it seems to slow down this task."
                    cp ~{input_bam} ~{input_bai} /scratch/
                    rm ~{input_bam} ~{input_bai}
                    my_input_bam=/scratch/~{input_bam_base}
                    my_input_bai=/scratch/~{input_bai_base}
                fi
            fi

            cp ~{reference} /scratch/
            my_reference=/scratch/~{reference_base}
            rm ~{reference}

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{reference_dict} /scratch/
            my_reference_dict=/scratch/~{reference_dict_base}
            rm ~{reference_dict}

            cp ~{recal_table} /scratch/
            my_recal_table=/scratch/~{recal_table_base}
            rm ~{recal_table}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_bam=/scratch/~{out_bam_base}
            my_out_bai=/scratch/~{out_bai_base}
        fi

        #ref_dir=$(dirname ~{reference})
        #mv ~{reference_index} ~{reference_dict} $ref_dir

        time \
        gatk \
            --java-options "-Xmx24g" \
            ApplyBQSR \
            -R $my_reference \
            -I $my_input_bam \
            -L ~{sep=" -L " interval_list} \
            -O $my_out_bam \
            -bqsr $my_recal_table \
            --verbosity INFO \
            -add-output-sam-program-record

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_bam ~{out_bam}
                cp $my_out_bai ~{out_bai}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_bam_base} ~{output_dir}
                gsutil cp ~{out_bai_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_bam = "~{out_bam}"
        String output_bai = "~{out_bai}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.3.0"
        cpu: 4
        memory: "50 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "100:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

# TODO I think this task is using way more resources than it needs
task MergeRecalibrated {
    input {
        Array[File] input_bams
        Array[File] input_bais
        String log_dir
        String task_name
        String out_bam
        String out_bai
        String out_bambai
        String out_md5

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String sa_key_base = basename(sa_key)

    String out_bam_base = basename(out_bam)
    String out_bai_base = basename(out_bai)
    String out_bambai_base = basename(out_bambai)
    String out_md5_base = basename(out_md5)

    Int disk_size = ceil(size(input_bams[0], "GB") * length(input_bams) * 3 + 100)

    command <<<
        set -eo pipefail

        echo "In MergeRecalibrated"
        env

        my_input_bam_string="~{sep=' -I ' input_bams}"
        my_sa_key=~{sa_key}

        my_out_bam=~{out_bam_base}
        my_out_bai=~{out_bai_base}
        my_out_bambai=~{out_bambai_base}
        my_out_md5=~{out_md5_base}
        echo "my_out_md5: Given md5 is ~{out_md5}; base is $my_out_md5"

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
                    cp ~{sa_key} /scratch/
                    rm ~{sa_key}
                    #my_sa_key=/scratch/~{sa_key_base}
                fi
            fi

            # On local environment copy the data files to the scratch space

            #echo ~{sep=' -I ' input_bams}
            #cp ~{sep=' ' input_bams} /scratch/
            #rm ~{sep=' ' input_bams}
            #rm ~{sep=' ' input_bais}
            #ls -li /scratch/
            #file_list=(/scratch/*.bam)
            #echo "${file_list[@]}"
            #echo "${file_list[*]}"
            #separator=" -I "
            #regex="$( printf "${separator}%s" "${file_list[@]}" )"
            #my_input_bam_string="${regex:${#separator}}"
            #echo my_input_bam_string is: $my_input_bam_string

            my_out_bam=~{out_bam}
            my_out_bai=~{out_bai}
            my_out_bambai=~{out_bambai}
            my_out_md5=~{out_md5}
            echo "my_out_md5: changed to $my_out_md5"
        fi

        time \
        gatk \
            --java-options "-Xmx20g" \
            MergeSamFiles \
            -O $my_out_bam \
            -I $my_input_bam_string \
            --ASSUME_SORTED true \
            --SORT_ORDER coordinate \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --USE_THREADING true

        # The CNV workflow requires the index be in the format .bai, not .bam.bai; we will output both
        #cp ~{out_bai_base} ~{sub(out_bai_base, ".bam.bai$", ".bai")}
        cp $my_out_bai $my_out_bambai

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_out_bam ~{out_bam}
                #cp $my_out_bai ~{out_bai}
                #cp $my_out_bambai ~{out_bambai}
                #cp $my_out_md5 ~{out_md5}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_bam_base} ~{output_dir}
                gsutil cp ~{out_bai_base} ~{output_dir}
                gsutil cp ~{sub(out_bai_base, ".bai$", ".bam.bai")} ~{output_dir}
                gsutil cp ~{out_md5_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_bam = "~{out_bam}"
        String output_bai = "~{out_bai}"
        String output_bambai = "~{out_bambai}"
        String output_md5 = "~{out_md5}"
    }

    # n1-standard-16; not sure this much memory is required (cpu fit to memory)
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:4.1.3.0"
        cpu: 16
        memory: "56 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "80:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task Idxstats {
    input {
        File bam
        File bai
        String log_dir
        String task_name
        String out_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
     }

    String bam_base = basename(bam)
    String bai_base = basename(bai)
    String sa_key_base = basename(sa_key)

    String out_file_base = basename(out_file)

    Int disk_size = ceil(size(bam, "GB") * 1.5 + 20)

    command <<<
        set -eo pipefail

        echo "In Idxstats"
        env

        my_bam=~{bam}
        my_bai=~{bai}
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
                    #cp ~{bam} ~{bai} /scratch/
                    #rm ~{bam} ~{bai}
                    #my_bam=/scratch/~{bam_base}
                    #my_bai=/scratch/~{bai_base}
                fi
            fi

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_file=/scratch/~{out_file_base}
        fi

        time \
        samtools idxstats \
            $my_bam \
            > $my_out_file

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_file ~{out_file}
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

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/samtools:1.9"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task Barchart {
    input {
        File idx_stats
        String sample
        String log_dir
        String task_name
        String out_pdf

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String idx_stats_base = basename(idx_stats)
    String sa_key_base = basename(sa_key)

    String out_pdf_base = basename(out_pdf)

    command {
        set -eo pipefail

        echo "In Barchart"
        env

        my_idx_stats=~{idx_stats_base}
        my_sa_key=~{sa_key}

        my_out_pdf=~{out_pdf_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{idx_stats} /scratch/
            my_idx_stats=/scratch/~{idx_stats_base}
            rm ~{idx_stats}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_pdf=/scratch/~{out_pdf_base}
        fi

        time \
        Rscript /opt/scripts/barchart.R \
            --input $my_idx_stats \
            --sample ~{sample} \
            --output $my_out_pdf

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_pdf ~{out_pdf}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_pdf_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_pdf = "~{out_pdf}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 HDD"
        preemptible: 2
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task Flagstat {
    input {
        File bam
        File bai
        String log_dir
        String task_name
        String out_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String bam_base = basename(bam)
    String bai_base = basename(bai)
    String sa_key_base = basename(sa_key)

    String out_file_base = basename(out_file)

    Int disk_size = ceil(size(bam, "GB") * 1.5 + 20)

    command <<<
        set -eo pipefail

        echo "In Flagstat"
        env

        my_bam=~{bam}
        my_bai=~{bai}
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
                    #cp ~{bam} ~{bai} /scratch/
                    #rm ~{bam} ~{bai}
                    #my_bam=/scratch/~{bam_base}
                    #my_bai=/scratch/~{bai_base}
                fi
            fi

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_file=/scratch/~{out_file_base}
        fi

        time \
        sambamba flagstat \
            $my_bam \
            > $my_out_file

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_file ~{out_file}
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

    # TODO may need more memory/cores (samtools didn't take advantage of it, sambamba m)
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/bwa-sambamba:0.7.8"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task DepthOfCoverage {
    input {
        File bam
        File bai
        File reference
        File reference_index
        File reference_dict
        File gatk_intervals
        String log_dir
        String task_name
        String out_prefix

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String bam_base = basename(bam)
    String bai_base = basename(bai)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String gatk_intervals_base = basename(gatk_intervals)
    String sa_key_base = basename(sa_key)

    String out_prefix_base = basename(out_prefix)

    Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 1.5 + 20)

    command <<<
        set -eo pipefail

        echo "In DepthOfCoverage"
        env

        my_bam=~{bam}
        my_bai=~{bai}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_gatk_intervals=~{gatk_intervals}
        my_sa_key=~{sa_key}

        my_out_prefix=~{out_prefix_base}

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
                    #cp ~{bam} ~{bai} /scratch/
                    #rm ~{bam} ~{bai}
                    #my_bam=/scratch/~{bam_base}
                    #my_bai=/scratch/~{bai_base}
                fi
            fi

            cp ~{reference} /scratch/
            my_reference=/scratch/~{reference_base}
            rm ~{reference}

            cp ~{reference_index} /scratch/
            my_reference_index=/scratch/~{reference_index_base}
            rm ~{reference_index}

            cp ~{reference_dict} /scratch/
            my_reference_dict=/scratch/~{reference_dict_base}
            rm ~{reference_dict}

            cp ~{gatk_intervals} /scratch/
            my_gatk_intervals=/scratch/~{gatk_intervals_base}
            rm ~{gatk_intervals}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_prefix=/scratch/~{out_prefix_base}
        fi

        #ref_dir=$(dirname ~{reference})
        #mv ~{reference_index} ~{reference_dict} $ref_dir

        time \
        java -Xmx16g -jar "$GATK" \
            -T DepthOfCoverage \
            -o $my_out_prefix \
            -I $my_bam \
            -R $my_reference \
            --omitDepthOutputAtEachBase \
            -L $my_gatk_intervals

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        ls $my_out_prefix* > output_files.txt
        
        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_prefix* ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}

                while IFS= read -r fname || [[ -n "$fname" ]]; do
                    gsutil cp $fname ~{output_dir}
                done < output_files.txt
                sed -i "s;^;~{output_dir};g" output_files.txt

                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        Array [String] output_files = read_lines("output_files.txt")
        String output_files_prefix = out_prefix
    }

    # TODO may need more memory; previous iteration used -Xmx4g
    # Runs >24h - no preemptible
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk:3.8"
        cpu: 2
        memory: "32 GB"
        disks: "local-disk " + disk_size + " HDD"
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task CollectInsertSizeMetrics {
    input {
        File bam
        File bai
        String log_dir
        String task_name
        String out_metrics
        String out_histogram

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String bam_base = basename(bam)
    String bai_base = basename(bai)
    String sa_key_base = basename(sa_key)

    String out_metrics_base = basename(out_metrics)
    String out_histogram_base = basename(out_histogram)

    Int disk_size = ceil(size(bam, "GB") * 1.2 + 20)

    command <<<
        set -eo pipefail

        echo "In CollectInsertSizeMetrics"
        env

        my_bam=~{bam}
        my_bai=~{bai}
        my_sa_key=~{sa_key}

        my_out_metrics=~{out_metrics_base}
        my_out_histogram=~{out_histogram_base}

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
                    #cp ~{bam} ~{bai} /scratch/
                    #rm ~{bam} ~{bai}
                    #my_bam=/scratch/~{bam_base}
                    #my_bai=/scratch/~{bai_base}
                fi
            fi

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_metrics=/scratch/~{out_metrics_base}
            my_out_histogram=/scratch/~{out_histogram_base}
        fi

        time \
        gatk \
            --java-options "-Xmx8g" \
            CollectInsertSizeMetrics \
            -I $my_bam \
            -O $my_out_metrics \
            -H $my_out_histogram \
            --VALIDATION_STRINGENCY LENIENT

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_metrics ~{out_metrics}
                cp $my_out_histogram ~{out_histogram}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_metrics_base} ~{output_dir}
                gsutil cp ~{out_histogram_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_metrics = "~{out_metrics}"
        String output_histogram = "~{out_histogram}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gatk-r:4.1.3.0_3.6.1"
        cpu: 4
        memory: "15 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "4:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task FASTQC {
    input {
        File bam
        Int? num_threads
        String log_dir
        String task_name
        String out_zip
        String out_html

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String bam_base = basename(bam)
    String sa_key_base = basename(sa_key)

    String fastqc_output_dir = if (environment == "LOCAL") then output_dir else "."

    String out_zip_base = basename(out_zip)
    String out_html_base = basename(out_html)

    Int n_threads = if (defined(num_threads)) then select_first([num_threads]) else 2
    Int disk_size = ceil(size(bam, "GB") * 1.2 + 20)

    command <<<
        set -eo pipefail

        echo "In FASTQC"
        env

        my_bam=~{bam}
        my_sa_key=~{sa_key}

        my_out_zip=~{out_zip_base}
        my_out_html=~{out_html_base}

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
                    #cp ~{bam} /scratch/
                    #rm ~{bam}
                    #my_bam=/scratch/~{bam_base}
                fi
            fi

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_zips=/scratch/~{out_zip_base}
            my_out_html=/scratch/~{out_html_base}
        fi
        
        time \
        fastqc \
            --format bam \
            --outdir ~{fastqc_output_dir} \
            --threads ~{n_threads} \
            $my_bam

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                #cp $my_out_zips $my_out_html ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_zip_base} ~{output_dir}
                gsutil cp ~{out_html_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_zip = "~{out_zip}"
        String output_html = "~{out_html}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/fastqc:0.11.8"
        cpu: n_threads
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}
