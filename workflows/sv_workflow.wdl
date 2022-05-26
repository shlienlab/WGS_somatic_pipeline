version 1.0
import "common_tasks.wdl" as CommonTasks

workflow SV {
    input {
        String output_dir
        String patient
        String sample
        File tumor_bam
        File tumor_bai
        File tumor_gatk_depth_of_coverage
        File normal_bam
        File normal_bai
        Array[String] sv_types
        Map[String, String] delly_sv_types
        File reference
        File reference_index
        File reference_dict
        File Delly_header
        File delly_exclude_file
        File centromeres
        File dustmaker
        Map[String, File] panel_of_normals_beds
        File cosmic_cancer_genes_rda
        File ucsc_info_tar
        File circos_config_karyotype
        File circos_config_ideogram_conf
        File circos_config_image_conf
        File circos_config_colors_fonts_patterns_conf
        File circos_config_housekeeping_conf
        File circos_config_ticks_conf

        String repo_dir
        File sa_key
        String environment
    }

    String normal = "normal"
    String tumor = "tumor"

    # Directories
    String log_dir = "~{output_dir}/wdl_logs+run_info"
    String progress_dir = "~{log_dir}/progress"

    # Task Names
    String Delly_normal_task_name_base = "~{patient}.~{sample}.delly.~{normal}"
    String Delly_tumor_task_name_base = "~{patient}.~{sample}.delly.~{tumor}"
    String VCF2Tab_normal_task_name_base = "~{patient}.~{sample}.vcf2tab.~{normal}"
    String VCF2Tab_tumor_task_name_base = "~{patient}.~{sample}.vcf2tab.~{tumor}"
    String Filtering_task_name_base = "~{patient}.~{sample}.filter1"
    String FilterRearrangements_task_name_base = "~{patient}.~{sample}.filter2"
    String PairToPair_task_name_base = "~{patient}.~{sample}.p2p"
    String PONFiltering_task_name_base = "~{patient}.~{sample}.pon"
    scatter (sv_type in sv_types) {
        String Delly_normal_task_names = "~{Delly_normal_task_name_base}.~{sv_type}"
        String Delly_tumor_task_names = "~{Delly_tumor_task_name_base}.~{sv_type}"
        String VCF2Tab_normal_task_names = "~{VCF2Tab_normal_task_name_base}.~{sv_type}"
        String VCF2Tab_tumor_task_names = "~{VCF2Tab_tumor_task_name_base}.~{sv_type}"
        String Filtering_task_names = "~{Filtering_task_name_base}.~{sv_type}"
        String FilterRearrangements_task_names = "~{FilterRearrangements_task_name_base}.~{sv_type}"
        String PairToPair_task_names = "~{PairToPair_task_name_base}.~{sv_type}"
        String PONFiltering_task_names = "~{PONFiltering_task_name_base}.~{sv_type}"
    }
    String MergeMutationTabs_task_name = "~{patient}.~{sample}.merge"
    String CancerGeneLists_task_name = "~{patient}.~{sample}.lists"
    String CircosConfigs_task_name = "~{patient}.~{sample}.config"
    String CircosPlots_task_name = "~{patient}.~{sample}.plot"
    String FilterGraph_task_name = "~{patient}.~{sample}.FilterGraph"
    Array[String] task_names = flatten([Delly_normal_task_names, Delly_tumor_task_names, VCF2Tab_normal_task_names, VCF2Tab_tumor_task_names, Filtering_task_names, FilterRearrangements_task_names, PairToPair_task_names, PONFiltering_task_names, [MergeMutationTabs_task_name, CancerGeneLists_task_name, CircosConfigs_task_name, CircosPlots_task_name, FilterGraph_task_name]])

    if (environment=="LOCAL") {
        call CommonTasks.FindWorkDir as SVFindWorkDir {
            input:
                sa_key = sa_key,
                patient=patient,
                sample=sample,
                my_workflow='SV',
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

    scatter (sv_type in sv_types) {
        String Delly_normal_task_name = "~{Delly_normal_task_name_base}.~{sv_type}"
        String Delly_normal_output_vcf = "~{output_dir}/~{normal}_~{sample}.~{sv_type}.vcf"
        if (CheckProgress.task_completion[Delly_normal_task_name] == "false") {
            call Delly as Delly_normal {
                input:
                    input_bam = normal_bam,
                    input_bai = normal_bai,
                    delly_sv_type = delly_sv_types[sv_type],
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    Delly_header = Delly_header,
                    exclude_file = delly_exclude_file,
                    log_dir = log_dir,
                    task_name = Delly_normal_task_name,
                    out_vcf = Delly_normal_output_vcf,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String Delly_tumor_task_name = "~{Delly_tumor_task_name_base}.~{sv_type}"
        String Delly_tumor_output_vcf = "~{output_dir}/~{tumor}_~{sample}.~{sv_type}.vcf"
        if (CheckProgress.task_completion[Delly_tumor_task_name] == "false") {
            call Delly as Delly_tumor {
                input:
                    input_bam = tumor_bam,
                    input_bai = tumor_bai,
                    delly_sv_type = delly_sv_types[sv_type],
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    Delly_header = Delly_header,
                    exclude_file = delly_exclude_file,
                    log_dir = log_dir,
                    task_name = Delly_tumor_task_name,
                    out_vcf = Delly_tumor_output_vcf,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String VCF2Tab_normal_task_name = "~{VCF2Tab_normal_task_name_base}.~{sv_type}"
        String VCF2Tab_normal_output_tab = "~{output_dir}/~{normal}_~{sample}.~{sv_type}.tab"
        String VCF2Tab_normal_output_filter_info = "~{output_dir}/~{sample}.normal.~{sv_type}.filter_info.txt"
        if (CheckProgress.task_completion[VCF2Tab_normal_task_name] == "false") {
            call VCF2Tab as VCF2Tab_normal {
                input:
                    vcf = select_first([Delly_normal.output_vcf, Delly_normal_output_vcf]),
                    bam = normal_bam,
                    bai = normal_bai,
                    sample_type = normal,
                    ucsc_info_tar = ucsc_info_tar,
                    log_dir = log_dir,
                    task_name = VCF2Tab_normal_task_name,
                    out_tab = VCF2Tab_normal_output_tab,
                    filter_info = VCF2Tab_normal_output_filter_info,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String VCF2Tab_tumor_task_name = "~{VCF2Tab_tumor_task_name_base}.~{sv_type}"
        String VCF2Tab_tumor_output_tab = "~{output_dir}/~{tumor}_~{sample}.~{sv_type}.tab"
        String VCF2Tab_tumor_output_filter_info = "~{output_dir}/~{sample}.tumor.~{sv_type}.filter_info.txt"
        if (CheckProgress.task_completion[VCF2Tab_tumor_task_name] == "false") {
            call VCF2Tab as VCF2Tab_tumor {
                input:
                    vcf = select_first([Delly_tumor.output_vcf, Delly_tumor_output_vcf]),
                    bam = tumor_bam,
                    bai = tumor_bai,
                    sv_type = sv_type,
                    sample_type = tumor,
                    ucsc_info_tar = ucsc_info_tar,
                    log_dir = log_dir,
                    task_name = VCF2Tab_tumor_task_name,
                    out_tab = VCF2Tab_tumor_output_tab,
                    filter_info = VCF2Tab_tumor_output_filter_info,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String Filtering_task_name = "~{Filtering_task_name_base}.~{sv_type}"
        String Filtering_output_tab = "~{output_dir}/~{sample}.unfiltered.~{sv_type}.tab"
        if (CheckProgress.task_completion[Filtering_task_name] == "false") {
            call Filtering {
                input:
                    tumor_tab = select_first([VCF2Tab_tumor.output_tab, VCF2Tab_tumor_output_tab]),
                    tumor_bam = tumor_bam,
                    tumor_bai = tumor_bai,
                    gatk_path = tumor_gatk_depth_of_coverage,
                    normal_tab = select_first([VCF2Tab_normal.output_tab, VCF2Tab_normal_output_tab]),
                    normal_bam = normal_bam,
                    normal_bai = normal_bai,
                    reference = reference,
                    reference_index = reference_index,
                    reference_dict = reference_dict,
                    centromeres = centromeres,
                    dustmaker = dustmaker,
                    log_dir = log_dir,
                    task_name = Filtering_task_name,
                    out_tab = Filtering_output_tab,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String FilterRearrangements_task_name = "~{FilterRearrangements_task_name_base}.~{sv_type}"
        String FilterRearrangements_output_tab = "~{output_dir}/~{sample}.almost_filtered.~{sv_type}.tab"
        String FilterRearrangements_output_bed = "~{output_dir}/~{sample}.Ftumours.~{sv_type}.bed"
        String FilterRearrangements_output_filter_info = "~{VCF2Tab_tumor_output_filter_info}.filter_rearrangements"
        Array [String] FilterRearrangements_files_to_delete = [VCF2Tab_tumor_output_filter_info]
        if (CheckProgress.task_completion[FilterRearrangements_task_name] == "false") {
            call FilterRearrangements {
                input:
                    sample = sample,
                    tab_file = select_first([Filtering.output_tab, Filtering_output_tab]),
                    filter_info = select_first([VCF2Tab_tumor.output_filter_info, VCF2Tab_tumor_output_filter_info]),
                    sv_type = sv_type,
                    log_dir = log_dir,
                    task_name = FilterRearrangements_task_name,
                    out_tab = FilterRearrangements_output_tab,
                    out_bed = FilterRearrangements_output_bed,
                    out_filter_info = FilterRearrangements_output_filter_info,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = FilterRearrangements_files_to_delete
            }
        }

        String PairToPair_task_name = "~{PairToPair_task_name_base}.~{sv_type}"
        String PairToPair_output_bed = "~{output_dir}/~{sample}.all_both.~{sv_type}.bed"
        if (CheckProgress.task_completion[PairToPair_task_name] == "false") {
            call PairToPair {
                input:
                    Ftumors_bed = select_first([FilterRearrangements.output_bed, FilterRearrangements_output_bed]),
                    normals_bed = panel_of_normals_beds[sv_type],
                    log_dir = log_dir,
                    task_name = PairToPair_task_name,
                    out_bed = PairToPair_output_bed,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String PONFiltering_task_name = "~{PONFiltering_task_name_base}.~{sv_type}"
        String PONFiltering_output_tab = "~{output_dir}/~{sample}.filtered.~{sv_type}.tab"
        String PONFiltering_output_filter_info = "~{FilterRearrangements_output_filter_info}.PONFiltering"
        Array [String] PONFiltering_files_to_delete = [FilterRearrangements_output_filter_info]
        if (CheckProgress.task_completion[PONFiltering_task_name] == "false") {
            call PONFiltering {
                input:
                    PairToPair_bed = select_first([PairToPair.output_bed, PairToPair_output_bed]),
                    almost_filtered_tab = select_first([FilterRearrangements.output_tab, FilterRearrangements_output_tab]),
                    filter_info = select_first([FilterRearrangements.output_filter_info, FilterRearrangements_output_filter_info]),
                    sv_type = sv_type,
                    log_dir = log_dir,
                    task_name = PONFiltering_task_name,
                    out_tab = PONFiltering_output_tab,
                    out_filter_info = PONFiltering_output_filter_info,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = PONFiltering_files_to_delete
            }
        }
    }

    String MergeMutationTabs_output_tab = "~{output_dir}/~{sample}.final_rearrangements.tab"
    if (CheckProgress.task_completion[MergeMutationTabs_task_name] == "false") {
        Array[File] MergeMutationTabs_inputs = if (length(select_all(PONFiltering.output_tab)) == length(sv_types)) then select_all(PONFiltering.output_tab) else PONFiltering_output_tab
        call MergeMutationTabs {
            input:
                final_sv_tabs = MergeMutationTabs_inputs,
                log_dir = log_dir,
                task_name = MergeMutationTabs_task_name,
                out_tab = MergeMutationTabs_output_tab,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir
        }
    }

    if (select_first([MergeMutationTabs.no_variants, "true"]) == "false") {
        String CancerGeneLists_output_tab = "~{output_dir}/~{sample}.cancer_rearrangements.tab"
        String CancerGeneLists_output_CF_tab = "~{output_dir}/~{sample}.CFformat_rearrangements.tab"
        if (CheckProgress.task_completion[CancerGeneLists_task_name] == "false") {
            call CancerGeneLists {
                input:
                    merged_tab = select_first([MergeMutationTabs.output_tab, MergeMutationTabs_output_tab]),
                    cosmic_rda = cosmic_cancer_genes_rda,
                    log_dir = log_dir,
                    task_name = CancerGeneLists_task_name,
                    out_tab = CancerGeneLists_output_tab,
                    out_cf_tab = CancerGeneLists_output_CF_tab,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String CircosConfigs_output_circos_all_file = "~{output_dir}/~{sample}.circos_all.txt"
        String CircosConfigs_output_config_file = "~{output_dir}/~{sample}.circos_configuration.tab"
        if (CheckProgress.task_completion[CircosConfigs_task_name] == "false") {
            call CircosConfigs {
                input:
                    merged_tab = select_first([MergeMutationTabs.output_tab, MergeMutationTabs_output_tab]),
                    sample = sample,
                    circos_config_karyotype = circos_config_karyotype,
                    circos_config_ideogram_conf = circos_config_ideogram_conf,
                    circos_config_image_conf = circos_config_image_conf,
                    circos_config_colors_fonts_patterns_conf = circos_config_colors_fonts_patterns_conf,
                    circos_config_housekeeping_conf = circos_config_housekeeping_conf,
                    circos_config_ticks_conf = circos_config_ticks_conf,
                    log_dir = log_dir,
                    task_name = CircosConfigs_task_name,
                    out_circos_all_file = CircosConfigs_output_circos_all_file,
                    out_config_file = CircosConfigs_output_config_file,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir
            }
        }

        String CircosPlots_output_prefix = "~{output_dir}/~{sample}.circos"
        String CircosPlots_output_png = "~{CircosPlots_output_prefix}.png"
        String CircosPlots_output_svg = "~{CircosPlots_output_prefix}.svg"
        Array [String] CircosPlots_files_to_delete = [CircosConfigs_output_circos_all_file, CircosConfigs_output_config_file]
        if (CheckProgress.task_completion[CircosPlots_task_name] == "false") {
            call CircosPlots {
                input:
                    config = select_first([CircosConfigs.output_config_file, CircosConfigs_output_config_file]),
                    config_all = select_first([CircosConfigs.output_circos_all_file, CircosConfigs_output_circos_all_file]),
                    circos_config_karyotype = circos_config_karyotype,
                    circos_config_ideogram_conf = circos_config_ideogram_conf,
                    circos_config_image_conf = circos_config_image_conf,
                    circos_config_colors_fonts_patterns_conf = circos_config_colors_fonts_patterns_conf,
                    circos_config_housekeeping_conf = circos_config_housekeeping_conf,
                    circos_config_ticks_conf = circos_config_ticks_conf,
                    log_dir = log_dir,
                    task_name = CircosPlots_task_name,
                    out_prefix = CircosPlots_output_prefix,
                    out_png = CircosPlots_output_png,
                    out_svg = CircosPlots_output_svg,
                    environment = environment,
                    sa_key = sa_key,
                    progress_dir = progress_dir,
                    output_dir = output_dir,
                    files_to_delete = CircosPlots_files_to_delete
            }
        }
    }
    
    String FilterGraph_output_file = "~{output_dir}/~{sample}.filtering.png"
    Array [String] FilterGraph_files_to_delete = PONFiltering_output_filter_info
    if (CheckProgress.task_completion[FilterGraph_task_name] == "false") {
        Array[File] FilterGraph_inputs = if (length(select_all(PONFiltering.output_filter_info)) == length(sv_types)) then select_all(PONFiltering.output_filter_info) else PONFiltering_output_filter_info
        call FilterGraph {
            input:
                sample = sample,
                filter_info = FilterGraph_inputs,
                log_dir = log_dir,
                task_name = FilterGraph_task_name,
                out_file = FilterGraph_output_file,
                environment = environment,
                sa_key = sa_key,
                progress_dir = progress_dir,
                output_dir = output_dir,
                files_to_delete = FilterGraph_files_to_delete
        }
    }

    output {
        Array[String] normal_vcfs = if (length(select_all(Delly_normal.output_vcf)) == length(sv_types)) then select_all(Delly_normal.output_vcf) else Delly_normal_output_vcf
        Array[String] tumor_vcfs = if (length(select_all(Delly_tumor.output_vcf)) == length(sv_types)) then select_all(Delly_tumor.output_vcf) else Delly_tumor_output_vcf
        Array[String] normal_tabs = if (length(select_all(VCF2Tab_normal.output_tab)) == length(sv_types)) then select_all(VCF2Tab_normal.output_tab) else VCF2Tab_normal_output_tab
        Array[String] tumor_tabs = if (length(select_all(VCF2Tab_tumor.output_tab)) == length(sv_types)) then select_all(VCF2Tab_tumor.output_tab) else VCF2Tab_tumor_output_tab
        Array[String] unfiltered_tabs = if (length(select_all(Filtering.output_tab)) == length(sv_types)) then select_all(Filtering.output_tab) else Filtering_output_tab
        Array[String] almost_filtered_tabs = if (length(select_all(FilterRearrangements.output_tab)) == length(sv_types)) then select_all(FilterRearrangements.output_tab) else FilterRearrangements_output_tab
        Array[String] Ftumors_beds = if (length(select_all(FilterRearrangements.output_bed)) == length(sv_types)) then select_all(FilterRearrangements.output_bed) else FilterRearrangements_output_bed
        Array[String] PairToPair_beds = if (length(select_all(PairToPair.output_bed)) == length(sv_types)) then select_all(PairToPair.output_bed) else PairToPair_output_bed
        Array[String] filtered_tabs = if (length(select_all(PONFiltering.output_tab)) == length(sv_types)) then select_all(PONFiltering.output_tab) else PONFiltering_output_tab
        String final_rearrangements_tab = select_first([MergeMutationTabs.output_tab, MergeMutationTabs_output_tab])
        # String cancer_rearrangements_tab = select_first([CancerGeneLists.output_tab, CancerGeneLists_output_tab])
        # String CF_format_rearrangements_tab = select_first([CancerGeneLists.output_CF_tab, CancerGeneLists_output_CF_tab])
        # String circos_png = select_first([CircosPlots.output_png, CircosPlots_output_png])
        # String circos_svg = select_first([CircosPlots.output_svg, CircosPlots_output_svg])
        String filtering_graph = select_first([FilterGraph.output_file, FilterGraph_output_file])
    }
}

task Delly {
    input {
        File input_bam
        File input_bai
        String delly_sv_type
        File reference
        File reference_index
        File reference_dict
        File exclude_file
        File Delly_header
        String log_dir
        String task_name
        String out_vcf

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
    String exclude_file_base = basename(exclude_file)
    String Delly_header_base = basename(Delly_header)
    String sa_key_base = basename(sa_key)

    String out_vcf_base = basename(out_vcf)

    Int disk_size = ceil(size(input_bam, "GB") + size(reference, "GB") * 2 + 50)

    command <<<
        set -eo pipefail

        env

        my_input_bam=~{input_bam}
        my_input_bai=~{input_bai}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_exclude_file=~{exclude_file}
        my_Delly_header=~{Delly_header}
        my_sa_key=~{sa_key}

        my_out_vcf=~{out_vcf_base}

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

            cp ~{reference} ~{reference_index} ~{reference_dict} /scratch/
            rm ~{reference} ~{reference_index} ~{reference_dict}
            my_reference=/scratch/~{reference_base}
            my_reference_index=/scratch/~{reference_index_base}
            my_reference_dict=/scratch/~{reference_dict_base}

            cp ~{exclude_file} /scratch/
            my_exclude_file=/scratch/~{exclude_file_base}
            rm ~{exclude_file}

            cp ~{Delly_header} /scratch/
            my_Delly_header=/scratch/~{Delly_header_base}
            rm ~{Delly_header}

            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            my_out_vcf=/scratch/~{out_vcf_base}
        fi

        time \
        delly \
            -t ~{delly_sv_type} \
            -x $my_exclude_file \
            -g $my_reference \
            -o $my_out_vcf \
            $my_input_bam &&

        if [ ! -f $my_out_vcf ]; then
            head -c -1 $my_Delly_header > $my_out_vcf
            echo -e "\t$(basename $my_input_bam '.bam')" >> $my_out_vcf
        fi

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_vcf ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_vcf_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_vcf = "~{out_vcf}"
    }

    # may need more memory (13GB), disk size(*20??)
    # TODO if disk size not an issue, reduce prod delly disk size
    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/delly:0.7.1S"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        walltime: "200:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task VCF2Tab {
    input {
        File vcf
        File bam
        File bai
        String? sv_type
        String sample_type
        File ucsc_info_tar
        String log_dir
        String task_name
        String out_tab
        String filter_info

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String bam_base = basename(bam)
    String bai_base = basename(bai)
    String vcf_base = basename(vcf)
    String ucsc_info_tar_base = basename(ucsc_info_tar)
    String sa_key_base = basename(sa_key)

    String out_tab_base = if (environment == "LOCAL") then out_tab else basename(out_tab)
    String filter_info_base = if (environment == "LOCAL") then filter_info else basename(filter_info)

    String low_qual_flag = if (sample_type == "normal") then "--lowqual True" else ""
    Int disk_size = ceil(size(vcf, "GB") * 2 + size(bam, "GB") + 50)

    command <<<
        set -eo pipefail

        env

        my_bam=~{bam}
        my_bai=~{bai}
        my_vcf=~{vcf}
        my_ucsc_info_tar=~{ucsc_info_tar}
        my_sa_key=~{sa_key}

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

            cp ~{vcf} ~{ucsc_info_tar} ~{sa_key} /scratch/
            rm ~{vcf} ~{ucsc_info_tar} ~{sa_key}
            my_vcf=/scratch/~{vcf_base}
            my_ucsc_info_tar=/scratch/~{ucsc_info_tar_base}
            my_sa_key=/scratch/~{sa_key_base}

            cd /scratch
        fi

        tar -zxvf $my_ucsc_info_tar

        time \
        python3 -B /opt/scripts/vcf2tab.py \
            ~{low_qual_flag} \
            --vcf $my_vcf \
            --bam $my_bam \
            --ucsc_info_path . \
            --output_file ~{out_tab_base} \
            --filter_info_file ~{filter_info_base} \
            ~{'--sv_type ' + sv_type}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_tab_base} ~{output_dir}
                gsutil cp ~{filter_info_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    >>>

    output {
        String output_tab = "~{out_tab}"
        String output_filter_info = "~{filter_info}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-python-libraries:3.5"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk " + disk_size + " HDD"
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

# TODO not certain if gatk_path is necessary here as it is in the ssm filtering
task Filtering {
    input {
        File tumor_tab
        File tumor_bam
        File tumor_bai
        File gatk_path
        File normal_tab
        File normal_bam
        File normal_bai
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

    String tumor_tab_base = basename(tumor_tab)
    String tumor_bam_base = basename(tumor_bam)
    String tumor_bai_base = basename(tumor_bai)
    String gatk_path_base = basename(gatk_path)
    String normal_tab_base = basename(normal_tab)
    String normal_bam_base = basename(normal_bam)
    String normal_bai_base = basename(normal_bai)
    String reference_base = basename(reference)
    String reference_index_base = basename(reference_index)
    String reference_dict_base = basename(reference_dict)
    String centromeres_base = basename(centromeres)
    String dustmaker_base = basename(dustmaker)
    String sa_key_base = basename(sa_key)

    String out_tab_base = if (environment == "LOCAL") then out_tab else basename(out_tab)

    Int disk_size = ceil((size(tumor_tab, "GB") + size(tumor_bam, "GB") + size(normal_bam, "GB") + size(reference, "GB")) * 1.2 + 20)

    command <<<
        set -eo pipefail

        env

        my_tumor_tab=~{tumor_tab}
        my_tumor_bam=~{tumor_bam}
        my_tumor_bai=~{tumor_bai}
        my_gatk_path=~{gatk_path}
        my_normal_tab=~{normal_tab}
        my_normal_bam=~{normal_bam}
        my_normal_bai=~{normal_bai}
        my_reference=~{reference}
        my_reference_index=~{reference_index}
        my_reference_dict=~{reference_dict}
        my_centromeres=~{centromeres}
        my_dustmaker=~{dustmaker}
        my_sa_key=~{sa_key}

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

            cp ~{tumor_tab} ~{normal_tab} ~{gatk_path} /scratch/
            rm ~{tumor_tab} ~{normal_tab} ~{gatk_path}
            my_tumor_tab=/scratch/~{tumor_tab_base}
            my_normal_tab=/scratch/~{normal_tab_base}
            my_gatk_path=/scratch/~{gatk_path_base}

            cp ~{reference} ~{reference_index} ~{reference_dict} /scratch/
            rm ~{reference} ~{reference_index} ~{reference_dict}
            my_reference=/scratch/~{reference_base}
            my_reference_index=/scratch/~{reference_index_base}
            my_reference_dict=/scratch/~{reference_dict_base}

            cp ~{centromeres} ~{dustmaker} ~{sa_key} /scratch/
            rm ~{centromeres} ~{dustmaker} ~{sa_key}
            my_centromeres=/scratch/~{centromeres_base}
            my_dustmaker=/scratch/~{dustmaker_base}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        python3 -B /opt/filterpipeline/runFilters.py \
            --tumor_tab $my_tumor_tab \
            --tumor_bam $my_tumor_bam \
            --normal_tab $my_normal_tab \
            --normal_bam $my_normal_bam \
            --config /opt/filterpipeline/config.yaml \
            --reference $my_reference \
            --centromeres $my_centromeres \
            --dustmaker $my_dustmaker \
            --output_tab ~{out_tab_base} \
            --cFilter --sFilter --bedFilter --baffler \
            --gatk_path $my_gatk_path

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
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
        walltime: "72:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

# TODO may need to change this to be based on R 3.1.1
task FilterRearrangements {
    input {
        String sample
        File tab_file
        File filter_info
        String sv_type
        String log_dir
        String task_name
        String out_tab
        String out_bed
        String out_filter_info

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String tab_file_base = basename(tab_file)
    String sa_key_base = basename(sa_key)

    String filter_info_base = basename(filter_info)

    String out_tab_base = basename(out_tab)
    String out_bed_base = basename(out_bed)
    String out_filter_info_base =basename(out_filter_info)

    Int disk_size = ceil(size(tab_file, "GB")*2 + 20)

    command {
        set -eo pipefail

        env

        my_tab_file=~{tab_file}
        my_sa_key=~{sa_key}

        my_out_tab=~{out_tab_base}
        my_out_bed=~{out_bed_base}
        my_out_filter_info=~{out_filter_info_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{tab_file} ~{sa_key} /scratch/
            rm ~{tab_file} ~{sa_key}
            my_tab_file=/scratch/~{tab_file_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_ucsc_info_path=/scratch/

            my_out_tab=/scratch/~{out_tab_base}
            my_out_bed=/scratch/~{out_bed_base}
            my_out_filter_info=/scratch/~{out_filter_info_base}
        fi

        touch $my_out_bed

        cp ~{filter_info} .

        time \
        Rscript /opt/scripts/filtering_rearrangements.R \
            $my_tab_file \
            ~{sv_type} \
            ~{sample} \
            $my_out_tab \
            $my_out_bed \
            $my_out_filter_info

        cp ~{filter_info_base} $my_out_filter_info

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_tab ~{output_dir}
                cp $my_out_bed ~{output_dir}
                cp $my_out_filter_info ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_tab_base} ~{output_dir}
                gsutil cp ~{out_bed_base} ~{output_dir}
                gsutil cp ~{out_filter_info_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_tab = "~{out_tab}"
        String output_bed = "~{out_bed}"
        String output_filter_info = "~{out_filter_info}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task PairToPair {
    input {
        File Ftumors_bed
        File normals_bed
        String log_dir
        String task_name
        String out_bed

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String Ftumors_bed_base = basename(Ftumors_bed)
    String normals_bed_base = basename(normals_bed)
    String sa_key_base = basename(sa_key)

    String out_bed_base = basename(out_bed)

    Int disk_size = ceil((size(Ftumors_bed, "GB") + size(normals_bed, "GB")) * 2 + 20)

    command {
        set -eo pipefail

        env

        my_Ftumors_bed=~{Ftumors_bed}
        my_normals_bed=~{normals_bed}
        my_sa_key=~{sa_key}

        my_out_bed=~{out_bed_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{Ftumors_bed} ~{normals_bed} ~{sa_key} /scratch/
            rm ~{Ftumors_bed} ~{normals_bed} ~{sa_key}
            my_Ftumors_bed=/scratch/~{Ftumors_bed_base}
            my_normals_bed=/scratch/~{normals_bed_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_bed=/scratch/~{out_bed_base}
        fi

        time \
        pairToPair \
            -type both \
            -a $my_Ftumors_bed \
            -b $my_normals_bed \
            > $my_out_bed

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_bed ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_bed_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_bed = "~{out_bed}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/bedtools:2.27.1"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

# TODO might have to actually make this R 3.1.1
task PONFiltering {
    input {
        File PairToPair_bed
        File almost_filtered_tab
        File filter_info
        String sv_type
        String log_dir
        String task_name
        String out_tab
        String out_filter_info

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String PairToPair_bed_base = basename(PairToPair_bed)
    String almost_filtered_tab_base = basename(almost_filtered_tab)
    String sa_key_base = basename(sa_key)

    String filter_info_base = basename(filter_info)

    String out_tab_base = basename(out_tab)
    String out_filter_info_base = basename(out_filter_info)

    Int disk_size = ceil((size(PairToPair_bed, "GB") + size(almost_filtered_tab, "GB")) * 2 + 20)

    command {
        set -eo pipefail

        env

        my_PairToPair_bed=~{PairToPair_bed}
        my_almost_filtered_tab=~{almost_filtered_tab}
        my_sa_key=~{sa_key}

        my_out_tab=~{out_tab_base}
        my_out_filter_info=~{out_filter_info_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{PairToPair_bed} ~{almost_filtered_tab} ~{sa_key} /scratch/
            rm ~{PairToPair_bed} ~{almost_filtered_tab} ~{sa_key}
            my_PairToPair_bed=/scratch/~{PairToPair_bed_base}
            my_almost_filtered_tab=/scratch/~{almost_filtered_tab_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_tab=/scratch/~{out_tab_base}
            my_out_filter_info=/scratch/~{out_filter_info_base}
        fi

        cp ~{filter_info} .

        time \
        Rscript /opt/scripts/panel_of_normals.R \
            $my_PairToPair_bed \
            $my_almost_filtered_tab \
            $my_out_tab \
            ~{filter_info_base} \
            ~{sv_type}

        cp ~{filter_info_base} $my_out_filter_info

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_tab ~{output_dir}
                cp $my_out_filter_info ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_tab_base} ~{output_dir}
                gsutil cp ~{out_filter_info_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_tab = "~{out_tab}"
        String output_filter_info = "~{out_filter_info}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        walltime: "48:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task MergeMutationTabs {
    input {
        Array[File] final_sv_tabs
        String log_dir
        String task_name
        String out_tab

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String sa_key_base = basename(sa_key)

    String out_tab_base = if (environment == "LOCAL") then out_tab else basename(out_tab)

    Int disk_size = ceil(size(final_sv_tabs[0], "GB") * length(final_sv_tabs) * 2 + 20)

    command {
        set -eo pipefail

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            rm ~{sa_key}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        Rscript /opt/scripts/merge_mutation_tabs.R \
            ~{out_tab_base} \
            ~{sep=" " final_sv_tabs} \
        && \
        if [ `wc -l < ~{out_tab_base}` -eq 1 ]; then echo "true"; else echo "false"; fi

        # Cleanup
        EXIT_STATUS=$?
        # echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
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
    }

    output {
        String no_variants = read_string(stdout())
        String output_tab = "~{out_tab}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        walltime: "24:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task CancerGeneLists {
    input {
        File merged_tab
        File cosmic_rda
        String log_dir
        String task_name
        String out_tab
        String out_cf_tab

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String merged_tab_base = basename(merged_tab)
    String cosmic_rda_base = basename(cosmic_rda)
    String sa_key_base = basename(sa_key)

    String out_tab_base = basename(out_tab)
    String out_cf_tab_base = basename(out_cf_tab)

    Int disk_size = ceil((size(merged_tab, "GB") + size(cosmic_rda, "GB")) * 2 + 20)

    command {
        set -eo pipefail

        env

        my_merged_tab=~{merged_tab}
        my_cosmic_rda=~{cosmic_rda}
        my_sa_key=~{sa_key}

        my_out_tab=~{out_tab_base}
        my_out_cf_tab=~{out_cf_tab_base}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{merged_tab} ~{cosmic_rda} ~{sa_key} /scratch/
            rm ~{merged_tab} ~{cosmic_rda} ~{sa_key}
            my_merged_tab=/scratch/~{merged_tab_base}
            my_cosmic_rda=/scratch/~{cosmic_rda_base}
            my_sa_key=/scratch/~{sa_key_base}

            my_out_tab=/scratch/~{out_tab_base}
            my_out_cf_tab=/scratch/~{out_cf_tab_base}
        fi

        time \
        Rscript /opt/scripts/create_cancer_gene_lists.R \
            $my_merged_tab \
            $my_cosmic_rda \
            $my_out_tab \
            $my_out_cf_tab

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                cp $my_out_tab ~{output_dir}
                cp $my_out_cf_tab ~{output_dir}
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_tab_base} ~{output_dir}
                gsutil cp ~{out_cf_tab_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_tab = "~{out_tab}"
        String output_CF_tab = "~{out_cf_tab}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        walltime: "24:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task CircosConfigs {
    input {
        File merged_tab
        String sample
        File circos_config_karyotype
        File circos_config_ideogram_conf
        File circos_config_image_conf
        File circos_config_colors_fonts_patterns_conf
        File circos_config_housekeeping_conf
        File circos_config_ticks_conf
        String log_dir
        String task_name
        String out_circos_all_file
        String out_config_file

        String environment
        String output_dir
        String progress_dir
        File sa_key
    }

    String sa_key_base = basename(sa_key)

    String out_circos_all_file_base = if (environment == "LOCAL") then out_circos_all_file else basename(out_circos_all_file)
    String out_config_file_base = if (environment == "LOCAL") then out_config_file else basename(out_config_file)

    Int disk_size = ceil(size(merged_tab, "GB") * 3 + 20)

    command {
        set -eo pipefail

        env

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            rm ~{sa_key}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        mv ~{circos_config_karyotype} ~{circos_config_ideogram_conf} ~{circos_config_image_conf} ~{circos_config_colors_fonts_patterns_conf} ~{circos_config_housekeeping_conf} ~{circos_config_ticks_conf} .

        time \
        Rscript /opt/scripts/create_configs_circos.R \
            ~{merged_tab} \
            ~{sample} \
            ~{out_circos_all_file_base} \
            ~{out_config_file_base} \
            ~{basename(circos_config_karyotype)} \
            ~{basename(circos_config_ideogram_conf)} \
            ~{basename(circos_config_image_conf)} \
            ~{basename(circos_config_colors_fonts_patterns_conf)} \
            ~{basename(circos_config_housekeeping_conf)} \
            ~{basename(circos_config_ticks_conf)}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_circos_all_file_base} ~{output_dir}
                gsutil cp ~{out_config_file_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_config_file = "~{out_config_file}"
        String output_circos_all_file = "~{out_circos_all_file}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/shlienlab-libraries:3.4.0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "24:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task CircosPlots {
    input {
        File config
        File config_all
        File circos_config_karyotype
        File circos_config_ideogram_conf
        File circos_config_image_conf
        File circos_config_colors_fonts_patterns_conf
        File circos_config_housekeeping_conf
        File circos_config_ticks_conf
        String log_dir
        String task_name
        String out_prefix
        String out_png
        String out_svg

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String sa_key_base = basename(sa_key)

    String out_prefix_base = basename(out_prefix)

    String out_png_base = basename(out_png)
    String out_svg_base = basename(out_svg)

    String config_base = basename(config)
    Int disk_size = ceil(size(config, "GB") * 2 + 20)

    command {
        set -eo pipefail

        env

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            rm ~{sa_key}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        mv ~{circos_config_karyotype} ~{circos_config_ideogram_conf} ~{circos_config_image_conf} ~{circos_config_colors_fonts_patterns_conf} ~{circos_config_housekeeping_conf} ~{circos_config_ticks_conf} .
        mv ~{config} ~{config_all} .

        if [ -f ~{config_base} ]; then
            time \
            circos \
                -conf ~{config_base} \
                -outputfile ~{out_prefix_base}
        fi

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                cp ~{out_png_base} ~{output_dir}
                cp ~{out_svg_base} ~{output_dir}
                rm -f ~{config} ~{config_all}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_png_base} ~{output_dir}
                gsutil cp ~{out_svg_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
            fi
        else
            exit 1
        fi
    }

    output {
        String output_png = "~{out_png}"
        String output_svg = "~{out_svg}"
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/private/profyle/circos:0.69"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk " + disk_size + " HDD"
        walltime: "24:00:00"
        log_dir: log_dir
        task_name: task_name
        output_dir: output_dir
    }
}

task FilterGraph {
    input {
        String sample
        Array[File] filter_info
        String log_dir
        String task_name
        String out_file

        String environment
        String output_dir
        String progress_dir
        File sa_key

        Array [String] files_to_delete
    }

    String sa_key_base = basename(sa_key)

    String out_file_base = if (environment == "LOCAL") then out_file else basename(out_file)

    Int disk_size = ceil(size(filter_info[0], "GB") * length(filter_info) * 2 + 20)

    command {
        set -eo pipefail

        env

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            rm ~{sa_key}
            my_sa_key=/scratch/~{sa_key_base}
        fi

        time \
        Rscript /opt/scripts/sv_filtering_graph.R \
            ~{sample} \
            ~{out_file_base} \
            ~{sep=' ' filter_info}

        # Cleanup
        EXIT_STATUS=$?
        echo EXIT STATUS $EXIT_STATUS

        if [ $EXIT_STATUS -eq 0 ]; then
            if [ "~{environment}" = "LOCAL" ]; then
                printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
                rm ~{sep=' ' files_to_delete}
            elif [ "~{environment}" = "CLOUD" ]; then
                gcloud auth activate-service-account --key-file ~{sa_key}
                printf "~{task_name}\ttrue\n" > ./~{task_name}
                gsutil cp ~{out_file_base} ~{output_dir}
                gsutil cp ./~{task_name} ~{progress_dir}
                gsutil rm ~{sep=' ' files_to_delete}
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
        memory: "1 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 4
        walltime: "0:30:00"
        task_name: task_name
        log_dir: log_dir
        output_dir: output_dir
    }
}
