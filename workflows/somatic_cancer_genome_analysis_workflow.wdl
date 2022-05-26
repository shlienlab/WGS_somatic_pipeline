version 1.0

import "align_workflow.wdl" as AlignNamespace
import "gridss_workflow.wdl" as GRIDSSNamespace
import "ssm_workflow.wdl" as SSMNamespace
import "sv_workflow.wdl" as SVNamespace
import "cnv_workflow.wdl" as CNVNamespace
import "default_values.wdl" as DefaultValues

workflow SomaticCancerGenomeAnalysis {
    input {
        Boolean run_GRIDSS
        String random_number
        String environment

        # Output directories
        Map[String, String]? align_output_dirs
        Map[String, Map[String, String]] gridss_output_dirs
        Map[String, Map[String, String]]? ssm_output_dirs # Map[Tumor, Map[Normal, output_dir]]
        Map[String, Map[String, String]]? sv_output_dirs
        Map[String, Map[String, String]]? cnv_output_dirs

        # Patient specific
        String patient
        Array[String] normal_samples
        Array[String] tumor_samples
        Map[String, Array[Array[String]]]? fastqs
        Map[String, Array[Int]]? lanes
        Map[String, Array[String]]? uuids
        Map[String, String]? libraries
        Array[Pair[String, Pair[String, String]]]? normal_bams
        Array[Pair[String, Pair[String, String]]]? normal_bais
        Array[Pair[String, Pair[String, String]]]? tumor_bams
        Array[Pair[String, Pair[String, String]]]? tumor_bais
        String sex

        # Other required inputs
        String repo_dir
        File sa_key = Defaults.sa_key
        String cloud_bucket = Defaults.cloud_bucket

        # Alignment
        String sequencing_platform = Defaults.sequencing_platform
        String sequencing_platform_unit = Defaults.sequencing_platform_unit
        String sequencing_center = Defaults.sequencing_center
        String alignment_reference = Defaults.alignment_reference
        String alignment_reference_index = Defaults.alignment_reference_index
        String alignment_reference_dict = Defaults.alignment_reference_dict
        String alignment_reference_amb = Defaults.alignment_reference_amb
        String alignment_reference_ann = Defaults.alignment_reference_ann
        String alignment_reference_bwt = Defaults.alignment_reference_bwt
        String alignment_reference_pac = Defaults.alignment_reference_pac
        String alignment_reference_sa = Defaults.alignment_reference_sa
        String known_sites = Defaults.known_sites
        String known_sites_index = Defaults.known_sites_index
        String gatk_intervals = Defaults.gatk_intervals

        # GRIDSS
        String gridss_purple_linx_data = Defaults.gridss_purple_linx_data

        # SSM/SV
        String reference = Defaults.reference
        String reference_index = Defaults.reference_index
        String reference_dict = Defaults.reference_dict
        String centromeres = Defaults.centromeres
        String dustmaker = Defaults.dustmaker

        # SSM
        Array[String] intervals = Defaults.intervals
        String build_version = Defaults.build_version
        String cosmic = Defaults.cosmic
        String cosmic_index = Defaults.cosmic_index
        String dbsnp = Defaults.dbsnp
        String dbsnp_index = Defaults.dbsnp_index
        String annovar_humandb_tar = Defaults.annovar_humandb_tar
        String annovar_bedfile = Defaults.annovar_bedfile
        String cosmic_rda = Defaults.cosmic_rda
        String source = Defaults.source
        String pon = Defaults.pon
        Int pon_max = Defaults.pon_max
        Int min_dist_complex = Defaults.min_dist_complex
        Int max_number_of_flagged = Defaults.max_number_of_flagged
        Int normal_alt_count_max = Defaults.normal_alt_count_max
        String common_biallelic_germline_variants = Defaults.common_biallelic_germline_variants
        String common_biallelic_germline_variants_index = Defaults.common_biallelic_germline_variants_index

        # SV
        Array[String] sv_types = Defaults.sv_types
        Map[String, String] delly_sv_types = Defaults.delly_sv_types
        String Delly_header = Defaults.Delly_header
        String delly_exclude_file = Defaults.delly_exclude_file
        Map[String, String] panel_of_normals_beds = Defaults.panel_of_normals_beds
        String cosmic_cancer_genes_rda = Defaults.cosmic_cancer_genes_rda
        String ucsc_info_tar = Defaults.ucsc_info_tar
        String circos_config_karyotype = Defaults.circos_config_karyotype
        String circos_config_ideogram_conf = Defaults.circos_config_ideogram_conf
        String circos_config_image_conf = Defaults.circos_config_image_conf
        String circos_config_colors_fonts_patterns_conf = Defaults.circos_config_colors_fonts_patterns_conf
        String circos_config_housekeeping_conf = Defaults.circos_config_housekeeping_conf
        String circos_config_ticks_conf = Defaults.circos_config_ticks_conf

        # CNV
        String impute_info = Defaults.impute_info
        String impute_files_tar = Defaults.impute_files_tar
        String thougen_loci_tar = Defaults.thougen_loci_tar
        String ignore_contig = Defaults.ignore_contig
        String prob_loci = Defaults.prob_loci
        String gc_correct_tar = Defaults.gc_correct_tar
    }

    call DefaultValues.Defaults as Defaults {
        input:
            environment = environment,
            random_number = random_number,
            repo_dir = repo_dir
    }

    if (defined(align_output_dirs)) {
        Map[String, String] align_output_dirs_input = select_first([align_output_dirs])
        Map[String, Array[Array[String]]] fastqs_input = select_first([fastqs])
        Map[String, Array[Int]] lanes_input = select_first([lanes])
        Map[String, Array[String]] uuids_input = select_first([uuids])
        Map[String, String] libraries_input = select_first([libraries])

        scatter (sample in normal_samples) {
            String normal_align_output_dir = if (environment == "LOCAL") then align_output_dirs_input[sample] else cloud_bucket + align_output_dirs_input[sample]
            call AlignNamespace.Align as AlignNormals {
                input:
                    output_dir = normal_align_output_dir,
                    patient = patient,
                    sample = sample,
                    tumor_normal = 'normal',
                    lanes = lanes_input[sample],
                    fastqs = fastqs_input[sample],
                    uuids = uuids_input[sample],
                    library = libraries_input[sample],
                    platform = sequencing_platform,
                    platform_unit = sequencing_platform_unit,
                    center = sequencing_center,
                    reference = alignment_reference,
                    reference_index = alignment_reference_index,
                    reference_dict = alignment_reference_dict,
                    reference_amb = alignment_reference_amb,
                    reference_ann = alignment_reference_ann,
                    reference_bwt = alignment_reference_bwt,
                    reference_pac = alignment_reference_pac,
                    reference_sa = alignment_reference_sa,
                    known_sites = known_sites,
                    known_sites_index = known_sites_index,
                    gatk_intervals = gatk_intervals,
                    repo_dir = repo_dir,
                    sa_key = sa_key,
                    environment = environment
            }

            call WriteWorkflowCompleteIndicator as AlignNormalsWorkflowComplete {
                input:
                    dependency = AlignNormals.bam,
                    output_dir = normal_align_output_dir,
                    sa_key = sa_key,
                    environment = environment
            }
        }

        scatter (sample in tumor_samples) {
            String tumour_align_output_dir = if (environment == "LOCAL") then align_output_dirs_input[sample] else cloud_bucket + align_output_dirs_input[sample]
            call AlignNamespace.Align as AlignTumors {
                input:
                    output_dir = tumour_align_output_dir,
                    patient = patient,
                    sample = sample,
                    tumor_normal = 'tumor',
                    lanes = lanes_input[sample],
                    fastqs = fastqs_input[sample],
                    uuids = uuids_input[sample],
                    library = libraries_input[sample],
                    platform = sequencing_platform,
                    platform_unit = sequencing_platform_unit,
                    center = sequencing_center,
                    reference = alignment_reference,
                    reference_index = alignment_reference_index,
                    reference_dict = alignment_reference_dict,
                    reference_amb = alignment_reference_amb,
                    reference_ann = alignment_reference_ann,
                    reference_bwt = alignment_reference_bwt,
                    reference_pac = alignment_reference_pac,
                    reference_sa = alignment_reference_sa,
                    known_sites = known_sites,
                    known_sites_index = known_sites_index,
                    gatk_intervals = gatk_intervals,
                    repo_dir = repo_dir,
                    sa_key = sa_key,
                    environment = environment
            }

            call WriteWorkflowCompleteIndicator as AlignTumorsWorkflowComplete {
                input:
                    dependency = AlignTumors.bam,
                    output_dir = tumour_align_output_dir,
                    sa_key = sa_key,
                    environment = environment
            }
        }
    }

    Array[Pair[String, Pair[String, String]]] tumor_bam_sample_pairs = if (defined(AlignTumors.sample_bam_pair)) then select_first([AlignTumors.sample_bam_pair]) else select_first([tumor_bams])
    Array[Pair[String, Pair[String, String]]] normal_bam_sample_pairs = if (defined(AlignNormals.sample_bam_pair)) then select_first([AlignNormals.sample_bam_pair]) else select_first([normal_bams])

    if (defined(ssm_output_dirs)) {
        scatter (tumor_pair in tumor_bam_sample_pairs) {
            String ssm_tumor_sample = tumor_pair.left
            String ssm_tumor_bam = tumor_pair.right.left
            String ssm_tumor_bai = sub(ssm_tumor_bam, ".bam$", ".bam.bai")
            String ssm_tumor_gatk_depth_of_coverage = sub(ssm_tumor_bam, ".bam$", ".gatk.depthofcoverage.sample_interval_summary")

            scatter (normal_pair in normal_bam_sample_pairs) {
                String ssm_normal_sample = normal_pair.left
                String ssm_normal_bam = normal_pair.right.left
                String ssm_normal_bai = sub(ssm_normal_bam, ".bam$", ".bam.bai")

                String ssm_output_dir_base = select_first([ssm_output_dirs])[ssm_tumor_sample][ssm_normal_sample]
                String ssm_output_dir = if (environment == "LOCAL") then ssm_output_dir_base else cloud_bucket + ssm_output_dir_base
                call SSMNamespace.SSM as SSM {
                    input:
                        output_dir = ssm_output_dir,
                        patient = patient,
                        normal_sample = ssm_normal_sample,
                        sample = ssm_tumor_sample,
                        tumor_bam = ssm_tumor_bam,
                        tumor_bai = ssm_tumor_bai,
                        normal_bam = ssm_normal_bam,
                        normal_bai = ssm_normal_bai,
                        tumor_gatk_depth_of_coverage = ssm_tumor_gatk_depth_of_coverage,
                        intervals = intervals,
                        reference = reference,
                        reference_index = reference_index,
                        reference_dict = reference_dict,
                        build_version = build_version,
                        cosmic = cosmic,
                        cosmic_index = cosmic_index,
                        dbsnp = dbsnp,
                        dbsnp_index = dbsnp_index,
                        annovar_humandb_tar = annovar_humandb_tar,
                        annovar_bedfile = annovar_bedfile,
                        cosmic_rda = cosmic_rda,
                        centromeres = centromeres,
                        dustmaker = dustmaker,
                        source = source,
                        pon = pon,
                        pon_max = pon_max,
                        min_dist_complex = min_dist_complex,
                        max_number_of_flagged = max_number_of_flagged,
                        normal_alt_count_max = normal_alt_count_max,
                        common_biallelic_germline_variants = common_biallelic_germline_variants,
                        common_biallelic_germline_variants_index = common_biallelic_germline_variants_index,
                        repo_dir = repo_dir,
                        sa_key = sa_key,
                        environment = environment
                }

                call WriteWorkflowCompleteIndicator as SSMWorkflowComplete {
                    input:
                        dependency = SSM.filtered_vcf,
                        output_dir = ssm_output_dir,
                        sa_key = sa_key,
                        environment = environment
                }

                if (run_GRIDSS) {
                    String gridss_tumor_bam = tumor_pair.right.right
                    String gridss_tumor_bai = sub(gridss_tumor_bam, ".bam$", ".bam.bai")

                    String gridss_normal_bam = normal_pair.right.right
                    String gridss_normal_bai = sub(gridss_normal_bam, ".bam$", ".bam.bai")

                    String gridss_output_dir_base = gridss_output_dirs[ssm_tumor_sample][ssm_normal_sample]
                    String gridss_output_dir = if (environment == "LOCAL") then gridss_output_dir_base else cloud_bucket + gridss_output_dir_base
                    call GRIDSSNamespace.GRIDSS as GRIDSS {
                        input:
                            output_dir = gridss_output_dir,
                            patient = patient,
                            sample = ssm_tumor_sample,
                            tumor_bam = gridss_tumor_bam,
                            tumor_bai = gridss_tumor_bai,
                            normal_bam = gridss_normal_bam,
                            normal_bai = gridss_normal_bai,
                            gridss_purple_linx_data = gridss_purple_linx_data,
                            snv_vcf = SSM.filtered_vcf,
                            reference = reference,
                            reference_index = reference_index,
                            reference_dict = reference_dict,
                            reference_amb = alignment_reference_amb,
                            reference_ann = alignment_reference_ann,
                            reference_bwt = alignment_reference_bwt,
                            reference_sa = alignment_reference_sa,
                            reference_pac = alignment_reference_pac,
                            repo_dir = repo_dir,
                            sa_key = sa_key,
                            environment = environment
                    }

                    call WriteWorkflowCompleteIndicator as GRIDSSWorkflowComplete {
                        input:
                            dependency = GRIDSS.output_archive,
                            output_dir = gridss_output_dir,
                            sa_key = sa_key,
                            environment = environment
                    }
                }
            }
        }
    }

    if (defined(sv_output_dirs)) {
        scatter (tumor_pair in tumor_bam_sample_pairs) {
            String sv_tumor_sample = tumor_pair.left
            String sv_tumor_bam = tumor_pair.right.left
            String sv_tumor_bai = sub(sv_tumor_bam, ".bam$", ".bam.bai")
            String sv_tumor_gatk_depth_of_coverage = sub(sv_tumor_bam, ".bam$", ".gatk.depthofcoverage.sample_interval_summary")

            scatter (normal_pair in normal_bam_sample_pairs) {
                String sv_normal_sample = normal_pair.left
                String sv_normal_bam = normal_pair.right.left
                String sv_normal_bai = sub(sv_normal_bam, ".bam$", ".bam.bai")

                String sv_output_dir_base = select_first([sv_output_dirs])[sv_tumor_sample][sv_normal_sample]
                String sv_output_dir = if (environment == "LOCAL") then sv_output_dir_base else cloud_bucket + sv_output_dir_base
                call SVNamespace.SV as SV {
                    input:
                        output_dir = sv_output_dir,
                        patient = patient,
                        sample = sv_tumor_sample,
                        tumor_bam = sv_tumor_bam,
                        tumor_bai = sv_tumor_bai,
                        normal_bam = sv_normal_bam,
                        normal_bai = sv_normal_bai,
                        tumor_gatk_depth_of_coverage = sv_tumor_gatk_depth_of_coverage,
                        sv_types = sv_types,
                        delly_sv_types = delly_sv_types,
                        reference = reference,
                        reference_index = reference_index,
                        reference_dict = reference_dict,
                        Delly_header = Delly_header,
                        delly_exclude_file = delly_exclude_file,
                        centromeres = centromeres,
                        dustmaker = dustmaker,
                        panel_of_normals_beds = panel_of_normals_beds,
                        cosmic_cancer_genes_rda = cosmic_cancer_genes_rda,
                        ucsc_info_tar = ucsc_info_tar,
                        circos_config_karyotype = circos_config_karyotype,
                        circos_config_ideogram_conf = circos_config_ideogram_conf,
                        circos_config_image_conf = circos_config_image_conf,
                        circos_config_colors_fonts_patterns_conf = circos_config_colors_fonts_patterns_conf,
                        circos_config_housekeeping_conf = circos_config_housekeeping_conf,
                        circos_config_ticks_conf = circos_config_ticks_conf,
                        repo_dir = repo_dir,
                        sa_key = sa_key,
                        environment = environment
                }

                call WriteWorkflowCompleteIndicator as SVWorkflowComplete {
                    input:
                        dependency = SV.final_rearrangements_tab,
                        output_dir = sv_output_dir,
                        sa_key = sa_key,
                        environment = environment
                }
            }
        }
    }

    if (defined(cnv_output_dirs)) {
        scatter (tumor_pair in tumor_bam_sample_pairs) {
            String cnv_tumor_sample = tumor_pair.left
            String cnv_tumor_bam = tumor_pair.right.left
            # CNV pipeline requires .bai indexes, NOT .bam.bai
            # Not true, it requires .bam.bai
            String cnv_tumor_bai = sub(cnv_tumor_bam, ".bam$", ".bam.bai")

            scatter (normal_pair in normal_bam_sample_pairs) {
                String cnv_normal_sample = normal_pair.left
                String cnv_normal_bam = normal_pair.right.left
                String cnv_normal_bai = sub(cnv_normal_bam, ".bam$", ".bam.bai")

                String cnv_output_dir_base = select_first([cnv_output_dirs])[cnv_tumor_sample][cnv_normal_sample]
                String cnv_output_dir = if (environment == "LOCAL") then cnv_output_dir_base else cloud_bucket + cnv_output_dir_base
                call CNVNamespace.CNV as CNV {
                    input:
                        output_dir = cnv_output_dir,
                        patient = patient,
                        sample = cnv_tumor_sample,
                        tumor_bam = cnv_tumor_bam,
                        tumor_bai = cnv_tumor_bai,
                        normal_bam = cnv_normal_bam,
                        normal_bai = cnv_normal_bai,
                        sex = sex,
                        reference = reference,
                        reference_index = reference_index,
                        reference_dict = reference_dict,
                        impute_info = impute_info,
                        impute_files_tar = impute_files_tar,
                        thougen_loci_tar = thougen_loci_tar,
                        ignore_contig = ignore_contig,
                        prob_loci = prob_loci,
                        gc_correct_tar = gc_correct_tar,
                        repo_dir = repo_dir,
                        sa_key = sa_key,
                        environment = environment
                }

                call WriteWorkflowCompleteIndicator as CNVWorkflowComplete {
                    input:
                        dependency = CNV.copynumber_csv,
                        output_dir = cnv_output_dir,
                        sa_key = sa_key,
                        environment = environment
                }
            }
        }
    }
}

task WriteWorkflowCompleteIndicator {
    input {
        String dependency
        String output_dir
        File sa_key
        String environment
    }

    command <<<
        set -eo pipefail
        if [ "~{environment}" = "LOCAL" ]; then
            touch ~{output_dir}/workflow_complete_indicator
            EXIT_STATUS=$?
            # Check that workflow complete indicator was successful
            if [ $EXIT_STATUS -eq 0 ]; then
                # Check to see if execution_dirs file exits
                echo "Exit status is good"
                if [ -f "~{output_dir}/execution_dirs.txt" ]; then
                    # remove all directories that have been used
                    echo "The execution_dirs file is present"
                    for d in `cat ~{output_dir}/execution_dirs.txt`; do [ ! -d "$d" ] || rm -r $d; done
                    # now remove the execution_dirs file itself.
                    echo "Trying to delete the execution file now."
                    rm ~{output_dir}/execution_dirs.txt
                fi
            fi
        elif [ "~{environment}" = "CLOUD" ]; then
            gcloud auth activate-service-account --key-file ~{sa_key}
            touch workflow_complete_indicator
            gsutil cp workflow_complete_indicator ~{output_dir}
        else
            echo "ERORR: environment must be one of ['LOCAL', 'CLOUD']. ~{environment} is not a valid value." >&2
            exit 1
        fi
    >>>

    runtime {
        memory: "4 GB"
        walltime: "23:59:59"
        task_name: "WriteWorkflowCompleteIndicator"
        log_dir: "~{output_dir}/wdl_logs+run_info"
    }
}
