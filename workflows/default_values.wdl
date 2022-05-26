version 1.0

workflow Defaults {
    input {
        String environment
        String random_number
        String repo_dir
    }

    # General
    String cloud_resource_bucket = PATH_TO_CLOUD_RESOURCE_BUCKET
    String cloud_output_bucket = PATH_TO_CLOUD_OUTPUT_BUCKET

    # Alignment
    String default_sequencing_platform = "ILLUMINA"
    String default_sequencing_platform_unit = "NONE"
    String default_sequencing_center = "TCAG"
    Boolean default_fix_misencoded_quality_scores = false

    # SSM
    String default_source = "WGS"
    String default_build_version = "hg19"
    Array[String] default_intervals = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    Int default_pon_max = 2
    Int default_normal_alt_count_max = 2 # Value changed to 2 on bitbucket
    Int default_min_dist_complex = 0
    Int default_max_number_of_flagged = 2

    # SV
    Array[String] default_sv_types = ["del", "dup", "inv", "tra"]
    Map[String, String] default_delly_sv_types = {"del": "DEL", "dup": "DUP", "inv": "INV", "tra": "TRA"}


    if (environment == "LOCAL") {
        String local_resource_dir = "~{repo_dir}/resources"
        String local_sa_key = sub("~{local_resource_dir}/dummy_key/RANDOM_NUMBER/dummy_sa_key.json", "RANDOM_NUMBER", random_number)

        # Alignment
        String local_alignment_reference = sub("~{local_resource_dir}/bwa_reference/RANDOM_NUMBER/hs37d5.fa", "RANDOM_NUMBER", random_number)
        String local_alignment_reference_index = sub("~{local_resource_dir}/bwa_reference/RANDOM_NUMBER/hs37d5.fa.fai", "RANDOM_NUMBER", random_number)
        String local_alignment_reference_dict = sub("~{local_resource_dir}/bwa_reference/RANDOM_NUMBER/hs37d5.dict", "RANDOM_NUMBER", random_number)
        String local_alignment_reference_amb = sub("~{local_resource_dir}/bwa_reference/RANDOM_NUMBER/hs37d5.fa.amb", "RANDOM_NUMBER", random_number)
        String local_alignment_reference_ann = sub("~{local_resource_dir}/bwa_reference/RANDOM_NUMBER/hs37d5.fa.ann", "RANDOM_NUMBER", random_number)
        String local_alignment_reference_bwt = sub("~{local_resource_dir}/bwa_reference/RANDOM_NUMBER/hs37d5.fa.bwt", "RANDOM_NUMBER", random_number)
        String local_alignment_reference_pac = sub("~{local_resource_dir}/bwa_reference/RANDOM_NUMBER/hs37d5.fa.pac", "RANDOM_NUMBER", random_number)
        String local_alignment_reference_sa = sub("~{local_resource_dir}/bwa_reference/RANDOM_NUMBER/hs37d5.fa.sa", "RANDOM_NUMBER", random_number)

        String local_known_sites = "~{local_resource_dir}/alignment_known_sites_and_intervals/dbsnp_138.b37.vcf"
        String local_known_sites_index = "~{local_resource_dir}/alignment_known_sites_and_intervals/dbsnp_138.b37.vcf.idx"
        String local_gatk_intervals = "~{local_resource_dir}/alignment_known_sites_and_intervals/hs37d5.intervals"

        # GRIDSS
        String local_gridss_purple_linx_data = "~{local_resource_dir}/gridss/gpl_ref_data_37.gz"

        # SSM / SV
        String local_reference = sub("~{local_resource_dir}/ssm_and_sv_reference/RANDOM_NUMBER/hs37d5.fa", "RANDOM_NUMBER", random_number)
        String local_reference_index = sub("~{local_resource_dir}/ssm_and_sv_reference/RANDOM_NUMBER/hs37d5.fa.fai", "RANDOM_NUMBER", random_number)
        String local_reference_dict = sub("~{local_resource_dir}/ssm_and_sv_reference/RANDOM_NUMBER/hs37d5.dict", "RANDOM_NUMBER", random_number)
        String local_centromeres = sub("~{local_resource_dir}/ssm_and_sv_reference/RANDOM_NUMBER/centromeres.tab", "RANDOM_NUMBER", random_number)
        String local_dustmaker = sub("~{local_resource_dir}/ssm_and_sv_reference/RANDOM_NUMBER/dustmasker60.bed", "RANDOM_NUMBER", random_number)

        # SSM
        # I'm not sure the dbsnp files are being used anymore, it looks like they've been superceded by the gnomad biallelic sites
        # Doesn't look like the CosmicCodingMuts files are being used anymore, the aren't being passed into SSM tasks
        String local_dbsnp = "~{local_resource_dir}/ssm/source/dbsnp_138.b37.vcf"
        String local_dbsnp_index = "~{local_resource_dir}/ssm/source/dbsnp_138.b37.vcf.idx"
        String local_cosmic = "~{local_resource_dir}/ssm/source/CosmicCodingMuts.vcf.gz"
        String local_cosmic_index = "~{local_resource_dir}/ssm/source/CosmicCodingMuts.vcf.gz.tbi"
        String local_annovar_humandb_tar = sub("~{local_resource_dir}/ssm/RANDOM_NUMBER/humandb.tar.gz", "RANDOM_NUMBER", random_number)
        String local_annovar_bedfile = "SureSelect_All_Exon_50mb_with_annotation_HG19_BED.removeChrUn.bed"
        String local_cosmic_rda = sub("~{local_resource_dir}/ssm/RANDOM_NUMBER/cosmic.cancer.gene.census.rda", "RANDOM_NUMBER", random_number)
        String local_pon = sub("~{local_resource_dir}/ssm/RANDOM_NUMBER/pon_count.rda", "RANDOM_NUMBER", random_number)
        String local_common_biallelic_germline_variants = sub("~{local_resource_dir}/ssm/RANDOM_NUMBER/biallelic_af_only_gnomad.genomes.r2.1.1.sites.vcf.gz", "RANDOM_NUMBER", random_number)
        String local_common_biallelic_germline_variants_index = sub("~{local_resource_dir}/ssm/RANDOM_NUMBER/biallelic_af_only_gnomad.genomes.r2.1.1.sites.vcf.gz.tbi", "RANDOM_NUMBER", random_number)

        # SV
        String local_Delly_header = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/DellyHeader.txt", "RANDOM_NUMBER", random_number)
        String local_delly_exclude_file = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/human.hg19.excl.tsv", "RANDOM_NUMBER", random_number)
        Map[String, String] local_panel_of_normals_beds = {"del": sub("~{local_resource_dir}/sv/RANDOM_NUMBER/del_ewings_normals.bed", "RANDOM_NUMBER", random_number),
                                                              "dup": sub("~{local_resource_dir}/sv/RANDOM_NUMBER/dup_ewings_normals.bed", "RANDOM_NUMBER", random_number),
                                                              "inv": sub("~{local_resource_dir}/sv/RANDOM_NUMBER/inv_ewings_normals.bed", "RANDOM_NUMBER", random_number),
                                                              "tra": sub("~{local_resource_dir}/sv/RANDOM_NUMBER/tra_ewings_normals.bed", "RANDOM_NUMBER", random_number)}
        String local_cosmic_cancer_genes_rda = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/cosmic_cancer_gene_census.rda", "RANDOM_NUMBER", random_number)
        String local_ucsc_info_tar = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/sv.tar.gz", "RANDOM_NUMBER", random_number)
        String local_circos_config_karyotype = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/karyotype.human.hg19.txt", "RANDOM_NUMBER", random_number)
        String local_circos_config_ideogram_conf = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/ideogram.conf", "RANDOM_NUMBER", random_number)
        String local_circos_config_image_conf = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/image.conf", "RANDOM_NUMBER", random_number)
        String local_circos_config_colors_fonts_patterns_conf = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/colors_fonts_patterns.conf", "RANDOM_NUMBER", random_number)
        String local_circos_config_housekeeping_conf = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/housekeeping.conf", "RANDOM_NUMBER", random_number)
        String local_circos_config_ticks_conf = sub("~{local_resource_dir}/sv/RANDOM_NUMBER/ticks.conf", "RANDOM_NUMBER", random_number)

        # CNV
        String local_impute_info = sub("~{local_resource_dir}/cnv/RANDOM_NUMBER/impute_info.txt", "RANDOM_NUMBER", random_number)
        String local_impute_files_tar = sub("~{local_resource_dir}/cnv/RANDOM_NUMBER/impute.tar.gz", "RANDOM_NUMBER", random_number)
        String local_thougen_loci_tar = sub("~{local_resource_dir}/cnv/RANDOM_NUMBER/1000genomesloci.tar.gz", "RANDOM_NUMBER", random_number)
        String local_ignore_contig = sub("~{local_resource_dir}/cnv/RANDOM_NUMBER/ignore.contigs.txt", "RANDOM_NUMBER", random_number)
        String local_prob_loci = sub("~{local_resource_dir}/cnv/RANDOM_NUMBER/probloci_270415.txt", "RANDOM_NUMBER", random_number)
        String local_gc_correct_tar = sub("~{local_resource_dir}/cnv/RANDOM_NUMBER/battenberg_wgs_gc_correction_1000g_v3.tar.gz", "RANDOM_NUMBER", random_number)
    }

    if (environment == "CLOUD") {
        String cloud_sa_key = cloud_resource_bucket + "sa-key.json"

        # Alignment
        String cloud_alignment_reference = cloud_resource_bucket + "hs37d5.fa"
        String cloud_alignment_reference_index = cloud_resource_bucket + "hs37d5.fa.fai"
        String cloud_alignment_reference_dict = cloud_resource_bucket + "hs37d5.dict"
        String cloud_alignment_reference_amb = cloud_resource_bucket + "hs37d5.fa.amb"
        String cloud_alignment_reference_ann = cloud_resource_bucket + "hs37d5.fa.ann"
        String cloud_alignment_reference_bwt = cloud_resource_bucket + "hs37d5.fa.bwt"
        String cloud_alignment_reference_pac = cloud_resource_bucket + "hs37d5.fa.pac"
        String cloud_alignment_reference_sa = cloud_resource_bucket + "hs37d5.fa.sa"
        String cloud_known_sites = cloud_resource_bucket + "dbsnp_138.b37.vcf"
        String cloud_known_sites_index = cloud_resource_bucket + "dbsnp_138.b37.vcf.idx"
        String cloud_gatk_intervals = cloud_resource_bucket + "hs37d5.intervals"

        # GRIDSS
        String cloud_gridss_purple_linx_data = cloud_resource_bucket + "gpl_ref_data_37.gz"

        # SSM / SV
        String cloud_reference = cloud_resource_bucket + "hs37d5.fa"
        String cloud_reference_index = cloud_resource_bucket + "hs37d5.fa.fai"
        String cloud_reference_dict = cloud_resource_bucket + "hs37d5.dict"
        String cloud_centromeres = cloud_resource_bucket + "centromeres.tab"
        String cloud_dustmaker = cloud_resource_bucket + "dustmasker60.bed"

        # SSM
        String cloud_dbsnp = cloud_resource_bucket + "dbsnp_138.b37.vcf"
        String cloud_dbsnp_index = cloud_resource_bucket + "dbsnp_138.b37.vcf.idx"
        String cloud_cosmic = cloud_resource_bucket + "CosmicCodingMuts.vcf.gz"
        String cloud_cosmic_index = cloud_resource_bucket + "CosmicCodingMuts.vcf.gz.tbi"
        String cloud_annovar_humandb_tar = cloud_resource_bucket + "humandb.tar.gz"
        String cloud_annovar_bedfile = cloud_resource_bucket + "SureSelect_All_Exon_50mb_with_annotation_HG19_BED.removeChrUn.bed"
        String cloud_cosmic_rda = cloud_resource_bucket + "cosmic.cancer.gene.census.rda"
        String cloud_pon = cloud_resource_bucket + "pon_count.rda"
        String cloud_common_biallelic_germline_variants = cloud_resource_bucket + "biallelic_af_only_gnomad.genomes.r2.1.1.sites.vcf.gz"
        String cloud_common_biallelic_germline_variants_index = cloud_resource_bucket + "biallelic_af_only_gnomad.genomes.r2.1.1.sites.vcf.gz.tbi"

        # SV
        String cloud_Delly_header = cloud_resource_bucket + "DellyHeader.txt"
        String cloud_delly_exclude_file = cloud_resource_bucket + "human.hg19.excl.tsv"
        Map[String, String] cloud_panel_of_normals_beds = {"del": cloud_resource_bucket + "del_ewings_normals.bed", "dup": cloud_resource_bucket + "dup_ewings_normals.bed", "inv": cloud_resource_bucket + "inv_ewings_normals.bed", "tra": cloud_resource_bucket + "tra_ewings_normals.bed"}
        String cloud_cosmic_cancer_genes_rda = cloud_resource_bucket + "cosmic_cancer_gene_census.rda"
        String cloud_ucsc_info_tar = cloud_resource_bucket + "UCSC_info.tar.gz"
        String cloud_circos_config_karyotype = cloud_resource_bucket + "karyotype.human.hg19.txt"
        String cloud_circos_config_ideogram_conf = cloud_resource_bucket + "ideogram.conf"
        String cloud_circos_config_image_conf = cloud_resource_bucket + "image.conf"
        String cloud_circos_config_colors_fonts_patterns_conf = cloud_resource_bucket + "colors_fonts_patterns.conf"
        String cloud_circos_config_housekeeping_conf = cloud_resource_bucket + "housekeeping.conf"
        String cloud_circos_config_ticks_conf = cloud_resource_bucket + "ticks.conf"

        # CNV
        String cloud_impute_info = cloud_resource_bucket + "impute_info.txt"
        String cloud_impute_files_tar = cloud_resource_bucket + "impute.tar.gz"
        String cloud_thougen_loci_tar = cloud_resource_bucket + "1000genomesloci.tar.gz"
        String cloud_ignore_contig = cloud_resource_bucket + "ignore.contigs.txt"
        String cloud_prob_loci = cloud_resource_bucket + "probloci_270415.txt"
        String cloud_gc_correct_tar = cloud_resource_bucket + "battenberg_wgs_gc_correction_1000g_v3.tar.gz"
    }

    output {
        # General cloud-specific
        String sa_key = select_first([local_sa_key, cloud_sa_key])
        String cloud_bucket = cloud_output_bucket

        # Alignment
        String alignment_reference = select_first([local_alignment_reference, cloud_alignment_reference])
        String alignment_reference_index = select_first([local_alignment_reference_index, cloud_alignment_reference_index])
        String alignment_reference_dict = select_first([local_alignment_reference_dict, cloud_alignment_reference_dict])
        String alignment_reference_amb = select_first([local_alignment_reference_amb, cloud_alignment_reference_amb])
        String alignment_reference_ann = select_first([local_alignment_reference_ann, cloud_alignment_reference_ann])
        String alignment_reference_bwt = select_first([local_alignment_reference_bwt, cloud_alignment_reference_bwt])
        String alignment_reference_pac = select_first([local_alignment_reference_pac, cloud_alignment_reference_pac])
        String alignment_reference_sa = select_first([local_alignment_reference_sa, cloud_alignment_reference_sa])
        String known_sites = select_first([local_known_sites, cloud_known_sites])
        String known_sites_index = select_first([local_known_sites_index, cloud_known_sites_index])
        String gatk_intervals = select_first([local_gatk_intervals, cloud_gatk_intervals])
        String sequencing_platform = default_sequencing_platform
        String sequencing_platform_unit = default_sequencing_platform_unit
        String sequencing_center = default_sequencing_center
        Boolean fix_misencoded_quality_scores = default_fix_misencoded_quality_scores

        # GRIDSS
        String gridss_purple_linx_data = select_first([local_gridss_purple_linx_data, cloud_gridss_purple_linx_data])

        # SSM / SV
        String reference = select_first([local_reference, cloud_reference])
        String reference_index = select_first([local_reference_index, cloud_reference_index])
        String reference_dict = select_first([local_reference_dict, cloud_reference_dict])
        String centromeres = select_first([local_centromeres, cloud_centromeres])
        String dustmaker = select_first([local_dustmaker, cloud_dustmaker])

        # SSM
        String source = default_source
        String build_version = default_build_version
        Array[String] intervals = default_intervals
        String dbsnp = select_first([local_dbsnp, cloud_dbsnp])
        String dbsnp_index = select_first([local_dbsnp_index, cloud_dbsnp_index])
        String cosmic = select_first([local_cosmic, cloud_cosmic])
        String cosmic_index = select_first([local_cosmic_index, cloud_cosmic_index])
        String annovar_humandb_tar = select_first([local_annovar_humandb_tar, cloud_annovar_humandb_tar])
        String annovar_bedfile = select_first([local_annovar_bedfile, cloud_annovar_bedfile])
        String cosmic_rda = select_first([local_cosmic_rda, cloud_cosmic_rda])
        String pon = select_first([local_pon, cloud_pon])
        Int pon_max = default_pon_max
        Int normal_alt_count_max = default_normal_alt_count_max
        Int min_dist_complex = default_min_dist_complex
        Int max_number_of_flagged = default_max_number_of_flagged
        String common_biallelic_germline_variants = select_first([local_common_biallelic_germline_variants, cloud_common_biallelic_germline_variants])
        String common_biallelic_germline_variants_index = select_first([local_common_biallelic_germline_variants_index, cloud_common_biallelic_germline_variants_index])

        # SV
        Array[String] sv_types = default_sv_types
        Map[String, String] delly_sv_types = default_delly_sv_types
        String Delly_header = select_first([local_Delly_header, cloud_Delly_header])
        String delly_exclude_file = select_first([local_delly_exclude_file, cloud_delly_exclude_file])
        Map[String, String] panel_of_normals_beds = select_first([local_panel_of_normals_beds, cloud_panel_of_normals_beds])
        String cosmic_cancer_genes_rda = select_first([local_cosmic_cancer_genes_rda, cloud_cosmic_cancer_genes_rda])
        String ucsc_info_tar = select_first([local_ucsc_info_tar, cloud_ucsc_info_tar])
        String circos_config_karyotype = select_first([local_circos_config_karyotype, cloud_circos_config_karyotype])
        String circos_config_ideogram_conf = select_first([local_circos_config_ideogram_conf, cloud_circos_config_ideogram_conf])
        String circos_config_image_conf = select_first([local_circos_config_image_conf, cloud_circos_config_image_conf])
        String circos_config_colors_fonts_patterns_conf = select_first([local_circos_config_colors_fonts_patterns_conf, cloud_circos_config_colors_fonts_patterns_conf])
        String circos_config_housekeeping_conf = select_first([local_circos_config_housekeeping_conf, cloud_circos_config_housekeeping_conf])
        String circos_config_ticks_conf = select_first([local_circos_config_ticks_conf, cloud_circos_config_ticks_conf])

        # CNV
        String impute_info = select_first([local_impute_info, cloud_impute_info])
        String impute_files_tar = select_first([local_impute_files_tar, cloud_impute_files_tar])
        String thougen_loci_tar = select_first([local_thougen_loci_tar, cloud_thougen_loci_tar])
        String ignore_contig = select_first([local_ignore_contig, cloud_ignore_contig])
        String prob_loci = select_first([local_prob_loci, cloud_prob_loci])
        String gc_correct_tar = select_first([local_gc_correct_tar, cloud_gc_correct_tar])
    }
}
