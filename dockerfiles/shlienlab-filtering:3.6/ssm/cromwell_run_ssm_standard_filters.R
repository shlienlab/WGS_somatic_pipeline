### run_ssm_standard_filters.R ####################################################################
# Run the standard simple somatic filters on MuTect2 output.  Output will be separate for SNVs and
# indels.

### HISTORY #######################################################################################
# Version           Date            Developer                                Comments
#------------------------------------------------------------------------------------
#    0.02     2017-04-13             rdeborja        removed MT from dataframe due to
#                                                               clip filtering issues
#    0.03     2019-07-31        Drew Thompson                   move to pipeline repo,
#                                                         bring in external functions
#    0.04     2019-08-22        Drew Thompson              Save filter info for graph
#    0.05     2019-12-09      Lisa-Monique Edward     separate filtering and merging steps
#    0.06     2020-10-01      Lisa-Monique Edward     Store and read txts vs rdas
#    0.07     2020-11-23      Lisa-Monique Edward         Store filtered variants for
#                                                                          annotation
#    0.08     2020-12-12      Lisa-Monique Edward     Dynamic sample names and sources

### NOTES #########################################################################################
#

### PREAMBLE ######################################################################################
usage <- function() {
    usage.text <- '\nUsage: cromwell_run_ssm_standard_filters.R <cosmic rda> <WGS|WXS|CPANEL> <snv+indel output> <snv output> <pype output> <indel output> <filter info file> <input file> <vcf file> <tumor sample name> <normal sample name>\n\n'
    return(usage.text)
    }

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 12) stop(usage())
cosmic_rda = args[1]
source = args[2]
all_output = args[3]
snv_output = args[4]
pype_output = args[5]
indel_output = args[6]
filter_info = args[7]
annot_file = args[8]
annovar_file = args[9]
vcf_file = args[10]
tumour_sample = args[11]
normal_sample = args[12]

### LIBRARIES #####################################################################################
library("plyr")
library("dplyr")
library("vcfR")
library("stringr")
library("readr")

### FUNCTIONS #####################################################################################
get.mutect2.annotated.header <- function() {
  header <- c(
    'annovar_chr',
    'annovar_start',
    'annovar_end',
    'annovar_ref',
    'annovar_alt',
    'annovar_func',
    'annovar_gene',
    'annovar_exonic_func',
    'annovar_annotation',
    'annovar_ens_func',
    'annovar_ens_gene',
    'annovar_ens_exonic_func',
    'annovar_ens_annotation',
    'annovar_dbsnp',
    'annovar_1000g',
    'annovar_esp',
    'annovar_complete_genomics',
    'annovar_cosmic',
    'annovar_clinvar',
    'annovar_exac',
    'annovar_target',
    'normal_name',
    'tumour_name',
    'gatk_filter',
    'gatk_normal_lod',
    'gatk_tumour_lod',
    'gatk_normal_genotype',
    'gatk_normal_ref_count',
    'gatk_normal_alt_count',
    'gatk_normal_depth',
    'gatk_normal_allele_fraction',
    'gatk_tumour_genotype',
    'gatk_tumour_ref_count',
    'gatk_tumour_alt_count',
    'gatk_tumour_depth',
    'gatk_tumour_allele_fraction',
    'gatk_variant_length',
    'gatk_mutation_type'
    )
  return(header)
  }

cromwell_get.mutect2.data <- function(annovar_file, load_vcfs = TRUE, vcf_file = "", t_gt = "TUMOR", n_gt = "NORMAL") {
  # Annovar files to process
  
  # Two dataframes that will be vertically concatenated to give the result
  if (load_vcfs){
    if (!exists("vcf_file") | vcf_file == "") {
      stop("Please specify path to vcf for processing or set load_vcfs to FALSE.")
    }
    vcf_data <- data.frame()
  }
  
  # Load the annovar results
  annovar_data <- try(
    read.table(file = annovar_file, header = FALSE, as.is = TRUE, sep = '\t', quote = "\"", skip = 1),  
    silent = TRUE
  )
  # the "try" block will return a class of "try-error"
  if (class(annovar_data) == 'try-error') {
    if (grepl("no lines available in input", geterrmessage())) {
      annovar_data <- data.frame(matrix(NA, ncol=38, nrow=0))
    } else {
      stop("Error loading annovar_file")
    }
  }
  
  if (load_vcfs) {
    # Load the results from the VCF
    chr_vcf_data <- try(vcfR::vcfR2tidy(vcfR::read.vcfR(vcf_file, verbose=FALSE), verbose=FALSE))
    if (class(chr_vcf_data) == 'try-error') {
      # If any VCF files are missing, then we don't load any of them.
      stop("Error reading VCF file.")
    }
    else {
      t <- chr_vcf_data$gt[chr_vcf_data$gt$Indiv == t_gt,]
      g <- chr_vcf_data$gt[chr_vcf_data$gt$Indiv == n_gt,]
      
      # Clean up each of these dataframes
      for (gt in c('g', 't')) {
        split <- stringr::str_split_fixed(get(gt)$gt_GT_alleles, '/', 2)
        
        assign(gt, get(gt) %>%
          dplyr::select(-gt_GT_alleles, -ChromKey, -Indiv) %>% 
          dplyr::rename(gt_POS = POS) %>%
          dplyr::mutate(gt_REF = split[,1], gt_ALT = split[,2]))
        
      }
      # Change gt_ prefix to g_ or t_ prefix
      colnames(g) <- sub('^gt(?=_)', 'g', colnames(g), perl=TRUE)
      colnames(t) <- sub('^gt(?=_)', 't', colnames(t), perl=TRUE)
      # Add g and t to the growing data of VCF FORMAT data
      vcf_data <- rbind(vcf_data, cbind(g %>% dplyr::rename(gt_POS=g_POS), t %>% dplyr::select(-t_POS)))
    }
  }
  
  colnames(annovar_data) <- get.mutect2.annotated.header()
  if(load_vcfs)
    return(cbind(annovar_data, vcf_data))
  else
    return(annovar_data)
}

remove.multiple.genes <- function(data) {
  first.gene <- unlist(strsplit(x = data, split = '[(,;]', perl = TRUE))[1];
  return(first.gene);
  }

add.aminoacid.column <- function(data) {
  aa <- unlist(strsplit(x = data, split = "[:,]", perl = TRUE))[5]
  return(aa);
  }

filter.for.cosmic.cancer.gene.census <- function(data, cosmic_rda) {
  if (!exists('cosmic.cancer.gene.census')) {
    load(cosmic_rda);
    }
  return(data %in% cosmic.cancer.gene.census$gene);
  }

annotate.mutect2.data <- function(data=NULL, cosmic_rda) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  data$ensembl_gene <- apply(X=as.matrix(data$annovar_ens_gene), MARGIN=1, FUN=remove.multiple.genes)
  data$hgnc_gene <- apply(X=as.matrix(data$annovar_gene), MARGIN=1, FUN=remove.multiple.genes)
  data$cosmic_census <- filter.for.cosmic.cancer.gene.census(data=data$hgnc_gene, cosmic_rda)
  data$aa <- ""
  if (length(which(data$annovar_annotation != "")) > 0)
    data$aa <- apply(X=as.matrix(data$annovar_annotation), MARGIN=1, FUN=add.aminoacid.column)
  return(data)
  }


to_annot_table <- function(annots, fn) {
  annots %>% 
    select(annovar_chr, gt_POS) %>%
    rename(gt_start = gt_POS) %>%
    mutate(gt_end = gt_start,
      fn = fn)

}


filter_indel <- function(data=NULL, source='WGS', coverage=TRUE, indels=c('ins', 'del'), exac=0.001, filter_info) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  # only process PASS indel data (i.e. gatk_mutation_type is "ins" or "del")
  data <- data %>%
    filter(gatk_filter %in% c('PASS', 'clustered_events')) %>%
    filter(gatk_mutation_type %in% indels)

  cat("INDEL\tMutect2_filters\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  annots <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(annots) <- c("annovar_chr", "gt_start", "gt_end", "fn")
  
  # we haven't done a detailed analysis of depth of coverage for indels, regardless of the library type
  # we will use the same coverage cutoffs until further investigation
  if (coverage == TRUE) {
    if (source %in% c('WXS', 'WGS', 'CPANEL')) {
      indel_coverage <- to_annot_table(
        filter(data, gatk_normal_depth < 10 | gatk_tumour_depth < 20), 
        "SSM_INDEL_coverage"
      )
      annots <- rbind(annots, indel_coverage)

      data <- data %>%
        filter(gatk_normal_depth >= 10) %>%
        filter(gatk_tumour_depth >= 20)
    } else {
      stop("Invalid source, must be either WGS, WXS, or CPANEL")
    }
  }
  #data <- data %>% filter(gatk_normal_depth >= normal.depth & gatk_tumour_depth >= tumour.depth)

  cat("INDEL\tcoverage\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  # if ('annovar_clinvar' %in% names(data)) {
  #   indel_annotations <- data %>%
  #     filter(!is.na(annovar_dbsnp) & annovar_dbsnp != '' & annovar_dbsnp != 0 & !grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
  #     filter(!is.na(annovar_complete_genomics) & annovar_complete_genomics == '' & annovar_complete_genomics != 0 & !grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
  #     filter(!is.na(annovar_1000g) & annovar_1000g != '' & annovar_1000g != 0 & !grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
    
  #   data <- data %>%
  #     filter(is.na(annovar_dbsnp) | annovar_dbsnp == '' | annovar_dbsnp == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
  #     # fixed annovar_complete_genomics != ''
  #     filter(is.na(annovar_complete_genomics) | annovar_complete_genomics == '' | annovar_complete_genomics == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
  #     filter(is.na(annovar_1000g) | annovar_1000g == '' | annovar_1000g == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
  # } else {
  #   indel_annotations <- data %>%
  #     filter(dbsnp_site == 'DBSNP') %>%
  #     filter(!is.na(annovar_dbsnp) & annovar_dbsnp != '' & annovar_dbsnp != 0) %>%
  #     filter(!is.na(annovar_complete_genomics) & annovar_complete_genomics != '' & annovar_complete_genomics != 0) %>%
  #     filter(!is.na(annovar_1000g) & annovar_1000g != '' & annovar_1000g != 0)
    
  #   data <- data %>%
  #     filter(dbsnp_site != 'DBSNP') %>%
  #     filter(is.na(annovar_dbsnp) | annovar_dbsnp == '' | annovar_dbsnp == 0) %>%
  #     filter(is.na(annovar_complete_genomics) | annovar_complete_genomics == '' | annovar_complete_genomics == 0) %>%
  #     filter(is.na(annovar_1000g) | annovar_1000g == '' | annovar_1000g == 0)
  # }

  # # there is an exome specific filter that can be applied using the ExAC and ESP databases
  # if ('annovar_exac' %in% names(data) & source == 'WXS') {
  #   indel_annotations <- data %>%
  #     filter(annovar_exac > exac & !is.na(annovar_exac)) %>%
  #     filter(!is.na(annovar_esp) & annovar_esp != '' & annovar_esp != 0 & !grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
    
  #   data <- data %>%
  #     filter(annovar_exac <= exac | is.na(annovar_exac)) %>%
  #     filter(is.na(annovar_esp) | annovar_esp == '' | annovar_esp == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
  # }
  # annots <- rbind(annots, to_annot_table(indel_annotations, "SSM_INDEL_annotation_filters"))
  # cat("INDEL\tannotation_filters\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  return(list(data, annots))
  }

filter_snv <- function(data=NULL, coverage=TRUE, source='WGS', exac=0.001, vaf=0.0, filter_info) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  # check for valid entries for the source argument
  if (!(source %in% c('WGS', 'WXS', 'CPANEL')) & coverage == TRUE) stop("Invalid source argument, must be one of WGS, WXS, or CPANEL")

  # use MuTect's internal filter and keep those listed as "PASS"
  data <- data %>%
    filter(gatk_filter %in% c('PASS', 'clustered_events')) %>%
    filter(gatk_mutation_type == 'snv')

  cat("SNV\tMutect2_filters\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  annots <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(annots) <- c("annovar_chr", "gt_start", "gt_end", "fn")
  
  if (coverage == TRUE) {
    if (source == 'WXS') {
      snv_coverage <- filter(data, gatk_normal_depth < 20 | gatk_tumour_depth < 30)
      data <- data %>%
        filter(gatk_normal_depth >= 20) %>%
        filter(gatk_tumour_depth >= 30)
    } else if (source == 'WGS') {
      snv_coverage <- filter(data, gatk_normal_depth < 10 | gatk_tumour_depth < 10)
      data <- data %>%
        filter(gatk_normal_depth >= 10) %>%
        filter(gatk_tumour_depth >= 10)
    } else if (source == 'CPANEL') {
      snv_coverage <- filter(data, gatk_normal_depth < 50 | gatk_tumour_depth < 50)
      data <- data %>%
        filter(gatk_normal_depth >= 50) %>%
        filter(gatk_tumour_depth >= 50)
    } else {
      stop("Invalid source, must be either WGS, WXS, or CPANEL")
      }
  }
  annots <- rbind(annots, to_annot_table(snv_coverage, "SSM_SNV_coverage"))
  cat("SNV\tcoverage\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  # # check if clinvar was used to filter the data (older data sets may not have this)
  # if ('annovar_clinvar' %in% names(data)) {
  #   snv_annotations <- data %>%
  #     filter(!is.na(annovar_dbsnp) & annovar_dbsnp != '' & annovar_dbsnp != 0 & !grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
  #     filter(!is.na(annovar_complete_genomics) & annovar_complete_genomics == '' & annovar_complete_genomics != 0 & !grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
  #     filter(!is.na(annovar_1000g) & annovar_1000g != '' & annovar_1000g != 0 & !grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
    
  #   data <- data %>%
  #     filter(is.na(annovar_dbsnp) | annovar_dbsnp == '' | annovar_dbsnp == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
  #     # fixed annovar_complete_genomics != ''
  #     filter(is.na(annovar_complete_genomics) | annovar_complete_genomics == '' | annovar_complete_genomics == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
  #     filter(is.na(annovar_1000g) | annovar_1000g == '' | annovar_1000g == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
  # } else {
  #   snv_annotations <- data %>%
  #     filter(dbsnp_site == 'DBSNP') %>%
  #     filter(!is.na(annovar_dbsnp) & annovar_dbsnp != '' & annovar_dbsnp != 0) %>%
  #     filter(!is.na(annovar_complete_genomics) & annovar_complete_genomics != '' & annovar_complete_genomics != 0) %>%
  #     filter(!is.na(annovar_1000g) & annovar_1000g != '' & annovar_1000g != 0)
    
  #   data <- data %>%
  #     filter(dbsnp_site != 'DBSNP') %>%
  #     filter(is.na(annovar_dbsnp) | annovar_dbsnp == '' | annovar_dbsnp == 0) %>%
  #     filter(is.na(annovar_complete_genomics) | annovar_complete_genomics == '' | annovar_complete_genomics == 0) %>%
  #     filter(is.na(annovar_1000g) | annovar_1000g == '' | annovar_1000g == 0)
  #   }

  # # there is an exome specific filter that can be applied using the ExAC and ESP databases
  # if ('annovar_exac' %in% names(data) & source == 'WXS') {
  #   snv_annotations <- data %>%
  #     filter(annovar_exac > exac & !is.na(annovar_exac)) %>%
  #     filter(!is.na(annovar_esp) & annovar_esp != '' & annovar_esp != 0 & !grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
    
  #   data <- data %>%
  #     filter(annovar_exac <= exac | is.na(annovar_exac)) %>%
  #     filter(is.na(annovar_esp) | annovar_esp == '' | annovar_esp == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
  #   }
  # annots <- rbind(annots, to_annot_table(snv_annotations, "SSM_SNV_annotation_filters"))
  # cat("SNV\tannotation_filters\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  # filter variants based on their variant allele fraction
  snv_vaf <- data %>% filter(gatk_tumour_allele_fraction < vaf)
  data <- data %>% filter(gatk_tumour_allele_fraction >= vaf)
  
  annots <- rbind(annots, to_annot_table(snv_vaf, "SSM_SNV_VAF"))
  cat("SNV\tVAF\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  return(list(data, annots))
  }

### GET DATA ######################################################################################
data <- cromwell_get.mutect2.data(annovar_file, vcf_file = vcf_file, t_gt = tumour_sample, n_gt = normal_sample)

### PROCESS DATA ##################################################################################
# add additional annotations to the dataframe
if (nrow(data) > 0) {
  data <- annotate.mutect2.data(data=data, cosmic_rda)
} else {
  data <- data %>% mutate(ensembl_gene = "", hgnc_gene = "", cosmic_census = "", aa = "")
}

# currently there is a bug in downstream filtering causing a pre-filter step to remove the
# mitochondrial DNA from the output
data <- data %>% filter(annovar_chr != 'MT')
data <- data %>% filter(annovar_chr != 'M')

annots <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(annots) <- c("annovar_chr", "gt_start", "gt_end", "fn")

cat("SNV+INDEL\tunfiltered\t", nrow(data), "\n", sep = "", file = filter_info)

# separately filter the snv and indel data
cat("INDEL\tunfiltered\t", nrow(data[data$gatk_mutation_type %in% c('ins', 'del'), ]), "\n", sep = "", file = filter_info, append = TRUE)
indel_filter_output <- filter_indel(data=data, source=source, filter_info = filter_info)
data.indel.filtered <- indel_filter_output[[1]]
annots <- rbind(annots, indel_filter_output[[2]])

cat("SNV\tunfiltered\t", nrow(data[data$gatk_mutation_type == "snv", ]), "\n", sep = "", file = filter_info, append = TRUE)
snv_filter_output <- filter_snv(data=data, source=source, filter_info = filter_info)
data.snv.filtered <- snv_filter_output[[1]]
annots <- rbind(annots, snv_filter_output[[2]]) %>% mutate(annovar_chr = factor(annovar_chr, levels = c(1:22, "X", "Y"))) %>% arrange(annovar_chr, gt_start)

# pype modifications
pype.mutect.data <- data.snv.filtered %>%
  mutate(trinuc = 'none', mutation.type = 'none') %>% 
  select(Chromosome1=annovar_chr,
    Chromosome2=annovar_chr,
    Position1=annovar_start,
    Position2=annovar_end,
    Ref.Allele=annovar_ref,
    Alt.Allele=annovar_alt,
    Func=annovar_func,
    Gene=annovar_gene,
    VAF=gatk_tumour_allele_fraction,
    mutation.type,
    trinuc,
    n_ref_count=gatk_normal_ref_count,
    n_alt_count=gatk_normal_alt_count,
    t_ref_count=gatk_tumour_ref_count,
    t_alt_count=gatk_tumour_alt_count
  )

# save the dataframes
write.table(x=data, file=all_output, sep='\t', quote=FALSE, row.names=FALSE)
write.table(x=data.snv.filtered, file=snv_output, sep='\t', quote=FALSE, row.names=FALSE)
write.table(x=pype.mutect.data, file=pype_output, sep='\t', quote=FALSE, row.names=FALSE)
write.table(x=data.indel.filtered, file=indel_output, sep='\t', quote=FALSE, row.names=FALSE)
write_delim(annots, path=annot_file, delim="\t")

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################

### SESSION INFORMATION ###########################################################################
sessionInfo()

