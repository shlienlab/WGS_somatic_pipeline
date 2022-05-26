### run_snv_postfiltering.R #######################################################################
# Run the somatic SNV post-filtering Rscript on mutation data

### HISTORY #######################################################################################
# Version           Date            Developer                                Comments
#------------------------------------------------------------------------------------
#    0.01     2016-10-13             rdeborja                     initial development
#    0.02     2019-07-31        Drew Thompson                   move to pipeline repo,
#                                                         bring in external functions
#    0.03     2019-08-22        Drew Thompson              Save filter info for graph
#    0.04     2019-12-09  Lisa-Monique Edward           allow for chromosomes with no 
#                                                                 unfiltered variants
#    0.05     2019-12-09  Lisa-Monique Edward         Fixed edge case: Single variant 
#                                                        Allele = T, column parsed as 
#                                                       logical rather than character
#    0.05     2019-12-09  Lisa-Monique Edward         Fixed edge case: Single variant 
#                                                                         is.na fails 
#    0.06     2020-10-01      Lisa-Monique Edward     Store and read txts vs rdas
#    0.07     2020-11-23      Lisa-Monique Edward         Store filtered variants for
#                                                                          annotation


### NOTES #########################################################################################
#

### PREAMBLE ######################################################################################
library('getopt')

usage <- function() {
    usage.text <- '\nUsage: run_snv_postfiltering.R --snv test.dat --cfilter cfilter.out --min_dist_complex 0 --max_number_of_flagged 2 --pon pon --output_file output_file --pon_max 2 --normal_alt_count_max 2 --filter_info filter_info --annot_file annot_file\n\n'
    return(usage.text)
    }

params = matrix(
    c(
        'snv', 'r', 1, 'character',
        'cfilter', 'c', 1, 'character',
        'min_dist_complex', 'm', 1, 'integer',
        'max_number_of_flagged', 'x', 1, 'integer',
        'pon', 'p', 1, 'character',
        'pon_max', 'o', 1, 'character',
        'output_file', 's', 1, 'character',
        'normal_alt_count_max', 'n', 1, 'character',
        'filter_info', 'f', 1, 'character',
        'annot_file', 'a', 1, 'character'
        ),
    ncol = 4,
    byrow = TRUE
    )

opt = getopt(params)

# verify arguments
if (is.null(opt$snv)) { stop(usage()) }
if (is.null(opt$cfilter)) { stop(usage()) }
if (is.null(opt$output_file)) { stop(usage()) }
if (is.null(opt$normal_alt_count_max)) { stop(usage()) }
if (is.null(opt$pon)) { stop(usage()) }
if (is.null(opt$min_dist_complex)) {
    opt$min_dist_complex <- 0
    }
if (is.null(opt$max_number_of_flagged)) {
    opt$max_number_of_flagged <- 2
    }
if (is.null(opt$pon_max)) {
    opt$pon_max <- 2
    }
if (is.null(opt$filter_info)) { stop(usage()) }

### LIBRARIES #####################################################################################
library("plyr")
library("dplyr")
library('readr')
library('tidyr')

### FUNCTIONS #####################################################################################
add.snvid <- function(data=NULL) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  data <- data %>% dplyr::mutate(snvid = paste(sep='_', annovar_chr, annovar_start, annovar_end, annovar_ref, annovar_alt))
  return(data)
  }

flag.complex.overlap.loci <- function(data=NULL, minimum=0, far.dist=999999) {
  if (is.null(data)) stop("Mandatory argument is missing")

  # NA values in the distance_to_low_complexity_1 are far from the complex region
  # set NA values of distance_to_low_complexity_1 to far.dist
  data <- data %>% mutate(distance_to_low_complexity_1 = replace_na(distance_to_low_complexity_1, far.dist))
  data$dist_low_complexity <- 0

  # had to add a try block, if the error happens we will just return data and no entry will be listed as 'FLAGGED'
  flagged_data <- try(data[abs(data$distance_to_low_complexity_1) <= minimum,]$dist_low_complexity <- 'FLAGGED', silent = TRUE)

  return(data)
  }

# this can be used within an apply statement to return
# the count of "FLAGGED" elements in a vector
add.flagged.count.column <- function(x) {
  count = 0
  for(i in 1:length(x)) {
    if (x[i] == 'FLAGGED') {
      count = count + 1
      }
    else {
      next()
      }
    }
  return(count)
  }

to_annot_table <- function(annots, fn) {
  annots %>% 
    select(annovar_chr, gt_POS) %>%
    rename(gt_start = gt_POS) %>%
    mutate(gt_end = gt_start,
      fn = fn)

}


cromwell_run_snv_filter_postprocessing <- function(data=NULL, cfilter=NULL, min.dist.complex=0, max_number_of_flagged=2, pon=NULL, pon_max=2, normal_alt_count_max=2, hardfilter=TRUE, filter_info) {
  if (is.null(data)) stop("Mandatory argument data is missing")
  if (is.null(cfilter)) stop("Mandatory argument cfilter is missing")
  if (is.null(pon)) stop("Mandatory argument pon is missing")

  min.dist.complex = as.integer(min.dist.complex)
  if (is.na(min.dist.complex)) stop("Please provide a valid value for min.dist.complex")
  max_number_of_flagged = as.integer(max_number_of_flagged)
  if (is.na(max_number_of_flagged)) stop("Please provide a valid value for max_number_of_flagged")
  pon_max = as.integer(pon_max)
  if (is.na(pon_max)) stop("Please provide a valid value for pon_max")
  normal_alt_count_max = as.integer(normal_alt_count_max)
  if (is.na(normal_alt_count_max)) stop("Please provide a valid value for normal_alt_count_max")

  annots <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(annots) <- c("annovar_chr", "gt_start", "gt_end", "fn")
  
  data <- add.snvid(data=data)
  data.cfilter <- read.table(
    file=cfilter,
    header=TRUE,
    as.is=TRUE,
    sep='\t',
    quote="\"",
    colClasses=c("Ref.Allele"="character", "Alt.Allele"="character")
    )
  data.cfilter <- data.cfilter %>% dplyr::mutate(
    snvid=paste(
      sep='_',
      Chromosome1,
      Position1,
      Position1,
      Ref.Allele,
      Alt.Allele
      )
    )
  data.cfilter <- data.cfilter %>% dplyr::select(
    snvid,
    in_centromere,
    normal_coverage_threshold,
    unique_mapping,
    high_depth,
    distance_to_low_complexity_1,
    distance_to_low_complexity_2,
    multi_mapping
    )

  # change the distance_to_low_complexity_1 values to "FLAGGED" or 0
  data.cfilter <- flag.complex.overlap.loci(
    data=data.cfilter,
    minimum=min.dist.complex
    )
  data.cfilter$flagged_count <- apply(
    X=data.cfilter[,c('in_centromere', 'normal_coverage_threshold', 'unique_mapping', 'high_depth', 'multi_mapping', 'dist_low_complexity')],
    MARGIN=1,
    FUN=add.flagged.count.column
    )

  data <- dplyr::left_join(x=data, y=data.cfilter, by=c('snvid'))

  # threshold for flagged items
  if (hardfilter == TRUE) {
    snv_cfilter <- data %>% filter(flagged_count >= max_number_of_flagged)
    data <- data %>% filter(flagged_count < max_number_of_flagged)
    }
  annots <- rbind(annots, to_annot_table(snv_cfilter, "cFilter"))
  cat("SNV\tcFilter\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  # load the panel of normal count table to annotate the somatic SNV data
  load(pon)
  data <- left_join(x=data, y=data.frame(pon_count), by=c('snvid'))

  # NAs are introduced if a snvid has no entry in the pon_count vector
  data <- data %>% mutate(pon_count = replace_na(pon_count, 0))

  # threshold for recurrence count
  if (hardfilter == TRUE) {
    snv_pon <- data %>% filter(pon_count >= pon_max)
    data <- data %>% filter(pon_count < pon_max)
    }
  annots <- rbind(annots, to_annot_table(snv_pon, "SNV_Post_PON"))
  cat("SNV\tPON\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  if (hardfilter == TRUE) {
    snv_alt_count <- data %>% filter(gatk_normal_alt_count >= normal_alt_count_max)
    data <- data %>% filter(gatk_normal_alt_count < normal_alt_count_max)
    }
  annots <- rbind(annots, to_annot_table(snv_alt_count, "SNV_Post_Normal_Alt_Count"))
  cat("SNV\tnormal_alt_count\t", nrow(data), "\n", sep = "", file = filter_info, append = TRUE)

  # cleanup the dataframe, remove extraneous columns
  data$distance_to_low_complexity_2 <- NULL

  return(list(data, annots))
  }

### GET DATA ######################################################################################
snv <- read.table(
    file=opt$snv,
    header=TRUE,
    as.is=TRUE,
    sep='\t',
    quote="\"",
    colClasses=c("annovar_ref"="character", "annovar_alt"="character", "g_REF"="character", "t_REF"="character", "t_ALT"="character")
    )

### PROCESS DATA ##################################################################################
if (nrow(snv) > 0 ) {
  snv_postprocessing_output <- cromwell_run_snv_filter_postprocessing(
    data = snv,
    cfilter = opt$cfilter,
    pon = opt$pon,
    min.dist.complex = opt$min_dist_complex,
    max_number_of_flagged = opt$max_number_of_flagged,
    pon_max = opt$pon_max,
    normal_alt_count_max = opt$normal_alt_count_max,
    filter_info = opt$filter_info
    )
  data.filtered <- snv_postprocessing_output[[1]]
  annots <- snv_postprocessing_output[[2]]
  prev_annots <- read_tsv(opt$annot_file)
  annots <- rbind(prev_annots, annots) %>% 
    mutate(annovar_chr = factor(annovar_chr, levels = c(1:22, "X", "Y"))) %>% 
    arrange(annovar_chr, gt_start)
  write_delim(annots, path=opt$annot_file, delim="\t")
} else {
  data.filtered <- data.frame(matrix(ncol=(ncol(snv) + 10), nrow=0))
  colnames(data.filtered) <- c(colnames(snv), "snvid", "in_centromere", "normal_coverage_threshold", "unique_mapping", "high_depth", "distance_to_low_complexity_1", "multi_mapping", "dist_low_complexity", "flagged_count", "pon_count")
}

write_tsv(x=data.filtered, path=opt$output_file)

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################

### SESSION INFORMATION ###########################################################################
sessionInfo()

