### run_snv_postfiltering.R #######################################################################
# Run the somatic SNV post-filtering Rscript on mutation data

### HISTORY #######################################################################################
# Version           Date            Developer               Comments
# 0.01              2016-10-13      rdeborja                initial development

### NOTES #########################################################################################
#

### PREAMBLE ######################################################################################
library('getopt')

usage <- function() {
    usage.text <- '\nUsage: run_snv_postfiltering.R --rda test.dat --cfilter cfilter.out --min_dist_complex 0 --max_number_of_flagged 2 --pon pon --sample sample\n\n'
    return(usage.text)
    }

params = matrix(
    c(
        'rda', 'r', 1, 'character',
        'cfilter', 'c', 1, 'character',
        'min_dist_complex', 'm', 1, 'integer',
        'max_number_of_flagged', 'x', 1, 'integer',
        'pon', 'p', 1, 'character',
        'pon_max', 'o', 1, 'character',
        'sample', 's', 1, 'character'
        ),
    ncol = 4,
    byrow = TRUE
    )

opt = getopt(params)

# verify arguments
if (is.null(opt$rda)) { stop(usage()) }
if (is.null(opt$cfilter)) { stop(usage()) }
if (is.null(opt$sample)) { stop(usage()) }
if (is.null(opt$pon)) {
    opt$pon = '/hpf/largeprojects/adam/projects/icgc_tcga_datasets/ewing_sarcoma_tirode/ref/pon_count.rda'
    }
if (is.null(opt$min_dist_complex)) {
    opt$min_dist_complex <- 0
    }
if (is.null(opt$max_number_of_flagged)) {
    opt$max_number_of_flagged <- 2
    }
if (is.null(opt$pon_max)) {
    opt$pon_max <- 2
    }

### LIBRARIES #####################################################################################
library(ShlienLab.Core.SNV)
library(ShlienLab.Core.Filter)

### FUNCTIONS #####################################################################################

### GET DATA ######################################################################################
snv <- load(opt$rda)
snv <- get(snv)

### PROCESS DATA ##################################################################################
data.filtered <- ShlienLab.Core.Filter::run_snv_filter_postprocessing(
    data = snv,
    cfilter = opt$cfilter,
    pon = opt$pon,
    min.dist.complex = opt$min_dist_complex,
    max_number_of_flagged = opt$max_number_of_flagged,
    pon_max = opt$pon_max
    )
output.filtered.file <- paste(
    sep='_',
    opt$sample,
    'cfilter',
    'pon.rda'
    )
save(
     x=data.filtered,
     file=output.filtered.file
     )

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################

### SESSION INFORMATION ###########################################################################
sessionInfo()

