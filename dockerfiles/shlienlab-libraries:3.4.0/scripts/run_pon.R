### run_pon.R #####################################################################################
# A R script that wraps the PON annotation as a stand alone script.

### HISTORY #######################################################################################
# Version           Date            Developer               Comments
# 0.01              2017-04-26      rdeborja                initial development

### NOTES #########################################################################################
#

### PREAMBLE ######################################################################################
library('getopt')

usage <- function() {
    usage.text <- '\nUsage: run_pon.R --rda test.dat --pon pon.rda --output output.rda\n\n'
    return(usage.text)
    }

params = matrix(
    c(
        'rda', 'r', 1, 'character',
        'pon', 'p', 1, 'character',
        'output', 'o', 1, 'character'
        ),
    ncol = 4,
    byrow = TRUE
    )

opt = getopt(params)

# verify arguments
if(is.null(opt$rda)) { stop(usage()) }
if (is.null(opt$pon)) { stop(usage) }
if (is.null(opt$output)) { opt$output <- 'output.rda' }

### LIBRARIES #####################################################################################
library(ShlienLab.Core.Filter)

### FUNCTIONS #####################################################################################

### GET DATA ######################################################################################
data_in <- load(opt$rda)
data_in <- get(data_in)

pon <- load(opt$pon)
pon <- get(pon)

### PROCESS DATA ##################################################################################
# check if the input dataframe has a snvid column, add one if it doesnt
if ('snvid' %in% names(data_in) == FALSE) {
    data_in <- ShlienLab.Core.Filter::add.snvid(data = data_in)
    }

data_in <- dplyr::left_join(
    x = data_in,
    y = data.frame(pon),
    by = c('snvid')
    )
data_in[is.na(data_in$pon_count),]$pon_count <- 0

data.filtered <- data_in
save(data.filtered, file=opt$output)

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################

### SESSION INFORMATION ###########################################################################
sessionInfo()
