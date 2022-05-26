### create_pon.R ##################################################################################
# Create a panel of normal (PON) file containing a set of normal mutation calls.

### HISTORY #######################################################################################
# Version           Date            Developer               Comments
# 0.01              2017-04-24      rdeborja                initial development\

### NOTES #########################################################################################
#

### PREAMBLE ######################################################################################
library('getopt')

usage <- function() {
    usage.text <- '\nUsage: create_pon.R --path test.dat --recurrence 2 --samples samples.txt\n\n'
    return(usage.text)
    }

params = matrix(
    c(
        'path', 'p', 1, 'character',
        'recurrence', 'r', 1, 'integer',
        'samples', 's', 1, 'character'
        ),
    ncol = 4,
    byrow = TRUE
    )

opt = getopt(params)

# the default arguments
if(is.null(opt$path)) {
    opt$path <- '/hpf/largeprojects/adam/projects/icgc_tcga_datasets/ewing_sarcoma_tirode/EGAD00001001051/PON_SNV/'
    }
if (is.null(opt$recurrence)) {
    opt$recurrence <- 2
    }
if (is.null(opt$samples)) { stop(usage()) }

samples <- read.table(
    file = opt$samples,
    header = FALSE,
    as.is = TRUE,
    sep = '\t',
    quote = "\""
    )
colnames(samples) <- 'name'
suffix <- '_annotated_filtered.rda'

### LIBRARIES #####################################################################################
library(devtools)
load_all(pkg='~/ShlienLab.Core.Filter')

### FUNCTIONS #####################################################################################

### GET DATA ######################################################################################
# get the list of files to use for each iteration
files <- vector()
for(i in 1:nrow(samples)) {
    files[i] <- paste(sep='',
        opt$path, '/',
        samples$name[i],
        suffix
        )
    }

for (i in 2:length(files)) {
    output <- paste(sep='_',
        i,
        'pon_count.rda'
        )

    pon_data <- create.panel.of.normal.dataframe(files=files[1:i])
    pon_count <- get.recurrent.snvids(data=pon_data)
    pon_count <- pon_count %>% filter(pon_count >= opt$recurrence)

    save(count, file=output)
    }

### PROCESS DATA ##################################################################################

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################

### SESSION INFORMATION ###########################################################################
sessionInfo()
