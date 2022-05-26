### run_ssm_standard_filters.R ####################################################################
# Run the standard simple somatic filters on MuTect2 output.  Output will be separate for SNVs and
# indels.

### HISTORY #######################################################################################
# Version           Date            Developer               Comments
# 0.01              2017-04-10      rdeborja                initial development
# 0.02              2017-04-13      rdeborja                removed MT from dataframe due to
#                                                           clip filtering issues

### NOTES #########################################################################################
#

### PREAMBLE ######################################################################################
library('getopt')

usage <- function() {
    usage.text <- '\nUsage: run_ssm_standard_filters.R --path </path/to/directory/containing/files> --sample <sample name> --source <WGS|WXS|CPANEL>\n\n'
    return(usage.text)
    }

params = matrix(
    c(
        'path', 'p', 1, 'character',
        'sample', 's', 1, 'character',
        'source', 'c', 1, 'character'
        ),
    ncol = 4,
    byrow = TRUE
    )

opt = getopt(params)

# verify arguments
if (is.null(opt$path)) { stop(usage()) }

output <- paste(sep='.', paste(sep='_', opt$sample, 'annotated'), 'rda')
snv_filtered <- paste(sep='.', paste(sep='_', opt$sample, 'annotated_filtered_snv'), 'rda')
indel_filtered <- paste(sep='.', paste(sep='_', opt$sample, 'annotated_filtered_indel'), 'rda')

### LIBRARIES #####################################################################################
library(ShlienLab.Core.SSM)

### FUNCTIONS #####################################################################################

### GET DATA ######################################################################################
data <- get.mutect2.data(path=opt$path)

### PROCESS DATA ##################################################################################
# add additional annotations to the dataframe
data <- ShlienLab.Core.SSM::annotate.mutect2.data(data=data)

# currently there is a bug in downstream filtering causing a pre-filter step to remove the
# mitochondrial DNA from the output
data <- data %>% filter(annovar_chr != 'MT')
data <- data %>% filter(annovar_chr != 'M')

# separately filter the snv and indel data
data.snv.filtered <- ShlienLab.Core.SSM::filter_snv(data=data, source='WGS')
data.indel.filtered <- ShlienLab.Core.SSM::filter_indel(data=data, source='WGS')

# save the dataframes
save(data, file=output)
save(data.snv.filtered, file=snv_filtered)
save(data.indel.filtered, file=indel_filtered)

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################

### SESSION INFORMATION ###########################################################################
sessionInfo()

