### convert_mutect_to_pype.R ######################################################################
# A R script that converts MuTect annotated R dataframes to a tab separated text file
# this is Pype compatible.

### HISTORY #######################################################################################
# Version           Date            Developer               Comments
# 0.01              2016-09-27      rdeborja                initial development

### NOTES #########################################################################################
#

### PREAMBLE ######################################################################################
library('getopt')

usage <- function() {
    usage.text <- '\nUsage: convert_mutect_to_pype.R --file test.rda --output output.txt\n\n'
    return(usage.text)
    }

params = matrix(
    c('file', 'f', 1, 'character', 'output', 'o', 1, 'character'),
    ncol = 4,
    byrow = TRUE
    )

opt = getopt(params)

# verify arguments
if(is.null(opt$file)) { stop(usage()) }

### LIBRARIES #####################################################################################
library(ShlienLab.Core.SNV)
library(devtools)
#load_all(pkg='/home/rdeborja/ShlienLab.Core.SNV')

### FUNCTIONS #####################################################################################
reorder.mutect2.data.to.pype <- function(data=NULL) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  return(
    data[
      c(
        "annovar_chr",
        "annovar_chr",
        "annovar_start",
        "annovar_end",
        "annovar_ref",
        "annovar_alt",
        "annovar_func",
        "annovar_gene",
        "gatk_tumour_allele_fraction",
        "mutation.type",
        "trinuc",
        "gatk_normal_ref_count",
        "gatk_normal_alt_count",
        "gatk_tumour_ref_count",
        "gatk_tumour_alt_count"
        )
      ]
    )
  }

### GET DATA ######################################################################################
somatic_data <- load(opt$file)
somatic_data <- get(somatic_data)

### PROCESS DATA ##################################################################################
somatic_data$trinuc <- 'none'
somatic_data$mutation.type <- 'none'
pype.mutect.data <- reorder.mutect2.data.to.pype(data=somatic_data)
colnames(pype.mutect.data) <- get.mutect.pype.header()
write.table(
    x=pype.mutect.data,
    file=opt$output,
    sep='\t',
    quote=FALSE,
    row.names=FALSE
    )

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################

### SESSION INFORMATION ###########################################################################
sessionInfo()
