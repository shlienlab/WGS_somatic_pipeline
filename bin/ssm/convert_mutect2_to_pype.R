### convert_mutect_to_pype.R ######################################################################
# A R script that converts MuTect annotated R dataframes to a tab separated text file
# this is Pype compatible.

### HISTORY #######################################################################################
# Version           Date            Developer                     Comments
#-------------------------------------------------------------------------
# 0.01        2016-09-27             rdeborja          initial development
# 0.02        2019-07-31        Drew Thompson        move to pipeline repo,
#                                                    bring in external functions
# 0.03        2019-12-09  Lisa-Monique Edward    allow for chromosomes with 
#                                                    no unfiltered variants

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
library(devtools)

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

get.mutect.pype.header <- function() {
  return(
    c(
      "Chromosome1",
      "Chromosome2",
      "Position1",
      "Position2",
      "Ref.Allele",
      "Alt.Allele",
      "Func","Gene",
      "VAF","mutation.type",
      "trinuc",
      "n_ref_count",
      "n_alt_count",
      "t_ref_count",
      "t_alt_count"
      )
    )
  }

### GET DATA ######################################################################################
somatic_data <- load(opt$file)
somatic_data <- get(somatic_data)

### PROCESS DATA ##################################################################################
if (nrow(somatic_data) > 0 ) {
	somatic_data$trinuc <- 'none'
	somatic_data$mutation.type <- 'none'
	pype.mutect.data <- reorder.mutect2.data.to.pype(data=somatic_data)
} else {
	pype.mutect.data <- data.frame(matrix(ncol=length(get.mutect.pype.header()), nrow=0))
}
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
