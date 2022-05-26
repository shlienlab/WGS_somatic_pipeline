### barchart.R ####################################################################################
# Create the alignment barchart.

### HISTORY #######################################################################################
# Version           Date            Developer                     Comments
#-------------------------------------------------------------------------
# 0.01                 ?                    ?          initial development
# 0.02        2019-07-31        Drew Thompson        move to pipeline repo,
#                                                    bring in external functions

### NOTES #########################################################################################
#

### PREAMBLE ######################################################################################
library('getopt')

usage <- function() {
    usage.text <-
        '\nUsage: barchart.R --input test.idxstats --sample samplename --output out.pdf\n\n'
    return(usage.text)
    }

params = matrix(
    c(
        'input', 'i', 1, 'character',
        'sample', 's', 1, 'character',
        'output', 'o', 1, 'character'
        ),
    ncol = 4,
    byrow = TRUE
    )

opt = getopt(params)

# verify arguments, the input filename should end in .idxstats
if (is.null(opt$input)) { stop(usage()) }
if (is.null(opt$output)) {
    filename <- unlist(strsplit(x=opt$input, split='\\.'))
    filename <- paste(
        sep='.',
        paste(collapse='.', filename[1:(length(filename) - 1)]),
        'alignment.barchart.pdf'
        )
} else {
  filename <- opt$output
  }

### LIBRARIES #####################################################################################
library("ggplot2")

### FUNCTIONS #####################################################################################
write.plot <- function(filename = NULL, plot = NULL, size.units = 'in', width = 8, height = 8, resolution = 1600) {
  if(is.null(plot)) stop("Mandatory argument plot is missing")
  ggsave(
    filename = filename,
    plot = plot,
    units = size.units,
    width = width,
    height = height,
    dpi = resolution
    );
  }

create.alignment.bar.chart <- function(statsfile=NULL, filename=NULL, sample=NULL) {
  if (is.null(statsfile)) stop("Mandatory argument statsfile is missing")
  if (is.null(sample)) stop("Mandatory argument sample is missing")
  if (is.null(filename)) {
    filename <- paste(sep='.', sample, 'alignment', 'piechart', 'pdf')
    filename <- paste(sep='/', getwd(), filename)
  }
  data <- read.table(
    file=statsfile,
    header=FALSE,
    as.is=TRUE,
    sep='\t',
    quote="\""
  )
  colnames(data) <- c(
    "chr",
    "length",
    "aligned",
    "unaligned"
  )
  
  # to ensure consistency between reference genomes, remove the chr prefix to the chromosome
  data$chr <- gsub(data$chr, pattern='^chr', replacement='')
  total.reads <- sum(data$aligned, data$unaligned)
  aligned.reads <- sum(data[grep(x=data$chr, pattern='^[0-9XYM]'),]$aligned)
  random.unknown.aligned.reads <- sum(data[grep(x=data$chr, pattern='^[GN]'),]$aligned)
  
  # different references use different decoy sequences
  decoy.aligned.reads <- sum(
    data[data$chr == 'hs37d5',]$aligned,
    data[grepl(x=data$chr, pattern='decoy'),]$aligned
    )
  
  unaligned.reads <- sum(data$unaligned)
  read.data <- c(
    aligned.reads,
    random.unknown.aligned.reads,
    decoy.aligned.reads,
    unaligned.reads
  )
  
  # reference http://www.statmethods.net/graphs/pie.html
  # for adding percentages to pie chart labels
  read.data.percentages <- round(read.data / total.reads * 100)
  chart.labels <- c('Aligned', 'Random and Unknown', 'Decoy', 'Unaligned')
  chart.labels <- paste(sep=' ', chart.labels, read.data.percentages)
  chart.labels <- paste(sep='', chart.labels, '%')
  chart.data <- data.frame(
    labels=chart.labels,
    slice=read.data,
    perc=read.data.percentages
  )
  chart.data$sample <- sample
  plot.object <- ggplot(data=chart.data, aes(x=sample, y=slice, fill=labels)) +
    geom_bar(stat='identity') +
    scale_fill_brewer(palette='Dark2') +
    theme_bw(base_size=20) +
    xlab('') +
    ylab('Count')
  write.plot(
    filename=filename,
    plot=plot.object
    )
  return.values <- list(
    plot=plot.object,
    output=filename,
    data=chart.data
    )
  return(return.values)
  }

### GET DATA ######################################################################################

### PROCESS DATA ##################################################################################

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################
create.alignment.bar.chart(
    statsfile=opt$input,
    filename=filename,
    sample=opt$sample
    )

### SESSION INFORMATION ###########################################################################
sessionInfo()

