### run_bicseq.R ####################################################################################
# Run the BICSeq pipeline at the command line using Rscript.

### HISTORY #######################################################################################
# Version           Date            Developer                     Comments
#-------------------------------------------------------------------------
# 0.01        2014-07-08             rdeborja          initial development
# 0.02        2019-07-31        Drew Thompson        move to pipeline repo,
#                                                    bring in external functions


### PREAMBLE ######################################################################################
library('getopt')

usage <- function() {
    usage.text <- '\nUsage: script.R --normal normal.bam --tumour tumour.bam\n\n';
    return(usage.text)
    }

params = matrix(
    c(
    	'tumour', 't', 1, 'character',
    	'normal', 'n', 1, 'character'
    	),
    ncol = 4,
    byrow = TRUE
    )

opt = getopt(params)

# verify arguments
if(is.null(opt)) { stop(usage()) }

### LIBRARIES #####################################################################################
library("BICseq")

### FUNCTIONS #####################################################################################
get.bicseq.chr <- function(reftype='b37') {
  valid_reftype = c('b37', 'ucsc')
  if (! reftype %in% valid_reftype) stop(paste(sep=' ', reftype, 'is an invalid reftype'))
  if (reftype == 'ucsc') {
    return(
      c(
        'chr1',
        'chr2',
        'chr3',
        'chr4',
        'chr5',
        'chr6',
        'chr7',
        'chr8',
        'chr9',
        'chr10',
        'chr11',
        'chr12',
        'chr13',
        'chr14',
        'chr15',
        'chr16',
        'chr17',
        'chr18',
        'chr19',
        'chr20',
        'chr21',
        'chr22',
        'chrX',
        'chrY'
        )
      )
   } else if (reftype == 'b37') {
    return(
      c(
        '1',
        '2',
        '3',
        '4',
        '5',
        '6',
        '7',
        '8',
        '9',
        '10',
        '11',
        '12',
        '13',
        '14',
        '15',
        '16',
        '17',
        '18',
        '19',
        '20',
        '21',
        '22',
        'X',
        'Y'
        )
      )
    }
  }

run.bicseq.cnv.pipeline <- function(tumour=NULL, normal=NULL, chr = get.bicseq.chr(), bin = 100, lambda = 2, winSize = 200, quant = 0.95, mult = 1) {
  if (is.null(tumour)) stop("Mandatory argument tumour is missing")
  if (is.null(normal)) stop("Mandatory argument normal is missing")
  
  bicseq <- list()
  bicseq$object <- BICseq(
    sample = tumour,
    reference = normal,
    seqNames = chr
    );
  segs <- BICseq::getBICseg(
    object = bicseq$object,
    bin = bin,
    lambda = lambda,
    winSize = winSize,
    quant = quant,
    mult = mult
    );
  bicseq$bins <- BICseq:::getRatios(bin(segs), what = 'bin');
  bicseq$segs <- segs;
  bicseq$seg.summary <- BICseq:::getSummary(segs, correction=TRUE);
  bicseq$seg.summary$seg_length <- bicseq$seg.summary$end - bicseq$seg.summary$start + 1
  return(bicseq);
  }

### GET DATA ######################################################################################

### PROCESS DATA ##################################################################################

bicseq <- run.bicseq.cnv.pipeline(
	normal = opt$normal,
	tumour = opt$tumour
	)
save(bicseq, file = 'bicseq_data.rda')

### ANALYSIS ######################################################################################

### PLOTTING ######################################################################################

