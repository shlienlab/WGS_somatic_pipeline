### filter_clips.R #############################################################################
# A R script that loads an a mutation dataframe and associated BAM file and creates a dataframe where all hard clippings are flagged

### HISTORY #######################################################################################
# Version           Date            Developer                                Comments
#------------------------------------------------------------------------------------
#    0.01     2016-12-08              chrissy                     initial development
#    0.02     2017-04-12             rdeborja    Rscript now copied to exec directory
#    0.03     2019-07-31        Drew Thompson                   move to pipeline repo,
#                                                         bring in external functions
#    0.04     2019-12-09  Lisa-Monique Edward           allow for chromosomes with no 
#                                                                 unfiltered variants


### NOTES #########################################################################################

### PREAMBLE ######################################################################################
library('getopt')

usage <- function() {
  usage.text <- '\nUsage: filter_clips.R --bam test.bam --rda test.rda --output_file output_file --offsets 10\n\n'
  return(usage.text)
}

params = matrix(
  c(
	'bam', 'b', 1, 'character',
  	'rda', 'r', 1, 'character',
  	'output_file', 's', 0, 'character',
	'offsets', 'o', 0, 'integer'
	),
  ncol = 4,
  byrow = TRUE
)

opt = getopt(params)

# verify arguments
if(is.null(opt$bam)) { stop(usage()) }
if(is.null(opt$rda)) { stop(usage()) }
if(is.null(opt$output_file)) { stop(usage()) }

### LIBRARIES #####################################################################################
library("plyr")
library("dplyr")
library("Rsamtools")
library("readr")

### FUNCTIONS #####################################################################################
get_which <- function(chr = NULL, position = NULL, offsets = NULL) {
  if(is.null(position)) stop("Mandatory argument position is missing")
  if(is.null(chr)) stop("Mandatory argument chr is missing")
  if(is.null(offsets)) stop("Mandatory argument offsets is missing")
  
  if (chr %in% c(1:22, 'X', 'Y')) {
    RL <- IRangesList(IRanges(position-offsets,position+offsets))
    names(RL) <- chr

    return(RL)
  }  
}

load_bam_reads <- function(bam_file = NULL, chr_mut = NULL, chr = NULL, position = NULL, operation = "S", offsets = NULL) {

  if(is.null(bam_file)) stop("Mandatory argument bam_file is missing")
  if(is.null(chr)) stop("Mandatory argument chr is missing")
  if(is.null(offsets)) stop("Mandatory argument offsets is missing")
  if(is.null(chr_mut) & is.null(position)) stop("Function requires one of either position or chr_mut argument")

  # If a position isn't specified, look through all chromosome mutation positions
  if (is.null(position)) {
    positions <- chr_mut$annovar_start
    whole_bam <- NULL

    for (pos in 1:length(positions)) {
      what <- c('pos', 'cigar', 'mapq', 'seq')
      which <- get_which(chr=as.character(chr), position=positions[pos], offsets=offsets)
      param <- ScanBamParam(what=what, which=which, tag="MD")
      bam <- as.data.frame(scanBam(bam_file, param = param))
      names(bam) <- c('pos', 'mapq', 'cigar', 'seq', 'md')
      whole_bam <- rbind(whole_bam, bam)
    }
  }

  # If a position is specified, do not go through for loop
  else {
    what <- c('pos', 'cigar', 'mapq', 'seq')
    which <- get_which(chr=as.character(chr), position=position, offsets=offsets)
    param <- ScanBamParam(what=what, which=which, tag="MD")
    whole_bam <- as.data.frame(scanBam(bam_file, param = param))
    names(whole_bam) <- c('pos', 'mapq', 'cigar', 'seq', 'md')
  }

  # filter only reads with soft clippings (or other operation) in CIGAR string
  s_cigars <- whole_bam %>% filter(grepl(operation, cigar))
  s_cigars$cigar <- sapply(s_cigars$cigar, as.character)
  return(s_cigars)
}

cigar_to_data_frame <- function(cigar = NULL, position = NULL)
{
  if(is.null(cigar)) stop("Mandatory argument cigar is missing")
  if(is.null(position)) stop("Mandatory argument position is missing")

  # get seperate lists of numbers (lengths) and characters (CIGAR operations)
  alpha <- unlist(strsplit(cigar, "[0-9]+"))
  alpha <- tail(alpha, length(alpha)-1)
  num <- as.numeric(unlist(strsplit(cigar, "[MIDNSHP=X]+")))

  # calculate each start and end position based on whole position and previous string
  start <- position
  if (length(alpha) > 1) {
    for (x in 2:length(alpha)) {
      start[x] <- start[x-1] + num[x-1]
    }
  }
  end <- start + num - 1

  # return entire data frame
  data <- as.data.frame(cbind(type=alpha, size=num, start=start, end=end))
  data$type <- sapply(data$type, as.character)
  return(data)
}

md_snv_position <- function(md = NULL, position = NULL, cigar = NULL, seq = NULL, operation = "S") {
  if(is.null(md)) stop("Mandatory argument md is missing")
  if(is.null(position)) stop("Mandatory argument position is missing")
  if(is.null(cigar)) stop("Mandatory argument cigar is missing")
  if(is.null(seq)) stop("Mandatory argument seq is missing")

  start_pos <- position
  pos_vector <- NULL
  alts <- NULL
  cigar_df <- cigar_to_data_frame(cigar = cigar, position = position)

  if (operation != "S") {
    # does md contain mutation, or is it a matched read (not including insertions)
    if (!grepl("[TCAG]", md)) {
      return("unmutated")
    }

    # does the md read contain ^ (check for deletions)
    if (grepl("[[:punct:]]", md)) {
      return("indel")
    }

    #do something to make sure there are not 2 mutations after one another
    alpha <- unlist(strsplit(md, "[0-9]+"))
    num <- unlist(strsplit(md, "[TCAG]+"))

    doubles <- alpha[nchar(alpha) > 1]
    for(d in doubles) {
      d_vector <- unlist(strsplit(d, ""))
      newd <- paste(d_vector, 0, sep="", collapse="")
      index <- grep(d, alpha)
      alpha[index] <- newd
    }

    while (length(alpha) > length(num)) {
      num <- c(num, "")
    }

    while (length(alpha) < length(num)) {
      alpha <- c(alpha, "")
    }

    if(alpha[1]=="") {
      newmd <- paste0(alpha, num)
    } else {
      newmd <- paste0(num, alpha)
    }
    md <- paste(newmd, collapse="")
    alpha <- unlist(strsplit(md, "[0-9]+"))
    num <- unlist(strsplit(md, "[TCAG]+"))

    # does the cigar string have the same length as the md read (check for insertion)
    clipped_df <- cigar_df[cigar_df$type != "H",]
    clipped_df <- clipped_df[clipped_df$type != "S",]
    cigar_length <- sum(as.numeric(as.character(clipped_df$size)))
    md_length <- sum(as.numeric(num), na.rm = TRUE) + length(alpha[alpha != ""])
    if (md_length != cigar_length) {
      return("indel")
    }

    # shave off the unnecessary front of alpha and num, and update starting position
    if (alpha[1] == "") {
      position <- position + as.numeric(num[1])
      num <- tail(num, length(num)-1)
      alpha <- tail(alpha, length(alpha)-1)
    } else {
      num <- tail(num, length(num)-1)
    }

    # calculate the position of each of the mutations
    for (mut in alpha) {
      pos_vector <- c(pos_vector, position)
      alts <- c(alts, substr(seq, position-start_pos+1, position-start_pos+1))
      position <- position + as.numeric(num[1]) + 1
      num <- tail(num, length(num)-1)
    }
  }

  # for soft clips, add every mutation within the soft clip
  else {
    # what happens with indels?

    # get all bases in soft clip at the beginning of read
    if (cigar_df[1,"type"]=="S") {
      clip_len <- cigar_df[1,"size"]
      for (base_pos in 1:as.numeric(as.character(clip_len))) {
        alts <- c(alts, substr(seq, base_pos, base_pos))
        pos_vector <- c(pos_vector, start_pos + base_pos - 1)
      }
    }

    # get all bases in soft clip at the end of read
    if (cigar_df[nrow(cigar_df),"type"]=="S") {
      clip_len <- cigar_df[nrow(cigar_df),"size"]
      clip_start <- sum(as.numeric(as.character(cigar_df$size))) - as.numeric(as.character(clip_len)) + 1
      for (base_pos in clip_start:(clip_start + as.numeric(as.character(clip_len)))) {
        alts <- c(alts, substr(seq, base_pos, base_pos))
        pos_vector <- c(pos_vector, start_pos + base_pos - 1)
      }
    }
  }

  #return a data frame of the
  mutation_df <- as.data.frame(cbind(alt=alts, position=pos_vector))
  return(mutation_df)
}

get_cigar_type <- function(operation = NULL) {
  if(is.null(operation)) stop("Mandatory argument 'operation' is missing")

  if (operation=="M") return("alignment_matches")
  if (operation=="I") return("insertions")
  if (operation=="D") return("deletions")
  if (operation=="N") return("skipped_regions")
  if (operation=="S") return("soft_clips")
  if (operation=="H") return("hard_clips")
  if (operation=="P") return("padding")
  if (operation=="=") return("sequence_matches")
  if (operation=="X") return("sequence_mismatches")
}

get_cigar_mutations <- function(bam_file=NULL, mut_file=NULL, chr=NULL, position=NULL, operation="S", offsets=10) {

  if(is.null(bam_file)) stop("Mandatory argument bam_file is missing")
  if(is.null(mut_file)) stop("Mandatory argument mut_file is missing")
  if(is.null(chr)) stop("Mandatory argument chr is missing")

  #load mutations file if file directory or data frame
  if (is.data.frame(mut_file)) {
    muts <- mut_file
  } else {
    snv <- load(mut_file)
    muts <- get(snv)
  }

  # update offsets to be non zero only if soft clips are called (clipped at the front)
  if (operation == "S") {
    offsets = offsets
  } else {
    offsets = 0
  }

  # filter only needed chromosome if chr included, change position to useable numbers
  chr_mut <- muts %>% filter(annovar_chr == as.character(chr))
  chr_mut$annovar_start <- as.numeric(as.character(chr_mut$annovar_start))

  # load bam file.
  s_reads <- load_bam_reads(bam_file=bam_file, chr_mut=chr_mut, chr=chr, position=position, operation=operation, offsets=offsets)

  # get unique soft clipping positions from reads
  s_clips_pos <- NULL
  if (nrow(s_reads) > 0) {
    for (read in 1:nrow(s_reads)) {
      s_frame <- cigar_to_data_frame(cigar=s_reads$cigar[read], position=s_reads$pos[read])

      # if the operation is hard clipping, mark the entire cigar read as hard clipped
      if (operation=="H") {
        H_read_pos <- data.frame(start=s_frame[1,"start"], end=s_frame[nrow(s_frame), "end"])
        s_clips_pos <- rbind(s_clips_pos, H_read_pos)
      } else {
        if (operation=="S" & s_frame[1,"type"]=="S") {
          s_reads[read,"pos"] <- s_reads[read,"pos"] - as.numeric(as.character(s_frame[1,"size"]))
          s_frame <- cigar_to_data_frame(cigar=s_reads$cigar[read], position=s_reads$pos[read])
        }
        s_frame <- s_frame %>% filter(type==operation)
        s_clips_pos <- rbind(s_clips_pos, s_frame[,c('start', 'end')])
      }
    }
    s_clips_pos$start <- as.numeric(as.character(s_clips_pos$start))
    s_clips_pos$end <- as.numeric(as.character(s_clips_pos$end))
  }

  # select only the mutations in the chromosome which occur in soft clippings
  f <- function(pos) any(ifelse(s_clips_pos$start<=pos & pos<=s_clips_pos$end, TRUE, FALSE))
  chr_mut$in_soft_clips <- sapply(chr_mut$annovar_start, f)

  chr_mut$index <- 1:nrow(chr_mut)

  # go through each positive, check where the mutation occurs (from md_snv_position)
  # and compare against the given position, to see if they are at the same position
  s_muts <- chr_mut %>% filter(in_soft_clips == TRUE)

  if (nrow(s_muts) > 0) {
    mut_in_hard_read <- 0

    # go through each read to see if mutation occurs in correct location
    for (read in 1:nrow(s_reads)) {
      mutation_pos <- md_snv_position(md=as.character(s_reads[read, "md"]),
                                      position=as.numeric(as.character(s_reads[read, "pos"])),
                                      cigar=as.character(s_reads[read, "cigar"]),
                                      seq=as.character(s_reads[read, "seq"]),
                                      operation=operation)

      if (is.data.frame(mutation_pos)) { #if position in right place, keep true
        if (nrow(merge(data.frame("alt"=as.character(s_muts[1,"annovar_alt"]), "position"=position), mutation_pos,
                       by=c("alt", "position"))) > 0) {
          mut_in_hard_read <- mut_in_hard_read + 1
        }
      } else if (mutation_pos == 'indel') { #indel, keep as true, just in case
        mut_in_hard_read <- mut_in_hard_read + 1
      }
    }
    chr_mut$in_soft_clips[s_muts$index[1]] <- mut_in_hard_read
  }

  # append column to initial data frame
  cigar_type <- paste("in", get_cigar_type(operation), sep="_")

  muts$in_S <- ifelse(muts[,"annovar_chr"]==chr, chr_mut$in_soft_clips, 0)

  # rename in_S column
  names(muts) <- c(head(names(muts), -1), eval(parse(text="cigar_type")))
  return(muts)
}

count_rda_clips <- function(bam = NULL, rda = NULL, operation = "S", offsets = 10) {
  if(is.null(bam)) stop("Mandatory argument bam is missing")
  if(is.null(rda)) stop("Mandatory argument rda is missing")

  if (is.data.frame(rda)) {
    muts <- rda
  } else {
    muts <- read_tsv(rda)
  }

  clips_df <- NULL
  if (nrow(muts) > 0 ) {
    for (row in 1:nrow(muts)) {
      new_row <- get_cigar_mutations(bam_file=bam, mut_file=muts[row,], chr=unlist(muts[row,"annovar_chr"], use.names=F),
                                     position=unlist(muts[row,"annovar_start"], use.names=F), operation=operation, offsets=offsets)
      clips_df <- rbind(clips_df, new_row)
    }
  } else {
    clips_df <- data.frame(matrix(ncol=(ncol(muts) + 1), nrow=0))
    if (operation == "S") {
      colnames(clips_df) <- c(colnames(muts), "in_soft_clips")
    } else {
      colnames(clips_df) <- c(colnames(muts), "in_hard_clips")
    }
    
  }
  
  return(clips_df)
}

filter_clips <- function(rda = NULL, bam = NULL) {
  hard.filtered <- count_rda_clips(bam = bam, rda = rda, operation = "H")
  soft.filtered <- count_rda_clips(bam = bam, rda = hard.filtered, offsets = 30)
  return(soft.filtered)
}

### GET DATA ######################################################################################

### PROCESS DATA ##################################################################################

### ANALYSIS ######################################################################################
data.filtered <- filter_clips(rda=opt$rda, bam=opt$bam)
save(x=data.filtered, file=opt$output_file)

### PLOTTING ######################################################################################

### SESSION INFORMATION ###########################################################################
sessionInfo()
