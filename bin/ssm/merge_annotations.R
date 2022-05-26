### merge_annotations.R ####################################################################
# Gather filtering results across chromosomes.  Output will be separate for SNVs and
# indels.

### HISTORY #######################################################################################
# Version           Date            Developer                                Comments
#------------------------------------------------------------------------------------
#    0.01     2019-12-09      Lisa-Monique Edward     separate filtering and merging steps from
#                                                     cromwell_run_standard_filters.R
#    0.02     2020-10-01      Lisa-Monique Edward     Store and read txts vs rdas
#    0.03     2020-11-23      Lisa-Monique Edward         Store filtered variants for
#                                                                          annotation

### NOTES #########################################################################################
#

### PREAMBLE ######################################################################################
usage <- function() {
    usage.text <- '\nUsage: merge_annotations.R <snv output> <indel output> <filter info output> <target positions output> <snvFile1,snvFile2,...,snvFileN> <indelFile1,indelFile2,...,indelFileN> <filterFile1,filterFile2,...,filterFileN>\n\n'
    return(usage.text)
    }

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 9) stop(usage())
snv_output = args[1]
indel_output = args[2]
target_tsv = args[3]
filter_info = args[4]
annotations = args[5]
snvs = strsplit(args[6], ",")[[1]]
indels = strsplit(args[7], ",")[[1]]
filter_files =  strsplit(args[8], ",")[[1]]
annot_files = strsplit(args[9], ",")[[1]]

### LIBRARIES #####################################################################################
library("plyr")
library("dplyr")
library("stringr")
library("readr")

### FUNCTIONS #####################################################################################
merge.txt <- function(files) {
  txt_data <- data.frame()
  
  for (f in files) {
    data <- try(
      read_tsv(file = f, col_names = TRUE, col_types = cols(.default = "?", annovar_chr = col_factor(levels = intervals), annovar_ref = "c", annovar_alt = "c", annovar_esp = "d", gt_REF = "c", gt_ALT = "c", in_centromere = "c", normal_coverage_threshold = "c", unique_mapping = "c", high_depth = "c", dist_low_complexity = "c", multi_mapping = "c")),
      silent = TRUE
    )
    # the "try" block will return a class of "try-error"
    if (class(data) == 'try-error')
      next()
    if (length(txt_data) == 0)
      txt_data <- data
    else
      txt_data <- rbind(txt_data, data)
  }
  
  return(txt_data)
}


merge.rda <- function(rda_files) {
  rda_data <- data.frame()
  for (rda_file in rda_files) {
    env <- new.env()
    rf <- try(
      load(rda_file, env)[1]
    )
    if (class(env[[rf]]) == 'try-error')
      next()
    rda_data <- rbind(rda_data, env[[rf]])
  }

  return(rda_data)
}

merge.filter_info <- function(filter_info_files) {
  filter_data <- data.frame()
  for (filter_file in filter_info_files) {
      filter_file_data <- try(
        read.table(file = filter_file, header = FALSE, as.is = TRUE, sep = '\t', quote = "\""),silent = TRUE
      )
      # the "try" block will return a class of "try-error"
      if (class(filter_file_data) == 'try-error')
        next()
      filter_data <- rbind(filter_data, filter_file_data)
  }
  filter_data <- aggregate(V3~V1+V2, filter_data, sum) %>% 
    mutate(V1 = factor(V1, levels = c("SNV+INDEL", "INDEL", "SNV"))) %>% 
    arrange(V1, desc(V3))
  
  return(filter_data)
}

merge.annotations <- function(annot_files, intervals) {
  annots <- data.frame()
  for (annot_file in annot_files) {
    annot_file_data <- try(
      read.table(file = annot_file, header = TRUE, as.is = TRUE, sep = '\t', quote = "\""),silent = TRUE
    )
    if (class(annot_file_data) == 'try-error')
        next()
      annots <- rbind(annots, annot_file_data)
  }
  annots <- annots %>% 
    mutate(annovar_chr = factor(annovar_chr, levels = intervals), gt_start = as.integer(gt_start), gt_end = as.integer(gt_end)) %>% 
    arrange(annovar_chr, gt_start)
  return(annots)
}


### GET DATA ######################################################################################
intervals <- c(seq(1, 22), "X", "Y")
snv_data <- merge.txt(snvs) %>% mutate(gt_POS = as.numeric(gt_POS))
indel_data <- merge.txt(indels) %>% mutate(gt_POS = as.numeric(gt_POS))
filter_data <- merge.filter_info(filter_files)
annot_data <- merge.annotations(annot_files, intervals = intervals)


## GET TARGETS ######################################################################################
targets <- bind_rows(select(snv_data, annovar_chr, gt_POS), select(indel_data, annovar_chr, gt_POS))
targets <- as_tibble(targets) %>% 
  rename(chrom = annovar_chr, pos = gt_POS) %>%
  mutate(chrom = factor(chrom, levels = intervals), pos = as.integer(pos)) %>%
  arrange(chrom, pos)


### WRITE OUTPUT ######################################################################################

write.table(snv_data, file = snv_output, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(indel_data, file = indel_output, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(filter_data, file = filter_info, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(targets, file = target_tsv, row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(annot_data, file = annotations, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

### SESSION INFORMATION ###########################################################################
sessionInfo()

