# Version              Date            Developer                                Comments
#---------------------------------------------------------------------------------------
#     0.0                 ?                    ?                     initial development
#     0.1        2018-12-17        Drew Thompson                         updated for WDL
#   0.1.1        2019-07-31        Drew Thompson        directly pass in output filename
#   0.1.2        2019-08-22        Drew Thompson              Save filter info for graph

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
    stop("Must provide the pairToPair result, the almost_filtered file, the output filename, the filter info file, and the SV type.")
}

pair_to_pair_result = args[1]
almost_filtered_file = args[2]
output_file = args[3]
filter_info = args[4]
sv_type = args[5]

variant_count = NULL
if (file.info(pair_to_pair_result)$size == 0) {
    sample.sv.dataset = read.table(file=almost_filtered_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    variant_count = nrow(sample.sv.dataset)
    file.copy(almost_filtered_file, output_file)
} else if (file.exists(almost_filtered_file)) {
    to.be.removed <- read.table(file=pair_to_pair_result, header=FALSE, sep="\t")
    to.be.removed$key <- paste(sep=".", to.be.removed$V1, to.be.removed$V2, to.be.removed$V3, to.be.removed$V4, to.be.removed$V5, to.be.removed$V6, to.be.removed$V7)
    sample.sv.dataset = read.table(file=almost_filtered_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    completely.filtered <- sample.sv.dataset[!sample.sv.dataset$key %in% to.be.removed$key,]
    variant_count = nrow(completely.filtered)
    write.table(completely.filtered, file=output_file, sep="\t", row.names=FALSE, quote=FALSE)
} else {
    variant_count = 0
}

cat(sv_type, "\tPON\t", variant_count, "\n", sep = "", file = filter_info, append = TRUE)