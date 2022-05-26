# Version        Date              Developer                                         Comments
#--------------------------------------------------------------------------------------------
#     0.0        ?                 ?                                      initial development
#     0.1        2018-12-17        Drew Thompson                              updated for WDL
#   0.1.1        2019-05-21        Drew Thompson         skip input files if they don't exist

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Must provide the merged output file as the first argument, and then all final mutation type tab files")
}

merged_data = NULL
for (i in 2:length(args)) {
	this_file = args[i]
	if (file.exists(this_file)) {
		data = read.table(this_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
		merged_data = rbind(merged_data, data)
	}
}
write.table(merged_data, file = args[1], sep = "\t", row.names = FALSE, quote = FALSE)
