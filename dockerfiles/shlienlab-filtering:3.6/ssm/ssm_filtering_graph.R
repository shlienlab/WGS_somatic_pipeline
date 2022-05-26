# Version              Date            Developer                          Comments
#---------------------------------------------------------------------------------
#     0.1        2019-08-22        Drew Thompson               initial development

library(ggplot2)
library(plyr)

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 3) { stop("Please provide the sample id, filter info tsv file, and the name of the output graph.") }

sample = args[1]
filter_info_file = args[2]
output_file = args[3]

filter_info = read.table(filter_info_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

make_filter_graph <- function(sample, filter_info, output_file) {
    colnames(filter_info) = c("Variant.Type", "Stage", "Variants")
    filter_info$Stage = paste(sep = ": ", filter_info$Variant.Type, filter_info$Stage)
    filter_info$Stage = factor(filter_info$Stage, levels = filter_info$Stage)
    plot = ggplot(data = filter_info, aes(x = reorder(Stage, desc(Stage)), y = Variants, fill = Variant.Type)) +
            geom_bar(position = "dodge", stat = "identity") +
            coord_flip() +
            geom_text(aes(label = Variants), hjust = "inward") +
            labs(title = sample, x = "Filtering Stage", y = "Number of Variants Remaining", fill = "Variant Types") +
            theme(plot.title = element_text(hjust = 0.5))
    ggsave(output_file, plot = plot)
}

make_filter_graph(sample, filter_info, output_file)
