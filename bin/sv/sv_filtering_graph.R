# Version              Date            Developer                          Comments
#---------------------------------------------------------------------------------
#     0.1        2019-08-22        Drew Thompson               initial development

library(ggplot2)
library(plyr)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3) { stop("Please provide the sample id, the name of the output graph, and one or more filter info tsv files") }

sample = args[1]
output_file = args[2]

filter_info = NULL
for (i in 3:length(args)) {
    tmp_info_file = args[i]
    tmp_info = read.table(tmp_info_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(tmp_info) = c("Type", "Stage", "Variants")
    filter_info = rbind(filter_info, tmp_info)
}

filter_info$Stage = factor(filter_info$Stage, levels = unique(filter_info$Stage))

plot = ggplot(data = filter_info, aes(x = reorder(Stage, desc(Stage)), y = Variants, fill = Type)) +
    geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
    coord_flip() +
    geom_text(aes(label = Variants), position=position_dodge(width = 0.9), hjust = "inward") +
    labs(title = sample, x = "Filtering Stage", y = "Number of Variants Remaining", fill = "Variant Type") +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_legend(reverse=TRUE))

ggsave(output_file, plot = plot)