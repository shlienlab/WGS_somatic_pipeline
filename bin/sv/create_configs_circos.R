# Version              Date            Developer                              Comments
#-------------------------------------------------------------------------------------
#     0.0                 ?                    ?                   initial development
#     0.1        2018-12-17        Drew Thompson                       updated for WDL
#   0.1.1        2019-07-31        Drew Thompson        provide output files as output
#   0.1.2        2019-08-07        Drew Thompson       don't fail if no rearrangements,
#                                                      remove hardcoded values

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 10) {
    stop("Must provide the final merged tab file, the sample id, the circos_all output file, the output config file, and config files")
}

final_file = args[1]
sample_id = args[2]
output_circos_all_file = args[3]
output_config_tab_file = args[4]
karyotype = args[5]
ideogram.conf = args[6]
image.conf = args[7]
colors_fonts_patterns.conf = args[8]
housekeeping.conf = args[9]
ticks.conf = args[10]

create.coordinates <- function(data, sample, output_file) {
  for(x in 1:nrow(data)) {
    print(paste(sep="_", unique(as.character(data$Sample)), x))
    split1 <- as.data.frame(strsplit(data$Chromosome1[x], fixed=T, split="chr"))
    data$chr1[x] <- as.character(split1[2,])
    split2 <- as.data.frame(strsplit(data$Chromosome2[x], fixed=T, split="chr"))
    data$chr2[x] <- as.character(split2[2,])
  }
  
  data$cchrom1 <- paste(sep="", "hs", data$chr1)
  data$cchrom2 <- paste(sep="", "hs", data$chr2)
  data$cpos1_1 <- data$Position1
  data$cpos1_2 <- data$Position1
  data$cpos2_1 <- data$Position2
  data$cpos2_2 <- data$Position2
  sample.of.interest <- data[data$Sample == sample,]
  sample.of.interest$Color <- "color=black,thickness=5p"
  
  circos.coordinates.all <- sample.of.interest[,c("cchrom1","cpos1_1","cpos1_2", "cchrom2","cpos2_1","cpos2_2","Color")]
  
  write.table(circos.coordinates.all, file=output_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

all.rearrangements <- read.table(file=final_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)

if (nrow(all.rearrangements) == 0) {
  print("No rearrangements")
  quit(save = "no", status = 0)
}

create.coordinates(all.rearrangements, sample_id, output_circos_all_file)

info1 <- paste0("
karyotype = ", karyotype, "

<ideogram>
  show_label = yes
label_font = bold
label_size = 25
label_parallel = no
label_case = upper
label_with_tag = yes
label_radius = 1r + 10p

<spacing>
default = 0.005r
</spacing>
  
  radius=0.9r
thickness = 20p
fill = yes
stroke_thickness = 1
stroke_color = black

<<include ", ideogram.conf, ">>
  
</ideogram>
  
###############################################
  
<image>
<<include ", image.conf, ">>
 radius* = 500
</image>
  
<<include ", colors_fonts_patterns.conf, ">>
  
<<include ", housekeeping.conf, ">>
  
<ticks>
<<include ", ticks.conf, ">>
</ticks>
  
<links>
  
<link>
file = ", basename(output_circos_all_file), " 
color = black_a5
radius=0.95r
thickness = 5
bezier_radius = 0.1r
</link>
  
</links>
")

write.table(info1, file=output_config_tab_file, sep="\t", col.names=F, row.names=F, quote=F)
