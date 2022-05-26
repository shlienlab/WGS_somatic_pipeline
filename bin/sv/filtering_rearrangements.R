# Version              Date            Developer                                                     Comments
#------------------------------------------------------------------------------------------------------------
#     0.0                 ?                    ?                                          initial development
#     0.1        2018-12-17        Drew Thompson                                              updated for WDL
#   0.1.1        2019-05-21        Drew Thompson              pass in sample instead of getting from filename
#   0.1.2        2019-06-24        Drew Thompson        check if variants remain before low complexity filter
#   0.1.3        2019-07-31        Drew Thompson                                     pass in output filenames
#   0.1.4        2019-08-22        Drew Thompson                                   Save filter info for graph

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
    stop("Must provide the unfiltered tab file, the mutation type, the sample name, the name for the output tab file, the name for the output bed file, the name of the filtering info file.")
}

unfiltered_file = args[1]
mutation_type = args[2]
sample = args[3]
output_tab = args[4]
output_bed = args[5]
filter_info = args[6]

apply.filter <- function(data, filter_info, sv_type) {
    somatic.filter <- data[data$BedFilter != "FLAGGED" & data$delly_passed_sFilter != "False",]
    cat(sv_type, "\tsomatic_filter\t", nrow(somatic.filter), "\n", sep = "", file = filter_info, append = TRUE)

    coverage <- somatic.filter[somatic.filter$normal_coverage_threshold != "FLAGGED",]
    cat(sv_type, "\tnormal_coverage_threshold\t", nrow(coverage), "\n", sep = "", file = filter_info, append = TRUE)

    tier1 <- coverage[coverage$delly_tier == "TIER1",]
    cat(sv_type, "\tdelly_tier1\t", nrow(tier1), "\n", sep = "", file = filter_info, append = TRUE)

    mapq <- tier1[tier1$MAPQ > 30,]
    cat(sv_type, "\tMAPQ\t", nrow(mapq), "\n", sep = "", file = filter_info, append = TRUE)

    mapq$distance_to_low_complexity_1 <- as.character(mapq$distance_to_low_complexity_1)
    mapq$distance_to_low_complexity_2 <- as.character(mapq$distance_to_low_complexity_2)

    if (nrow(mapq) == 0) {
        cat(sv_type, "\tduplicates\t0\n", sep = "", file = filter_info, append = TRUE)
        cat(sv_type, "\tcFilter\t0\n", sep = "", file = filter_info, append = TRUE)
        mapq$key <- character()
        mapq$duplicated <- logical()
        return(mapq)
    }

    for (i in 1:nrow(mapq)) {
        if(is.na(mapq$distance_to_low_complexity_1[i]) == "TRUE") {
            mapq$distance_to_low_complexity_1[i] <- "FAR"
        }
        else if(mapq$distance_to_low_complexity_1[i] == "0") {
            mapq$distance_to_low_complexity_1[i] <- "OVERLAP"
        }
        if(is.na(mapq$distance_to_low_complexity_2[i]) == "TRUE") {
            mapq$distance_to_low_complexity_2[i] <- "FAR"
        }
        else if(mapq$distance_to_low_complexity_2[i] == "0") {
            mapq$distance_to_low_complexity_2[i] <- "OVERLAP"
        }
    }

    mapq$key <- paste(sep=".", mapq$Chromosome1, mapq$Chromosome2, mapq$Position1, mapq$Position2, mapq$CT, mapq$Sample)
    mapq$duplicated <- duplicated(mapq$key)
    remove.dups <- mapq[mapq$duplicated == "FALSE",]

    cat(sv_type, "\tduplicates\t", nrow(remove.dups), "\n", sep = "", file = filter_info, append = TRUE)

#flagged by any two
    UM.HD <- remove.dups[remove.dups$unique_mapping == "FLAGGED" & remove.dups$high_depth == "FLAGGED",]
    UM.LCR <- remove.dups[remove.dups$unique_mapping =="FLAGGED" & remove.dups$distance_to_low_complexity_1 == "OVERLAP" & remove.dups$distance_to_low_complexity_2 == "OVERLAP",]
    UM.MM <- remove.dups[remove.dups$unique_mapping == "FLAGGED" & remove.dups$multi_mapping == "FLAGGED",]
    HD.MM <- remove.dups[remove.dups$high_depth == "FLAGGED" & remove.dups$multi_mapping == "FLAGGED",]
    LCR.MM <- remove.dups[remove.dups$distance_to_low_complexity_1 == "OVERLAP" & remove.dups$distance_to_low_complexity_2 == "OVERLAP" & remove.dups$multi_mapping == "FLAGGED",]
    HD.LCR <- remove.dups[remove.dups$high_depth == "FLAGGED" & remove.dups$distance_to_low_complexity_1 == "OVERLAP" & remove.dups$distance_to_low_complexity_2 == "OVERLAP",]

    remove.flagged.calls <- remove.dups[!remove.dups$key %in% c(UM.HD$key,UM.LCR$key,UM.MM$key, HD.MM$key, LCR.MM$key, HD.LCR$key),]

    cat(sv_type, "\tcFilter\t", nrow(remove.flagged.calls), "\n", sep = "", file = filter_info, append = TRUE)

    return(remove.flagged.calls)
}

## Tumour Bedfiles from filtered output - for concatenated normals
create.tumour.bedfiles <- function(filtered.tumour, output_tab, output_bed) {
    if (nrow(filtered.tumour) > 0) {
        filtered.tumour$Chromosome1 <- as.character(filtered.tumour$Chromosome1)
        filtered.tumour$Chromosome2 <- as.character(filtered.tumour$Chromosome2)
        for (i in 1:nrow(filtered.tumour)) {
            filtered.tumour$Chromosome1[i] <- paste(sep="", "chr", filtered.tumour$Chromosome1[i])
            filtered.tumour$Chromosome2[i] <- paste(sep="", "chr", filtered.tumour$Chromosome2[i])
        }
        Chr1Pos1 <- filtered.tumour[,c("Chromosome1","Position1","Breakpoint1","Sample")]
        Chr2Pos2 <- filtered.tumour[,c("Chromosome2","Position2","Breakpoint2","Sample")]

        bed.window1 <- function(window) {
            for (i in 1:nrow(Chr1Pos1)) {
                Chr1Pos1$start[i] <- Chr1Pos1$Position1[i] - window
                Chr1Pos1$end[i] <- Chr1Pos1$Position1[i] + window
            }
            return(Chr1Pos1)
        }

        bed.window2 <- function(window) {
            for (i in 1:nrow(Chr2Pos2)) {
                Chr2Pos2$start[i] <- Chr2Pos2$Position2[i] - window
                Chr2Pos2$end[i] <- Chr2Pos2$Position2[i] + window
            }
            return(Chr2Pos2)
        }
        almost.bed1 <- bed.window1(500)
        almost.bed2 <- bed.window2(500)
        bedfile1 <- almost.bed1[,c("Chromosome1","start","end")]
        bedfile2 <- almost.bed2[,c("Chromosome2","start","end","Sample")]
        colnames(bedfile1) <- c("chrom1","start1","end1")
        colnames(bedfile2) <- c("chrom2","start2","end2","Sample")
        combined.tumour <- cbind(bedfile1, bedfile2) 
        filtered.tumour$key <- paste(sep=".", combined.tumour$chrom1, combined.tumour$start1, combined.tumour$end1, combined.tumour$chrom2, combined.tumour$start2, combined.tumour$end2, combined.tumour$Sample)

        write.table(combined.tumour, file = output_bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
    write.table(filtered.tumour, file = output_tab, sep="\t", quote=FALSE, row.names=FALSE)
}

dataset <- read.table(file=unfiltered_file, header=TRUE, sep="\t", fill=TRUE, stringsAsFactors = FALSE)
dataset$Sample <- as.character(sample)
dataset$SVTYPE <- mutation_type

filtered_dataset = apply.filter(dataset, filter_info, mutation_type)
create.tumour.bedfiles(filtered_dataset, output_tab, output_bed)
