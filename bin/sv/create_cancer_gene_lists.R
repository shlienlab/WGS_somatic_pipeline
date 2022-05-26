# Version              Date            Developer                        Comments
#-------------------------------------------------------------------------------
#     0.0                 ?                    ?             initial development
#     0.1        2018-12-17        Drew Thompson                 updated for WDL
#   0.1.1        2019-07-31        Drew Thompson        provide output filenames

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("Must provide the final merged tab file, the path to a cosmic cancer gene census R data file, and the output tab file name, and the output CF tab file name.")
}

data_file = args[1]
cosmic_path = args[2]
output_file = args[3]
output_cf_file = args[4]
load(cosmic_path)

pype.to.chainfinder <- function(data){
  data$num <- row.names(data)
  for (i in 1:nrow(data)){
    semicolon.strsplit <- strsplit(as.character(data$Breakpoint1[i]), fixed=TRUE, split=";")
    row.numbahs <- as.data.frame(semicolon.strsplit)
    if(nrow(row.numbahs) > 1){
      inter.dfrm <- t(as.data.frame(semicolon.strsplit))
      open.b.strsplit_1 <- strsplit(inter.dfrm[1], fixed=TRUE, split="(")
      open.b.1.dfrm <- as.data.frame(open.b.strsplit_1)
      
      if(nrow(open.b.1.dfrm) > 1){
        open.b.strsplit_2 <- strsplit(inter.dfrm[2], fixed=TRUE, split="(")
        close.b.1.strsplit <- strsplit(as.character(open.b.1.dfrm[2,]), fixed=TRUE, split=")")
        open.b.2.dfrm <- as.data.frame(open.b.strsplit_2)
        close.b.2.strsplit <- strsplit(as.character(open.b.2.dfrm[2,]), fixed=TRUE, split=")")
        close.b.1.split.dfrm <- as.data.frame(close.b.1.strsplit, fixed=TRUE, split=",")
        comma.b.1.split <- strsplit(as.character(close.b.1.split.dfrm[1,]), fixed=TRUE, split=",")
        close.b.2.split.dfrm <- as.data.frame(close.b.2.strsplit, fixed=TRUE, split=",")
        comma.b.2.split <- strsplit(as.character(close.b.2.split.dfrm[1,]), fixed=TRUE, split=",")
        comma.dfrm1 <- as.data.frame(comma.b.1.split)
        distance1 <- as.numeric(as.character(comma.dfrm1[2,]))
        strandedness1 <- strsplit(as.character(comma.dfrm1[1,]), fixed=TRUE, split="")
        stranded.dfrm1 <- as.data.frame(strandedness1)
        strand1 <- as.character(stranded.dfrm1[1,])
        closest.gene.dfrm1 <- as.character(stranded.dfrm1[2:nrow(stranded.dfrm1),])
        closest.gene1 <- paste(closest.gene.dfrm1, collapse="")
        comma.dfrm2 <- as.data.frame(comma.b.2.split)
        strandedness2 <- strsplit(as.character(comma.dfrm2[1,]), fixed=TRUE, split="")
        distance2 <- as.numeric(as.character(comma.dfrm2[2,]))
        stranded.dfrm2 <- as.data.frame(strandedness2)
        strand2 <- as.character(stranded.dfrm2[1,])
        closest.gene.dfrm2 <- as.character(stranded.dfrm2[2:nrow(stranded.dfrm2),])
        closest.gene2 <- paste(closest.gene.dfrm2, collapse="")
        
        if(distance1 > distance2){
          data$GeneA[i] <- closest.gene2
          data$str1[i] <- strand2
          data$func1[i] <- "intergenic"
        } else {
          data$GeneA[i] <- closest.gene1
          data$str1[i] <- strand1
          data$func1[i] <- "intergenic"
        }
      } else {
        closest.gene_nb1 <- strsplit(as.character(inter.dfrm[1,1]), fixed=TRUE, split="")
        closest.gene_nb2 <- strsplit(as.character(inter.dfrm[1,2]), fixed=TRUE, split="")   
        cg_nb1 <- t(as.data.frame(closest.gene_nb1))
        nb1_strand <- as.character(cg_nb1[1,1])
        nb1_gene <- as.character(paste(sep="", cg_nb1[1,2:ncol(cg_nb1)], collapse=""))
        cg_nb2 <- t(as.data.frame(closest.gene_nb2))
        nb2_strand <- as.character(cg_nb1[1,1])
        nb2_gene <- as.character(paste(sep="", cg_nb2[1,2:ncol(cg_nb2)], collapse=""))
        data$GeneA[i] <- nb1_gene
        data$str1[i] <- nb1_strand
        data$func1[i] <- "intergenic"
      }
    } else {
      gene.split <- strsplit(as.character(semicolon.strsplit), fixed=TRUE, split="")
      gene.dfrm <- t(as.data.frame(gene.split))
      strand.gen <- as.character(gene.dfrm[1,1])
      gene.gene <- paste(gene.dfrm[1,2:ncol(gene.dfrm)], collapse="")
      data$GeneA[i] <- gene.gene
      data$str1[i] <- strand.gen
      data$func1[i] <- "genic"
    }
  }
  
  for (i in 1:nrow(data)){
    semicolon.strsplit <- strsplit(as.character(data$Breakpoint2[i]), fixed=TRUE, split=";")
    row.numbahs <- as.data.frame(semicolon.strsplit)
    if(nrow(row.numbahs) > 1){
      inter.dfrm <- t(as.data.frame(semicolon.strsplit))
      open.b.strsplit_1 <- strsplit(inter.dfrm[1], fixed=TRUE, split="(")
      open.b.1.dfrm <- as.data.frame(open.b.strsplit_1)
      
      if(nrow(open.b.1.dfrm) > 1){
        open.b.strsplit_2 <- strsplit(inter.dfrm[2], fixed=TRUE, split="(")
        close.b.1.strsplit <- strsplit(as.character(open.b.1.dfrm[2,]), fixed=TRUE, split=")")
        open.b.2.dfrm <- as.data.frame(open.b.strsplit_2)
        close.b.2.strsplit <- strsplit(as.character(open.b.2.dfrm[2,]), fixed=TRUE, split=")")
        close.b.1.split.dfrm <- as.data.frame(close.b.1.strsplit, fixed=TRUE, split=",")
        comma.b.1.split <- strsplit(as.character(close.b.1.split.dfrm[1,]), fixed=TRUE, split=",")
        close.b.2.split.dfrm <- as.data.frame(close.b.2.strsplit, fixed=TRUE, split=",")
        comma.b.2.split <- strsplit(as.character(close.b.2.split.dfrm[1,]), fixed=TRUE, split=",")
        comma.dfrm1 <- as.data.frame(comma.b.1.split)
        distance1 <- as.numeric(as.character(comma.dfrm1[2,]))
        strandedness1 <- strsplit(as.character(comma.dfrm1[1,]), fixed=TRUE, split="")
        stranded.dfrm1 <- as.data.frame(strandedness1)
        strand1 <- as.character(stranded.dfrm1[1,])
        closest.gene.dfrm1 <- as.character(stranded.dfrm1[2:nrow(stranded.dfrm1),])
        closest.gene1 <- paste(closest.gene.dfrm1, collapse="")
        comma.dfrm2 <- as.data.frame(comma.b.2.split)
        strandedness2 <- strsplit(as.character(comma.dfrm2[1,]), fixed=TRUE, split="")
        distance2 <- as.numeric(as.character(comma.dfrm2[2,]))
        stranded.dfrm2 <- as.data.frame(strandedness2)
        strand2 <- as.character(stranded.dfrm2[1,])
        closest.gene.dfrm2 <- as.character(stranded.dfrm2[2:nrow(stranded.dfrm2),])
        closest.gene2 <- paste(closest.gene.dfrm2, collapse="")
        
        if(distance1 > distance2){
          data$GeneB[i] <- closest.gene2
          data$str2[i] <- strand2
          data$func2[i] <- "intergenic"
        } else {
          data$GeneB[i] <- closest.gene1
          data$str2[i] <- strand1
          data$func2[i] <- "intergenic"
        }
      } else {
        closest.gene_nb1 <- strsplit(as.character(inter.dfrm[1,1]), fixed=TRUE, split="")
        closest.gene_nb2 <- strsplit(as.character(inter.dfrm[1,2]), fixed=TRUE, split="")   
        cg_nb1 <- t(as.data.frame(closest.gene_nb1))
        nb1_strand <- as.character(cg_nb1[1,1])
        nb1_gene <- as.character(paste(sep="", cg_nb1[1,2:ncol(cg_nb1)], collapse=""))
        cg_nb2 <- t(as.data.frame(closest.gene_nb2))
        nb2_strand <- as.character(cg_nb1[1,1])
        nb2_gene <- as.character(paste(sep="", cg_nb2[1,2:ncol(cg_nb2)], collapse=""))
        data$GeneB[i] <- nb1_gene
        data$str2[i] <- nb1_strand
        data$func2[i] <- "intergenic"
      }
    } else {
      gene.split <- strsplit(as.character(semicolon.strsplit), fixed=TRUE, split="")
      gene.dfrm <- t(as.data.frame(gene.split))
      strand.gen <- as.character(gene.dfrm[1,1])
      gene.gene <- paste(gene.dfrm[1,2:ncol(gene.dfrm)], collapse="")
      data$GeneB[i] <- gene.gene
      data$str2[i] <- strand.gen
      data$func2[i] <- "genic"
    }
  }
  
  for (i in 1:nrow(data)){
    if(data$str1[i] == "+"){
      data$str1[i] <- "0"
    } else {
      data$str1[i] <- "1"
    }
  }
  
  for (i in 1:nrow(data)){
    if(data$str2[i] == "+"){
      data$str2[i] <- "0"
    } else {
      data$str2[i] <- "1"
    }
  }
  
  chainfinder.shlien <- data[,c("Sample", "num", "Chromosome1", "str1", "Position1", "Chromosome2", "str2", "Position2", "GeneA", "GeneB", "func1", "func2")]
  colnames(chainfinder.shlien) <- c("sample", "num", "chr1", "str1", "pos1", "chr2", "str2", "pos2", "site1", "site2", "site1_region", "site2_region")
  return(chainfinder.shlien)
}


temp <- read.table(file=data_file, header=T, sep="\t", stringsAsFactors = F)
if(nrow(temp) > 0){
  chainfinder_format <- pype.to.chainfinder(temp)
  write.table(chainfinder_format, file=output_cf_file, sep="\t", row.names=F, quote=F)
    
  cancer_genes <- chainfinder_format[chainfinder_format$site1 %in% cosmic.cancer.gene.census$gene | chainfinder_format$site2 %in% cosmic.cancer.gene.census$gene,]
  write.table(cancer_genes, file=output_file, sep="\t", row.names=F, quote=F)
} else {
  print("No rearrangements")
}
