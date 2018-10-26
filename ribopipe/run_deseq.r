source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

setwd("./") #change this
library(DESeq2)
args = commandArgs(trailingOnly=TRUE)

#import counts_data
count_table <- read.table(args[1],sep=",",header=TRUE,row.names=1)
x_start <- startsWith(colnames(count_table), "X")

#create conditions dataframe
sample_table <- read.table(text=readLines(args[3], warn = FALSE), header=TRUE, sep=",")
if (TRUE %in% x_start) {
  sample_table[["X"]] <- paste("X", sample_table[["X"]], sep="")
  }

#run DESeq2 analysis on data
if (args[6] == 'False') {

  if (args[4] == 'True') {
    if (args[5] == 'riboseq') {
      quit(status=1) # coming soon
    } else if (args[5] == 'rnaseq') {
      quit(status=1) # coming soon
    } else {
      quit(status=1)
    }
  } else {
    if (args[5] == 'riboseq') {
      dds <- DESeqDataSetFromMatrix(countData = count_table, colData = sample_table, design = ~type+condition+type:condition)
    } else if (args[5] == 'rnaseq') {
      dds <- DESeqDataSetFromMatrix(countData = count_table, colData = sample_table, design = ~type)
    } else {
      quit(status=1)
    }
  }
} else {
  dds <- DESeqDataSetFromMatrix(countData = count_table, colData = sample_table, design = args[6])
}




dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]

#write output to new file
write.table(as.data.frame(resOrdered),file=args[2],sep=",",col.names=T,row.names=T)
