library(deMULTIplex)
library(ShortRead)
library(tsne)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(umap)
source("~/R_startup_source.R")


align_multiseq <- function(barcode.path,cells.path,R1.path,R2.path){
  bar.ref <- read.table(barcode.path,stringsAsFactors = F)$V2
  cell.id.vec <- read.table(cells.path,stringsAsFactors = F)$V1
  # 
  # ## Pre-process MULTI-seq sample barcode FASTQs
  readTable <- MULTIseq.preProcess(
    R1 = R1.path,
    R2 = R2.path,
    cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8)
  )
  bar.table <- MULTIseq.align(readTable, cell.id.vec, bar.ref)
  
}

quantile_sweep <- function(bar.table){
  
  bar.table_sweep.list <- list()
  n <- 0
  for (q in seq(0.01, 0.99, by=0.02)) {
    # print(bar.table[,1])
    n <- n + 1
    bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
    names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
  }
  
  threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
  
  p <- ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
    geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
  print(p)
  ## Finalize round 1 classifications, remove negative cells
  round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
  neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
  bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
  return(list(bar.table,round1.calls,neg.cells))
}

dir.create("MULTIseq/outs/")


bar.table.ishi <- align_multiseq(barcode.path = "MULTIseq/multiseq_barcodes.txt",
                                 cells.path = "MULTIseq/ishi_cell_barcodes.tsv",
                                 R1.path = 'MULTIseq/17473X1_191121_A00421_0130_AHH7HGDRXX_S1_L001_R1_001.fastq.gz',
                                 R2.path = 'MULTIseq/17473X1_191121_A00421_0130_AHH7HGDRXX_S1_L001_R2_001.fastq.gz')
# Perform MULTI-seq sample barcode alignment
bar.table.t47d <- align_multiseq(barcode.path = "MULTIseq/multiseq_barcodes.txt",
                                 cells.path = "MULTIseq/t47d_cell_barcodes.tsv",
                                 R1.path='MULTIseq/17899X1_200416_A00421_0186_BHLJ7GDRXX_S1_L002_R1_001.fastq.gz',
                                 R2.path='MULTIseq/17899X1_200416_A00421_0186_BHLJ7GDRXX_S1_L002_R2_001.fastq.gz')

### Run Repeated barcode sweeps until there are no negative calls left
pdf("MULTIseq/outs/Ishi_MULTIseq_barcode_sweep.pdf")
bar.table.ishi <- bar.table.ishi[,1:4]
while(!is.na(table(bar.table.ishi[[2]])["Negative"])){
  print("another round")
 bar.table.ishi <- quantile_sweep(bar.table.ishi)
}
dev.off()
 
ishi.multiseq.calls <- bar.table.ishi[[2]]
save(ishi.multiseq.calls,file="MULTIseq/outs/ishi.multiseq.calls_3-17-22.Robj")

pdf("MULTIseq/outs/T47D_MULTIseq_barcode_sweep.pdf")
bar.table.t47d <- bar.table.t47d[,1:4]
while(!is.na(table(bar.table.t47d[[2]])["Negative"])){
  print("another round")
  bar.table.t47d <- quantile_sweep(bar.table.t47d)
}
dev.off()

t47d.multiseq.calls <- bar.table.t47d.round3[[2]]
save(t47d.multiseq.calls,file="MULTIseq/outs/t47d.multiseq.calls_3-17-22.Robj")







