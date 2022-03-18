make_differential_expression_wilcox <- function(seuratobject,...){
  
  #calls list_wilcox_pairwise to calculate differential expression using wilcox test
  seurat_list <- SplitObject(seuratobject,split.by ="timepoint")
  genelist <- as.list(row.names(seuratobject@assays$RNA@data))
  
  wilcox_info <- pbapply::pbsapply(genelist,list_wilcox_pairwise,seurat_list,cl=8,simplify = F)
  wilcox_info.df <- do.call(rbind,wilcox_info)
  return(wilcox_info.df)
}

list_wilcox_pairwise <- function(gene,seurat_list,...){
  #pull out numeric vector from each timepoint for given gene
  d.0hr <- as.numeric(seurat_list$`0hr`@assays$RNA@data[gene,])
  d.2hr <-as.numeric(seurat_list$`2hr`@assays$RNA@data[gene,])
  d.4hr <-as.numeric(seurat_list$`4hr`@assays$RNA@data[gene,])
  d.8hr <-as.numeric(seurat_list$`8hr`@assays$RNA@data[gene,])
  
  #calculate wilcox test between sequential timepoints
  w.02 <- wilcox.test(d.0hr,d.2hr,...)$p.value
  w.04 <- wilcox.test(d.0hr,d.4hr,...)$p.value
  w.08 <- wilcox.test(d.0hr,d.8hr,...)$p.value
  
  #find difference in mean compared to 0hr determine whether a gene's expression is going up or down
  diff.2vs0 <- mean(d.2hr,na.rm=T)-mean(d.0hr,na.rm=T)
  diff.4vs0 <- mean(d.4hr,na.rm=T)-mean(d.0hr,na.rm=T)
  diff.8vs0 <- mean(d.8hr,na.rm=T)-mean(d.0hr,na.rm=T)
  
  #make dataframe with all info
  wilcox.info <- data.frame(time=c("0hr_vs_2hr","0hr_vs_4hr","0hr_vs_8hr"),
                            wilcox.p=c(w.02,w.04,w.08),
                            diff_late_minus_early=c(diff.2vs0,diff.4vs0,diff.8vs0),
                            gene=gene)
  wilcox.info$p_adj <- p.adjust(wilcox.info$wilcox.p,method="fdr")
  
  
  return(wilcox.info)
}



classifyTrajectories_earlyLate <- function(diffexpress.dataframe){
  print(head(diffexpress.dataframe %>% filter(p_adj <.05)))
  earlyUp <- diffexpress.dataframe %>% filter(p_adj < .05,time=="0hr_vs_2hr",diff_late_minus_early>0)
  gene_trajectory <- earlyUp %>% select(gene) %>% mutate(trajectory="Early Up")
  
  earlyDown <- diffexpress.dataframe %>% filter(p_adj < .05,time=="0hr_vs_2hr",
                                                diff_late_minus_early<0,!gene %in% gene_trajectory$gene) 
  earlydown.df <- data.frame(gene=unique(earlyDown$gene),trajectory="Early Down")
  gene_trajectory <- rbind(gene_trajectory,earlydown.df)
  
  lateUp <- diffexpress.dataframe %>% filter(p_adj < .05,time %in% c("0hr_vs_4hr","0hr_vs_8hr"),
                    diff_late_minus_early>0,!gene %in% gene_trajectory$gene) 
  lateup.df <- data.frame(gene=unique(lateUp$gene),trajectory="Late Up")
  gene_trajectory <- rbind(gene_trajectory,lateup.df)
  
  lateDown <- diffexpress.dataframe %>% filter(p_adj < .05,time %in% c("0hr_vs_4hr","0hr_vs_8hr"),
                    diff_late_minus_early<0,!gene %in% gene_trajectory$gene)
  latedown.df <- data.frame(gene=unique(lateDown$gene),trajectory="Late Down")
  gene_trajectory <- rbind(gene_trajectory,latedown.df)
  gene_trajectory$trajectory <- factor(gene_trajectory$trajectory,levels=c("Early Up","Late Up","Early Down","Late Down"))
  
  return(gene_trajectory)
  
}

