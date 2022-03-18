library(tools)
library(cowplot)
library(Seurat)
library(boot)
library(RColorBrewer)
library(Hmisc)
library(venn)
library(enrichR)
library(readr)
source('~/R_startup_source.R')
path.ishi = "/Users/matt/Documents/Research2/singleCell_e2TimeCourse/outs/filtered_feature_bc_matrix/"

#####f
###import and normalize data from output of CellRanger using Seurat
tc_ishi <- Read10X(data.dir = path.ishi)
tc_ishi <- CreateSeuratObject(tc_ishi)

##default is to o	Feature counts for each cell are divided by the total counts for 
#that cell and multiplied by the scale.factor=10000. This is then natural-log transformed using log1p
tc_ishi <- NormalizeData(tc_ishi)

tc_ishi <- PercentageFeatureSet(tc_ishi, pattern = "^MT-", col.name = "percent.mt")
tc_ishi <- SetIdent(tc_ishi,value = "orig.ident")

tc_ishi <- subset(tc_ishi, percent.mt < 7 & nFeature_RNA>2500 & nFeature_RNA < 6000 & nCount_RNA>10000 & nCount_RNA<35000)
tc_ishi <- ScaleData(tc_ishi)

# filter Seurat Object for genes with mean expresion levels > .01
tc_ishi.express <- rowMeans(tc_ishi@assays$RNA@data)
tc_ishihighexpress <- names(tc_ishi.express[tc_ishi.express>.01])
tc_ishi <- subset(tc_ishi,features=tc_ishihighexpress)


# join Seurat metadata with multiseq data to determine barcodes for each cell
load("~/Documents/Research2/singleCell_e2TimeCourse/MULTIseq/final_multiseq_classification.Robj",verbose = T)

tc_ishi@meta.data$cells <- gsub("-1","",row.names(tc_ishi@meta.data))

ishi_meta <- left_join(tc_ishi@meta.data,final_multiseq_classification,by="cells")
row.names(ishi_meta) <- row.names(tc_ishi@meta.data)

#change barcode to timepoint
barcode_timepoint_map <- c("Bar1"="8hr","Bar2"="4hr","Bar3"="2hr","Bar4"="0hr")
ishi_meta$timepoint <- barcode_timepoint_map[ishi_meta$calls]
ishi_meta$timepoint <- factor(ishi_meta$timepoint,levels = c("0hr","2hr","4hr","8hr"))

#replace metadata in initial frame
tc_ishi@meta.data <- ishi_meta
tc_ishi <- SetIdent(tc_ishi, value="timepoint")

#filter out cells with no mapped timepoint (or multiple timepoints)
validcells <- row.names(tc_ishi@meta.data[!is.na(tc_ishi@meta.data$timepoint),])
tc_ishi <- subset(tc_ishi,cells = validcells)
                  
#plot example gene
RidgePlot(tc_ishi,"TGFA")

write_rds(tc_ishi,file="~/Documents/Research2/papers/SingleCellPaper_code/tc_ishi.rds")


##### 
#T-47D initialization
path.t47d = "/Users/matt/Documents/Research2/singleCell_e2TimeCourse/T47D/filtered_feature_bc_matrix/"

tc_t47d <- Read10X(data.dir = path.t47d)
tc_t47d <- CreateSeuratObject(tc_t47d)

##default is to o	Feature counts for each cell are divided by the total counts for 
#that cell and multiplied by the scale.factor=10000. This is then natural-log transformed using log1p
tc_t47d <- NormalizeData(tc_t47d)

tc_t47d <- PercentageFeatureSet(tc_t47d, pattern = "^MT-", col.name = "percent.mt")
tc_t47d <- SetIdent(tc_t47d,value = "orig.ident")

tc_t47d <-  subset(tc_t47d, percent.mt < 20 & nFeature_RNA>3000 & nFeature_RNA < 6000 & nCount_RNA>10000 & nCount_RNA<40000)
#subset(tc_t47d,subset=nFeature_RNA > 3500 & nFeature_RNA < 5500 & percent.mt < 20)
tc_t47d <- ScaleData(tc_t47d)

# filter Seurat Object for genes with mean expresion levels > .01
tc_t47d.express <- rowMeans(tc_t47d@assays$RNA@data)
tc_t47dhighexpress <- names(tc_t47d.express[tc_t47d.express>.01])
tc_t47d <- subset(tc_t47d,features=tc_t47dhighexpress)



# join Seurat metadata with multiseq data to determine barcodes for each cell

tc_t47d@meta.data$cells <- gsub("-1","",row.names(tc_t47d@meta.data))
# write.table(tc_t47d@meta.data$cells,"~/Documents/Research2/singleCell_e2TimeCourse/MULTIseq/t47d_cell_barcodes.tsv",quote = F,sep="\t",
#             row.names = F,col.names = F)

load("~/Documents/Research2/singleCell_e2TimeCourse/T47D/MULTIseq/t47d_final_multiseq_classification.Robj",verbose = T)
t47d_meta <- left_join(tc_t47d@meta.data,t47d_final_multiseq_classification,by="cells")
row.names(t47d_meta) <- row.names(tc_t47d@meta.data)

#change barcode to timepoint
barcode_timepoint_map <- c("Bar1"="8hr","Bar2"="4hr","Bar3"="2hr","Bar4"="0hr")
t47d_meta$timepoint <- barcode_timepoint_map[t47d_meta$calls]
t47d_meta$timepoint <- factor(t47d_meta$timepoint,levels = c("0hr","2hr","4hr","8hr"))

#replace metadata in initial frame
tc_t47d@meta.data <- t47d_meta
tc_t47d <- SetIdent(tc_t47d, value="timepoint")

#filter out cells with no mapped timepoint (or multiple timepoints)
validcells <- row.names(tc_t47d@meta.data[!is.na(tc_t47d@meta.data$timepoint),])
tc_t47d <- subset(tc_t47d,cells = validcells)

#plot example gene
RidgePlot(tc_t47d,"PGR")

#####

### Calling gene trajectories - early, late, up and down
source("~/Documents/Research2/papers/SingleCellPaper_code/reading_and_initializing_10xData_functions.R")


### calling Ishikawa Trajectories
ishi_wilcox <- make_differential_expression_wilcox(tc_ishi)
ishi_trajectories <- classifyTrajectories_earlyLate(ishi_wilcox)
dim(ishi_trajectories)
table(ishi_trajectories$trajectory)
ishikawa_mean_trajmap <- ishi_trajectories$trajectory
names(ishikawa_mean_trajmap) <- ishi_trajectories$gene
load("~/Documents/Research2/singleCell_e2TimeCourse/code/ishi_geneStats_2-14-22.Robj",verbose=T)

### calling T-47D Trajectories
t47d_wilcox <- make_differential_expression_wilcox(tc_t47d)
t47d_trajectories <- classifyTrajectories_earlyLate(t47d_wilcox)
dim(t47d_trajectories)
table(t47d_trajectories$trajectory)
t47d_mean_trajmap <- t47d_trajectories$trajectory
names(t47d_mean_trajmap) <- ishi_trajectories$gene
load("~/Documents/Research2/singleCell_e2TimeCourse/code/t47d_geneStats_2-14-22.Robj",verbose=T)

t47d_geneStats <- left_join(t47d_geneStats %>% select(-trajectory),t47d_trajectories,by="gene") %>% group_by(gene) %>% 
  mutate(mean.scale=scale(mean))
ggplot(t47d_geneStats,aes(x=time,y=mean.scale)) + geom_boxplot()+facet_wrap(~trajectory)
VariableFeatures(tc_ishi) <- reg_genes_ishi
tc_ishi <- RunPCA(tc_ishi,features=reg_genes_ishi)
tc_ishi <- RunUMAP(tc_ishi, dims = c(1:37))

pdf("~/Documents/Research2/singleCell_e2TimeCourse/plots/Ishikawa_feature_count_percentMT.vlns.pdf",width=3,height=4)
VlnPlot(tc_ishi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + theme(text=element_text(size=8)) + xlab("Ishikawa")
VlnPlot(tc_t47d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + theme(text=element_text(size=8)) + xlab("T-47D")

dev.off()
