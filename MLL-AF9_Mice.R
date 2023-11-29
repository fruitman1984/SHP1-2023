##### HOXA9-MEIS1 -----
mat_dir_WT<-"/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/final_cellranger_outs/WT/WT/outs/filtered_feature_bc_matrix"
mat_dir_WC<-"/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/final_cellranger_outs/WC/WC/outs/filtered_feature_bc_matrix"
mat_dir_KO<-"/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/final_cellranger_outs/KO/KO/outs/filtered_feature_bc_matrix"
mat_dir_KC<-"/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/final_cellranger_outs/KC/KC/outs/filtered_feature_bc_matrix"

setwd("/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/analysis/output")

# QC and dim decision
WT_scoreJC <- sc_10x(WT_10x,"WT")
WC_scoreJC <- sc_10x(mat_dir_WC,"WC")
KO_scoreJC <- sc_10x(mat_dir_KO,"KO")
KC_scoreJC <- sc_10x(mat_dir_KC,"KC")

WT_umap <- redu_10x(sc_ScoreJS = WT_scoreJC,dim = 15,proj_name = "WT")
# read 10x data
WT_10x <- Read10X(mat_dir_WT)
WC_10x <- Read10X(mat_dir_WC)
KO_10x <- Read10X(mat_dir_KO)
KC_10x <- Read10X(mat_dir_KC)

# create expression matrix
WT_mat <- data.frame(WT_10x)
WC_mat <- data.frame(WC_10x)
KO_mat <- data.frame(KO_10x)
KC_mat <- data.frame(KC_10x)

# find cell numbers
length(names(WT_mat))
length(names(WC_mat))
length(names(KO_mat))
length(names(KC_mat))

# rename cell names
names(WT_mat)<-paste0("WT_",formatC(1:length(names(WT_mat)),flag = "0",width = 5))
names(WC_mat)<-paste0("WC_",formatC(1:length(names(WC_mat)),flag = "0",width = 5))
names(KO_mat)<-paste0("KO_",formatC(1:length(names(KO_mat)),flag = "0",width = 5))
names(KC_mat)<-paste0("KC_",formatC(1:length(names(KC_mat)),flag = "0",width = 5))

# merge four matrix
lsc_df<-cbind(WT_mat,WC_mat,KO_mat,KC_mat)
#write.table(lsc_df,"/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/analysis/WT_WC_KO_KC_10x_matrix.txt",sep = "\t",quote = F,row.names = T)

# Creat Seurat object
lsc_sc<-CreateSeuratObject(counts = lsc_df,project = "LSC",min.cells = 3,min.features = 200);lsc_sc
lsc_sc[["percent.mt"]] <- PercentageFeatureSet(object = lsc_sc, pattern = "^mt-")
p1<-VlnPlot(object = lsc_sc,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.5)
p1
ggsave(p1,filename = paste0(proj_name,"_01_","FeatureScatter.pdf"),width = 20,height = 5)
sc_fil <- subset(lsc_sc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5);sc_fil
sc_nor<-NormalizeData(sc_fil, normalization.method = "LogNormalize", scale.factor = 15000)
sc_fea <- FindVariableFeatures(sc_nor, selection.method = "vst", nfeatures = 2000)
sc_scale<-ScaleData(sc_fea)
sc_scale <- ScaleData(sc_scale, vars.to.regress = "percent.mt")
sc_pca<-RunPCA(sc_scale, features = VariableFeatures(object = sc_scale))
sc_JS <- JackStraw(sc_pca, num.replicate = 100)
sc_ScoreJS <- ScoreJackStraw(sc_JS, dims = 1:20)
sc_cluster <- FindNeighbors(sc_ScoreJS, dims = 1:15)
sc_cluster <- FindClusters(sc_cluster, resolution = 0.5)

# dim reduction, umap and tsne
sc_umap <- RunUMAP(sc_cluster, dims = 1:15)
sc_umap$cellgroup <- ifelse(sc_umap$seurat_clusters %in% c(1,5,6),'crLSC',
                            ifelse(sc_umap$seurat_clusters %in% c(0,2,3),'cyclingLSC',
                                   ifelse(sc_umap$seurat_clusters %in% c(4,8,11,7),'Blast','Immune')))
colrs <- c("#FF3F3F","#2D82B5","#FF5858","#FF7272","#9B9B9B","#54A5D7","#8ACBF7","#B6B6B6","#D1D1D1","#D6E9CA","#C1D8AC","#CCCCCC","#A8C97F","#A8BF93","#93B881","#769164")
# colrs2 <- colrs[c(4,1,5,6,7,2,3,8,9,10,11,12,13,14,15,16)] 
clst_umap<-DimPlot(sc_umap, reduction = "umap",pt.size = 1,label = T,
                   cols = colrs)
ggsave(clst_umap,filename = paste0("LSC","_09_UMAP_dim","15",".pdf"),width = 6, height = 4)

split.plot<-DimPlot(sc_umap, reduction = "umap",pt.size = 0.005,label = T,
                    #cells.highlight = names(WT_mat),
                    sizes.highlight = F,
                    split.by = "orig.ident",
                    cols = colrs2)
ggsave(split.plot,filename = "split.dim.plot.pdf",width = 20, height = 4.5)


###### MA9 Expression
library(Rmisc)
library(Seurat)
library(ggplot2)
library(dplyr)
sc.exp<- FetchData(sc_umap,vars = c("seurat_clusters","MA9"),slot = 'counts')
data_summary <- summarySE(sc.exp, measurevar="MA9", groupvars=c("seurat_clusters"))
ggplot(sc.exp,aes(x=seurat_clusters,y=MA9,fill=seurat_clusters))+
  geom_violin(trim = T,color="white",scale = "width")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",size=10),
        axis.text.y=element_text(size=16,face="plain"),
        axis.title.y=element_text(size = 20,face="plain"), 
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
        legend.text=element_text(colour="black", 
                                 size=16),
        legend.title=element_text( colour="black",
                                   size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylab("Expression")+xlab("")+
  geom_point(data = data_summary,aes(x=seurat_clusters, y=MA9),pch=19,position=position_dodge(0.3),size=3,col = "grey30")+
  geom_jitter(size = 0.05)+
  scale_fill_manual(values = colrs2)

###### proportion bar -----
cell.prop<-as.data.frame(prop.table(table(sc_umap$seurat_clusters, sc_umap$orig.ident),margin = 2))
colnames(cell.prop)<-c("cluster","sample","proportion")

ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title=NULL))+ scale_fill_manual(values=colrs2)


#### pheatmap 

# threemks <- read.table("/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/analysis/output/rLSC.pLSC.gene.txt")
# lsc_genes <- read.table('/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/genesets/lsc_gene.txt')[[1]]
# metabolic <- read.table('/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/genesets/metabolic.txt')[[1]]
# immune <- read.table('/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/genesets/immune.txt')[[1]]
# threemks <- c(lsc_genes,metabolic,immune)
threemks <- read.table("genes.txt")[[1]]
meta <- sc_umap@meta.data
wt.cell <- rownames(meta)[grep("WT",rownames(meta))]
wt.umap<-subset(sc_umap,cells = wt.cell)
sc.cr <- subset(wt.umap,idents = c("1","5","6"),features = threemks)
sc.nl <- subset(wt.umap,idents = c("0","2","3"),features = threemks)
sc.bc <- subset(wt.umap,idents = c("4","7","8","11"),features = threemks)
gap <- c(898,898+4889)
ht <- cbind(as.data.frame(sc.cr@assays$RNA@data)[,sample(1:898,898,replace = F)],as.data.frame(sc.nl@assays$RNA@data)[,sample(1:4889,4889,replace = F)],as.data.frame(sc.bc@assays$RNA@data)[,sample(1:2038,2038,replace = F)])
ht <- ht[threemks,]
color = colorRampPalette(c("#000079","#000093","#0000C6","#0000C6","white",
                           "#AE0000","#930000","#930000","#4D0000"))(1000)
p <- pheatmap(ht,
              cluster_cols = F,
              cluster_rows = F,
              scale = "row",
              show_colnames = F,
              # cellheight = 8,
              # cellwidth = 0.1,
              gaps_col = gap,
              color = color)
p
ggsave(p,filename = "/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/analysis/output/heatmap.3cluster.deg.pdf",width = 15, height = 15)


# pearson corelation -----
## gene matrix
pearson.exp = data.table(data.frame(sc_umap@assays$RNA@data),keep.rownames = T,key = "rn")
## gene sets
lsc.gene <- read.table("/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/analysis/output/pearson.geneset.txt")
lsc.gene <- as.vector(lsc.gene$V1)
lsc.gene <- lsc.gene[!duplicated(lsc.gene)]
lsc.gene
t.gene <- pearson.exp$rn

ht <- pearson.heatmap(pearson.mat = pearson.exp, lsc.gene = t.gene)

pearson.heatmap<-function(pearson.mat,lsc.gene, cluster.rows = T, cluster.cols = T){
  pearson.mat = pearson.mat[lsc.gene]
  pearson.ava = data.frame(rep(NA,length(lsc.gene)))
  rownames(pearson.ava) = lsc.gene
  
  for (i in clst.ls){
    q = data.frame(apply(pearson.mat[,..i],1,mean))
    pearson.ava = cbind(pearson.ava,q)
  }
  
  pearson.ava = pearson.ava[,-1]
  names(pearson.ava) = paste0("Clst.",c(0:15))
  pearson.ava
  pco = c()
  for (k in 1:16){
    for (j in 1:16){
      co = cor(pearson.ava[,k],pearson.ava[,j],method = "pearson")
      pco = c(pco,co)
    }
  }
  pco = data.frame(matrix(pco,nrow = 16,byrow = T))
  names(pco) = paste0("Clst.",c(0:15))
  rownames(pco) = paste0("Clst.",c(0:15))
  pco
  pco2 = pco[c(2,6,7,1,4,3,5,8,9,10,11,12,13,14,15,16),c(2,6,7,1,4,3,5,8,9,10,11,12,13,14,15,16)]
  pl = pheatmap(pco2,
                cluster_rows = cluster.rows,
                cluster_cols = cluster.cols)
  return(pl)
}


###### GSEA -----------
.libPaths("/home/wj/R/x86_64-pc-linux-gnu-library/reg")
library("utils")
library("tools")
#library("dplyr")
#library("GSEA")
# load GSEA package
pkgs<-dir("/media/user/sda/R/GSEA/GSEA_R-master/R")
for (i in pkgs){
  source(paste0('/media/user/sda/R/GSEA/GSEA_R-master/R/',i))
  print(paste("LOADING:",i))
}

cat("\n")
cat("Starting...\n")
cat("\n")
###
paths <- "/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/analysis/output/GSEA/input/"
# Input path to GCT or RNK formatted gene expression dataset (or drop file into R window):
inputds <- paste0(paths,"rLSC.vs.pLSC.wt.wc.gct")
# Input path to experiment CLS file (or drop file into R window):
inputcls <- paste0(paths,"rLSC.vs.pLSC.wt.wc.cls")
# Input path to GMT formatted gene set database (or drop file into R window):
gsdb <- paste0(paths,"genesets.gmt")

# Drop a directory into R window to use as the output folder or enter directory path:
outdir <- "/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/analysis/output/GSEA/out" 
# Enter a prefix to label output files:
outname <- "rplsc.wt.kc"

###
#Collapse data set to Gene Symbols?
collapsedataset <- FALSE
inputchip <- "NOCHIP"
collapsemode <- "NOCOLLAPSE"
# Select GSEA Permutation Type (recommended: Phenotype)": c("Phenotype", "gene_set")
reshuffetype <- "Phenotype"; permutation <- "sample.labels"
#reshuffetype <- "gene_set"; permutation <- "gene.labels"

# Use default number of permutations for significance testing? (default: 1000)
nperms <- 10
#  Use default maximum gene set size filter? (default: 500 genes)
maxsize <- 1000
# Use default minimum gene set size filter? (default: 15 genes)
minsize <- 15
# Use default Signal2Noise metric for ranking genes? "S2N" or "ttest"
rankmetric <- "S2N"

rankmethod <- "GSEA"
#rankmethod <- "preranked"

GSEA(
  # Input/Output Files 
  input.ds = inputds,                    # Input gene expression dataset file in GCT format
  input.cls = inputcls,                  # Input class vector (phenotype) file in CLS format
  gs.db = gsdb,                          # Gene set database in GMT format
  input.chip = inputchip,               # CHIP File
  output.directory      = outdir,        # Directory where to store output and results (default: "")
  #  Program parameters 
  doc.string            = outname,         # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
  reshuffling.type      = permutation,     # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
  nperm                 = as.integer(nperms),            # Number of random permutations (default: 1000)
  weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
  nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
  fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
  fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
  topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
  adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
  gs.size.threshold.min = as.integer(minsize),         # Minimum size (in genes) for database gene sets to be considered (default: 25)
  gs.size.threshold.max = as.integer(maxsize),         # Maximum size (in genes) for database gene sets to be considered (default: 500)
  reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
  preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
  random.seed           = as.integer(as.POSIXct(Sys.time())),            # Random number generator seed. (default: 123456)
  perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
  fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
  replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
  collapse.dataset      = collapsedataset, # Collapse dataset to gene symbols using a user provided chip file (default: F)
  collapse.mode         = collapsemode,
  save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
  use.fast.enrichment.routine = T,         # Use faster routine to compute enrichment for random permutations (default: T)
  gsea.type = rankmethod,                     # Select Standard GSEA (default) or preranked
  rank.metric = rankmetric
)
# Overlap and leading gene subset assignment analysis of the GSEA results

GSEA.Analyze.Sets(
  directory           = outdir,        # Directory where to store output and results (default: "")
  topgs = 20,                                                           # number of top scoring gene sets used for analysis
  height = 16,
  width = 16,
  gsea.type = rankmethod,
  doc.string = outname
)


######### monocle -------
.libPaths("/home/wj/R/x86_64-pc-linux-gnu-library/sc")
library(monocle)
sample_ann <-  sc_umap@meta.data 
gene_ann <- data.frame(gene_short_name = rownames(sc_umap@assays$RNA),
                       row.names =  rownames(sc_umap@assays$RNA))
#head(gene_ann)
pd <- new("AnnotatedDataFrame",data=sample_ann)
fd <- new("AnnotatedDataFrame",data=gene_ann)
ct=as.data.frame(sc_umap@assays$RNA@counts)

sc_cds <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds <- detectGenes(sc_cds, min_expr = 1) 
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F)

cds <- reduceDimension(cds, max_components = 2, num_dim = 20,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 5) 
plot_cell_clusters(cds, 1, 2 )
table(pData(cds)$cellgroup) 
colnames(pData(cds))

table(pData(cds)$cellgroup)
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~cellgroup")
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 

diff_celltype<- diff_celltype[order(diff_celltype$qval)]
ordering_genes <- row.names(diff_celltype[1:1000],)
#ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
#ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds_redu <- orderCells(cds)
saveRDS(cds_redu,file = "/media/user/sda/data/Xuxi_10xGenomics/2.cellranger_outs/analysis/output/xuxi_cds_redu.rds")


genes <- c('Kit','Cd93','Ikzf2','Flt3','Dpp4',
           'Cd47','Cd24a','Bcl2','Ccnd2','Cdkn3','Gpl1','Pfkp','Lyz2','Cxcr2','Ly6g','Itgam')
genes <- genes[genes %in% rownames(cds_redu )]
plot_pseudotime_heatmap(cds_redu [genes,],show_rownames = T,cluster_rows = F)

pt.matrix <- exprs(cds_redu)[match(genes,rownames(cds_redu)),order(cds_redu$Pseudotime)] %>% as.matrix()
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
# ccc <- colnames(pt.matrix)
# pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
# rownames(pt.matrix) <- genes
# colnames(pt.matrix) <- ccc
library(ComplexHeatmap)
cellInfo <- as.data.frame(subset(sc_umap@meta.data[colnames(pt.matrix),],select = 'cellgroup'))
annotation <- HeatmapAnnotation(df = cellInfo,col = list(cellgroup = c('crLSC' = "Blue", 'cyclingLSC' = "red", 'Blast' = "grey")))
pt.matrix <- pt.matrix[genes,]

Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = c('grey28','white','firebrick3'),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 10),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = F,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  use_raster = T,top_annotation  = annotation)

####### Go term --------
######DGE
fc <- FindMarkers(sc_umap,ident.1 = c('1','5','6'),ident.2 = c("0","3","4"),logfc.threshold = 0)
write.csv(fc,file = 'res_vs_sens_deg.csv')
######GO term
logFC_t=0.15
P.Value_t = 0.05
k1 = (fc$p_val_adj < P.Value_t)&(fc$avg_log2FC < -logFC_t)
k2 = (fc$p_val_adj < P.Value_t)&(fc$avg_log2FC > logFC_t)
table(k1)
table(k2)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
fc$change <- change
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)
#library(biomaRt)
fc$gene <- rownames(fc)
s2e <- bitr(fc$gene, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Mm.eg.db)
fc <- inner_join(fc,s2e,by=c("gene"="SYMBOL"))

#GO database analysis 
gene_up = fc[fc$change == 'up','ENTREZID'] 
gene_down = fc[fc$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = fc[,'ENTREZID']

H_df = msigdbr(species = "Mus musculus", category = "H") %>%  as.data.frame()
term2gene = data.frame(Term = H_df$gs_id, Gene = H_df$entrez_gene)
term2name = data.frame(Term = H_df$gs_id, Name = H_df$gs_name)

ENRICH = enricher(gene = gene_up,pvalueCutoff = 0.9,pAdjustMethod = 'fdr',qvalueCutoff = 0.9,TERM2GENE = term2gene,TERM2NAME = term2name)
write.csv(ENRICH@result,file = 'Hallmark_enrichment.csv')

ego_BP_up <- enrichGO(gene = gene_up,
                      OrgDb= org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
ego_BP_down <- enrichGO(gene = gene_down,
                        OrgDb= org.Mm.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        minGSSize = 1,
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.01,
                        readable = TRUE)
write.csv(ego_BP_up@result,file = 'ego_BP_up_s3g.csv',row.names = F)
write.csv(ego_BP_down@result,file = 'ego_BP_down_s3g.csv',row.names = F)

library(ggstance)
library(forcats)
ego_BP_selected <- read.csv('./ego_BP_selected_s3g.csv')
ego_BP_selected <- read.csv('./ego_up_s3g.csv')
ego_BP_selected$LogP <- abs(log10(ego_BP_selected$pvalue))
ego_BP_selected$Direction <- 'Up'
sortdf <- ego_BP_selected[order(ego_BP_selected$LogP),]
sortdf$Description <- factor(sortdf$Description, levels = sortdf$Description)
ggplot(sortdf, aes(LogP, fct_reorder(Description, LogP), fill=Direction), showCategory=26) + 
  geom_barh(stat='identity',width = 0.8) + 
  # scale_fill_continuous(low='blue', high='red') + 
  scale_fill_manual(values = '#55A0FB') +
  theme_minimal() + ylab("")+
  theme(axis.title = element_text(size = 13), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size=10,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 13.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        panel.grid = element_blank(),axis.line = element_line(),) + 
  scale_x_continuous(expand = c(0, 0), 
                     position = "top") + 
  scale_y_discrete(expand = expansion(add = c(0, 0.5)))

####### PTPN6 Expression -------
library(Rmagic)
reticulate::use_condaenv('Rmagic')
sc_umap <- magic(sc_umap)

DefaultAssay(sc_umap) <- 'MAGIC_RNA'
data<- FetchData(sc_umap,vars = c("LSC","Ptpn6"))
data <- subset(data,LSC %in% c('crLSC','cyclingLSC'))
cols <- c('crLSC'="#54A5D7",'cyclingLSC'='#FF5858','Blast'='#b8b8b8')
ggplot(data, aes(x=LSC,y=Ptpn6,colour = LSC)) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=10,hjust = 1,vjust=0.5))+
  geom_jitter(pch=19,cex=0.5, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0),aes(color = factor(LSC)))+scale_color_manual(values = cols) +
  NoLegend()+theme(plot.title = element_text(hjust = 0.5))
