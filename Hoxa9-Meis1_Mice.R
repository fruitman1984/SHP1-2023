setwd('/media/user/sdg/home/qiuguo/10X/SHP1/H9M')

library(Seurat)
library(ggplot2)
library(dplyr)
HM9_ctrl <- Read10X('control')
colnames(HM9_ctrl) <- paste0(colnames(HM9_ctrl),'_WT')
HM9_ctrl <- CreateSeuratObject(HM9_ctrl,names.field = 2,names.delim = '_',min.cells = 10,min.features = 200)
H9M_chemo <- Read10X('chemo')
colnames(H9M_chemo) <- paste0(colnames(H9M_chemo),'_WC')
H9M_chemo <- CreateSeuratObject(H9M_chemo,names.field = 2,names.delim = '_',min.cells = 10,min.features = 200)

sc_H9M <- merge(HM9_ctrl,H9M_chemo)
sc_H9M[["percent.mt"]] <- PercentageFeatureSet(object = sc_H9M, pattern = "^mt-")
p1<-VlnPlot(object = sc_H9M,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.5)
p1
ggsave(p1,filename = "Vln_QC.pdf",width = 20,height = 5)

gene.freq <- do.call("cbind", tapply(sc_H9M@meta.data$nFeature_RNA,sc_H9M@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(sc_H9M@meta.data$nCount_RNA,sc_H9M@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(sc_H9M@meta.data$percent.mt,sc_H9M@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))
rm(gene.freq,rna.freq,mt.freq)
View(freq.combine)

sc_H9M <- subset(sc_H9M, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 1000 & nCount_RNA < 30000  & percent.mt < 20)


sc_H9M.list <- SplitObject(sc_H9M, split.by = "orig.ident")
for (i in 1:length(sc_H9M.list)) {
  sc_H9M.list[[i]] <- NormalizeData(sc_H9M.list[[i]]) %>% FindVariableFeatures()
}
var.features <- SelectIntegrationFeatures(object.list = sc_H9M.list)
sc_H9M <- merge(sc_H9M.list[[1]], y = sc_H9M.list[2 : length(sc_H9M.list)],merge.data = T)
VariableFeatures(sc_H9M) <- var.features
sc_H9M <-  ScaleData(sc_H9M,vars.to.regress = c('nCount_RNA','percent.mt'),feature = rownames(rownames(sc_H9M)))
sc_H9M <- RunPCA(object = sc_H9M,  verbose = F)
ElbowPlot(sc_H9M,ndims = 30)
sc_H9M <- FindNeighbors(sc_H9M, reduction = "pca", dims = 1:15)
sc_H9M <- FindClusters(sc_H9M, resolution = c(0.6))
sc_H9M <- RunUMAP(sc_H9M, reduction = "pca", dims = 1:15, reduction.name = "umap")
sc_H9M@meta.data$orig.ident <- factor(sc_H9M@meta.data$orig.ident,levels = c('WT','WC'))
DimPlot(sc_H9M,split.by = 'orig.ident',label = T,raster=T,pt.size=2)
ggsave('Umap_sc_H9M.pdf',width = 8,height = 4)

markers <- FindAllMarkers(sc_H9M,only.pos = T,logfc.threshold = 0.5)
write.csv(markers,file = 'cluster_marker.csv',row.names = F)


col_cluster=c('#c53d43','#ec6d71','#c97586','#59b9c6','#2ca9e1','#007bbb','#1e50a2',
              "#b0ca71", "#a9c96c","#a2c867","#9bc762","#94c55d","#8dc358","#86c153","#7fbe4e","#78bc49","#71ba44",#"#6ab83f",
              "#f0f0f0", "#e8e8e8", "#e0e0e0","#d8d8d8","#d0d0d0", "#c8c8c8", "#c0c0c0", "#b8b8b8","#b0b0b0","#a8a8a8", # 10
              "#a0a0a0", "#989898", "#909090","#888888","#808080","#787878","#707070","#686868","#606060","#585858", # 20
              "#505050", "#484848") 

colrs2 <- colrs[c(4,1,5,6,7,2,3,8,9,10,11,12,13,14,15,16)]

sc_H9M@meta.data$cluster <- as.numeric(sc_H9M@meta.data$seurat_clusters)-1
sc_H9M@meta.data$cluster <- factor(sc_H9M@meta.data$cluster,levels = c(5,1,0,3,4,6,11,2,7,9,8,10,12,13,14))
colrs <- c("#54A5D7","#FF5858",
           "#d0d0d0", "#c0c0c0", "#b8b8b8","#a0a0a0", "#909090",
           "#C1D8AC","#A8BF93","#A8C97F","#b0ca71", "#a9c96c","#a2c867","#93B881","#769164")
names(colrs) <- c(5,1,0,3,4,6,11,2,7,9,8,10,12,13,14)
DimPlot(sc_H9M,split.by = 'orig.ident',group.by = 'cluster',label = T,raster=T,pt.size=3,cols = colrs)
ggsave('Umap_sc_H9M_cluster_size3.pdf',width = 8,height = 4)

DotPlot(sc_H9M,features = c('Kit','Sox4','Cd34','Cd93','Myb','Myc','Hoxa9','Meis1','Flt3','Ptpn6','Pfkp'),group.by = 'cluster')
ggsave('Dot_LSC_marker.pdf',width = 8,height = 5)

sc_H9M@meta.data$celltype <- ifelse(sc_H9M@meta.data$cluster == 5,'crLSC',
                                    ifelse(sc_H9M@meta.data$cluster == 1,'cycling LSC',
                                           ifelse(sc_H9M@meta.data$cluster %in% c(0,3,4,6,11),'Blast','TME')))


######## H9M GFP --------
sc_H9M <- readRDS('sc_H9M.rds')
GFP1 <- read.table('SC86-Hoxa9-MeisI-GFP.xls',sep = '\t',header = T,row.names = 1)
GFP2 <- read.table('SC94-Hoxa9-MeisI-GFP.xls',sep = '\t',header = T,row.names = 1)
colnames(GFP1) <- gsub('Chemo_','',colnames(GFP1))
colnames(GFP1) <- paste0(colnames(GFP1),'_WT')
colnames(GFP2) <- gsub('SC94_','',colnames(GFP2))
colnames(GFP2) <- paste0(colnames(GFP2),'_WC')

GFP <- cbind(GFP1,GFP2)

GFP <- GFP[, Cells(GFP)[Cells(GFP) %in% Cells(sc_H9M)]]
rownames(GFP) <- 'H9M'
GFP_t <- t(GFP) %>% as.data.frame()


meta <- sc_H9M@meta.data[,c("orig.ident","seurat_clusters")]
GFP_cluster <- meta[rownames(GFP_t),] %>% cbind(GFP_t)

###### Magic ----
library(reticulate)
use_condaenv('Rmagic')
library(Rmagic)
sc_H9M <- magic(sc_H9M)


####### Heatmap
sc_H9M <- readRDS('/media/user/sdg/home/qiuguo/10X/SHP1/H9M/sc_H9M.rds')

sc_H9M_lsc <- subset(sc_H9M, celltype %in% c('crLSC','cycling LSC','Blast'))
sc_H9M_lsc <- ScaleData(sc_H9M_lsc,assay = 'MAGIC_RNA',vars.to.regress = c('nCount_RNA','percent.mt'))

Idents(sc_H9M) <- 'celltype'
future::plan('multisession',workers=8)
celltype_marker <- FindAllMarkers(sc_H9M_lsc,only.pos = F)
write.csv(celltype_marker,file = 'celltype_marker.csv')

s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes


genes <- c('Sox4','Flt3','Cd93','Rbx1','Cd34','Spink2','Egfl7',
           'Ly6g','Itgam','Ltf','Ifitm2','Ifitm6','Lrg1','Mmp8',
           'Mki67','Top2a','Cdk12','Cdk8','Ccnl2','Ccnb1ip1',
           'Cox10','Pdk1','Tcirg1','Afg3l2','Aaas','Tpi1','Alas1',
           'Cd47','Cd52','Bst2','Mif','B2m','Cd274','Tapbp','Emc6','Ifnar2',
           'Ifit3','Fas','H2-T24','Cxcr4','Entpd1','Tnfsf14','Itsn2','Tnf','Ulbp1')

# color = colorRampPalette(c("#000079","#000093","#0000C6","#0000C6","white",
#                            "#AE0000","#930000","#930000","#4D0000"))(1000)

DoHeatmap(sc_H9M_lsc,assay = 'MAGIC_RNA', features = genes,size = 2,group.by = 'cluster',draw.lines = F,group.colors = colrs,raster = T,raster.dpi=c(3000,3000))+scale_fill_gradientn(colors = c('navy','white','firebrick3'))
ggsave('Heatmap1_H9M.pdf',width=6,height = 4,dpi = 'retina')

####### Heatmap
sc_H9M <- readRDS('/media/user/sdg/home/qiuguo/10X/SHP1/H9M/sc_H9M.rds')

sc_H9M <- ScaleData(sc_H9M,assay = 'MAGIC_RNA',vars.to.regress = c('nCount_RNA','percent.mt'))


sc_H9M_lsc <- subset(sc_H9M, cluster %in% c(1,5,6,11))
sc_H9M_lsc <- ScaleData(sc_H9M_lsc,assay = 'MAGIC_RNA',vars.to.regress = c('nCount_RNA'))
sc_H9M_lsc <- ScaleData(sc_H9M_lsc,assay = 'MAGIC_RNA')

Idents(sc_H9M) <- 'cluster'
future::plan('multisession',workers=8)
celltype_marker <- FindAllMarkers(sc_H9M_lsc,only.pos = F)
write.csv(celltype_marker,file = 'lsc2_celltype_marker.csv')

s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes


genes <- c('Sox4','Flt3','Cd93','Rbx1','Cd34','Spink2','Egfl7',
           'Cebpe','Itgam','Ltf','S100a8','S100a9','Cd177','Ngp',
           'Mki67','Top2a','Cdca2','Cdc7','Cdca8','Cdk11b','E2f2',
           'Adpgk','Hk3','Alas1','Pdk1','Dnm1l','Tcirg1','Nup107',#'Afg3l2','Tpi1',
           'Cd47','Cd52','Bst2','Mif','B2m','Cd274','Tapbp','Emc6','Ifnar2',
           'H2-T24','H2-Q4','Tnf','Cxcr2','Igf1r','Itsn2','Tnfsf13b')

# color = colorRampPalette(c("#000079","#000093","#0000C6","#0000C6","white",
#                            "#AE0000","#930000","#930000","#4D0000"))(1000)
colrs2 <- c('5'="#54A5D7",'1'="#FF5858",'6'="#a0a0a0",'11'= "#a0a0a1")
DoHeatmap(sc_H9M_lsc,assay = 'MAGIC_RNA', features = genes,size = 2,group.by = 'cluster',draw.lines = F,group.colors = colrs2,raster = T)+scale_fill_gradientn(colors = c('navy','white','firebrick3'))
ggsave('Heatmap1_H9M2.pdf',width=6,height = 4,dpi = 'retina')

##### GSVA -----
setwd('/media/user/sdg/home/qiuguo/10X/SHP1/H9M/GSEA')
de_gsva <- function(exprSet,meta,compare = NULL){
  allDiff = list()
  design <- model.matrix(~0+factor(meta))
  colnames(design)=levels(factor(meta))
  rownames(design)=colnames(exprSet)
  
  fit <- lmFit(exprSet,design)
  if(length(unique(meta))==2){
    if(is.null(compare)){
      stop("there are 2 Groups,Please set  compare value")
    }
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
    
  }else if(length(unique(meta))>2){
    for(g in colnames(design)){
      fm = ""
      for(gother in colnames(design)[which(!colnames(design) %in% g)]){
        fm = paste0(fm,"+",gother)
      } 
      
      fm = paste0(g,"VsOthers = ",g,"-(",substring(fm,2),")/",ncol(design)-1)
      contrast.matrix <- makeContrasts(contrasts = fm,levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix) 
      fit2 <- eBayes(fit2)
      allDiff[[g]]=topTable(fit2,adjust='fdr',coef=1,number=Inf)
    }
  }else{
    stop("error only have one group")
  }
  
  return(allDiff)
}
expr = as.matrix(sc_H9M_lsc@assays$RNA@data)
C2_df = msigdbr(species = "Mus musculus", category = "C2" ) %>% split(x = .$gene_symbol, f = .$gs_name)
KEGG_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") %>% split(x = .$gene_symbol, f = .$gs_name)
GO_df = msigdbr(species = "Mus musculus", category = "C5",subcategory = 'BP') %>% split(x = .$gene_symbol, f = .$gs_name)
H_df = msigdbr(species = "Mus musculus", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)
msgdH = msigdbr(species = "Mus musculus", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)
msgdReact = msigdbr(species = "Mus musculus", category = "C2",subcategory = "REACTOME") %>% split(x = .$gene_symbol, f = .$gs_name)
Hallmark_LSC <- gsva(expr, gset.idx.list = H_df, kcdf = "Gaussian", method = "zscore", parallel.sz = 12)

sc_H9M_lsc <- AddMetaData(sc_H9M_lsc,metadata = t(Hallmark_LSC))
sc_H9M_lsc$HALLMARK_G2M_CHECKPOINT
sc_H9M_lsc$HALLMARK_OXIDATIVE_PHOSPHORYLATION
VlnPlot(sc_H9M_lsc,features = c('HALLMARK_G2M_CHECKPOINT'),cols = colrs2,pt.size = 0)+NoLegend()
VlnPlot(sc_H9M_lsc,features = c('HALLMARK_OXIDATIVE_PHOSPHORYLATION'),cols = colrs2,pt.size = 0)+NoLegend()
VlnPlot(sc_H9M_lsc,features = c('HALLMARK_G2M_CHECKPOINT'),cols = colrs2,pt.size = 0)+NoLegend()
VlnPlot(sc_H9M_lsc,features = c('HALLMARK_G2M_CHECKPOINT'),cols = colrs2,pt.size = 0)+NoLegend()

######## Prop -----
cell.prop<-as.data.frame(prop.table(table(sc_H9M$cluster, sc_H9M$orig.ident),margin = 2))
colnames(cell.prop)<-c("cluster","sample","proportion")

ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill",width = 0.5)+
  ggtitle("")+
  theme_bw()+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), 
        axis.ticks.length=unit(0.25,'cm'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title=NULL))+ scale_fill_manual(values=colrs)
ggsave('Prop_H9M_bar.pdf',width = 3,height = 4)


DefaultAssay(sc_H9M_lsc) <- 'MAGIC_RNA'
data<- FetchData(sc_H9M_lsc,vars = c("celltype","Ptpn6"))
data <- subset(data,celltype %in% c('crLSC','cycling LSC'))
cols <- c('crLSC'="#54A5D7",'cycling LSC'='#FF5858','Blast'='#b8b8b8')
ggplot(data, aes(x=celltype,y=Ptpn6,colour = celltype)) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=10,hjust = 1,vjust=0.5))+
  geom_jitter(pch=19,cex=1, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0),aes(color = factor(celltype)))+scale_color_manual(values = cols) +
  NoLegend()+theme(plot.title = element_text(hjust = 0.5))

