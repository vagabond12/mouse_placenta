
pbmc=mouse.troph.combined

pbmc=mouse.combined
.libPaths('C:\\Users\\WHMX\\Documents\\R\\win-library\\3.6')

.libPaths()
rm(list=ls())
rm(pbmc2)
pbmc <- tr13580
pbmc<- pbmc2
pbmc2 <- pbmc
pbmc <- pbmc10
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(object = pbmc, ident =c("ExE ectoderm"),invert=F)])
Idents(object = pbmc) <- "doublet"
Idents(object = pbmc2) <- "C17E"
Idents(object = pbmc) <- "celltype"
Idents(object = pbmc2) <- "stage"
Idents(object = pbmc) <- "C17E"
pbmc <- SubsetData(object = pbmc, ident.remove =  c(6,7))
pbmc <- SubsetData(object = pbmc, ident.use = c("S3"))
pbmc <- SubsetData(object = pbmc, ident.use = c(6,7,"E1","E2","P1","P2","P3"))
pbmc <- SubsetData(object = pbmc, ident.use = c("E1","E2","P1"))
pbmc <- SubsetData(object = pbmc2, ident.use = c(6,7,"P2","P3"))
pbmc <- SubsetData(object = pbmc, ident.use = c("13"))
pbmc <- SubsetData(object = pbmc, ident.remove = c(10:12))
pbmc <- SubsetData(object = pbmc, ident.use =  c(0))

pbmc <- SubsetData(object = pbmc2, ident.use =  c(1:5))
pbmc <- SubsetData(object = pbmc2, ident.use =  c(1,"P1","P2","P3","E2",6:9))
pbmc <- SubsetData(object = pbmc2, ident.use =  c(1,"P1","P2","P3","E1",14:17))

pbmc <- CreateSeuratObject(pbmc$RNA@counts  , meta.data =pbmc@meta.data,min.cells = 3)

pbmc <- SubsetData(object = pbmc2, ident.use =c("E13.5"))
pbmc <- SubsetData(object = pbmc2, ident.use =  c("E07.5","E08.5","E09.5","E10.5"))
pbmc <- SubsetData(object = pbmc2, ident.use =  c("E07.5","E08.5"))
pbmc <- SubsetData(object = pbmc, ident.use =  c("S0", "S1", "S2"))
pbmc10 <- SubsetData(object = pbmc9, ident.use =  c("EPC0", "EPC1", "EPC2"))
pbmcx <-subset(pbmc, idents =1:5, invert = F)

pbmc.data <- as.matrix(x = pbmc$RNA@counts)
meta.data=pbmc@meta.data
pbmc <- CreateSeuratObject(counts = pbmc.data, meta.data =meta.data, min.cells = 3)
pbmc

gcm=pbmc
epc=pbmc
pbmc=epc
pbmc8=pbmc

pbmc <- merge(x = pbmc8, y =c(epc,epc2), project = "mpl-atlas")
pbmc

pbmc10$C17E<- Idents(pbmc)
pbmc$C17E<- Idents(epc)

Idents(pbmc9)<- Idents(pbmc)


pbmc.data=read.csv("mpl12/mtr15682-old.csv",check.names = FALSE, header=TRUE,row.names=1)     
meta.data=read.csv("mpl12/mtr15682-meta-gse.csv",check.names = FALSE, header=TRUE,row.names=1)
pbmc <- CreateSeuratObject(pbmc.data, meta.data =meta.data,min.cells = 3)
pbmc

pbmc@meta.data=meta.data

DimPlot(pbmc,label=F,label.size = 8)+ theme(text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))
x=DimPlot(pbmc,label=T,label.size = 8)+ theme(text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

DimPlot(pbmc,label=T,label.size = 8,group.by = "branch2")+ theme(text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))
DimPlot(pbmc,label=T,label.size = 8,group.by = "C17")+ theme(text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))


library(help=Seurat)
library(help=ggplot2)
library()
version
session_Info()

library(Seurat);library(Matrix);library(dplyr);library(ggplot2)
library(cowplot);library(scater);library(loomR)
library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);suppressMessages(library(SeuratDisk));suppressMessages(library(SeuratWrappers));suppressMessages(library(SeuratData));

pbmc.data <- Read10X(data.dir = "atlas")
meta.data=read.csv("atlas/atlas-meta.csv",check.names = FALSE, header=TRUE,row.names=1)
pbmc <- CreateSeuratObject(counts = pbmc.data, meta.data =meta.data, names.delim = ":", min.cells = 3)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = pbmc))


pancreas=pbmc
pancreas.list <- SplitObject(object = pancreas, split.by = "stage")
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list ,  anchor.features = 4000, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(object = pancreas.integrated)
pancreas.integrated <- RunPCA(object = pancreas.integrated, npcs = 30, verbose = FALSE)
DimPlot(pancreas.integrated)

pbmc= pancreas.integrated

ElbowPlot(object = pancreas.integrated)
pbmc <- RunUMAP(pancreas.integrated, reduction.use = "pca", dims = 1:10)

DimPlot(pbmc,reduction = "umap",  label=T,label.size = 8) + theme(text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

FeaturePlot(pbmc, c("Ctsq","Prl2c5","Cts8", "Prl6a1"),min.cutoff = "q01", max.cutoff = "q90")


pbmc <- ScaleData(object = pbmc)

pbmc <- RunPCA(pbmc,verbose = F)
DimPlot(object = pbmc)
ElbowPlot(object = pbmc)
DimPlot(object = pbmc, reduction = "pca", label=T)
pbmc <- RunTSNE(object = pbmc, dims = 1:7)
DimPlot(object = pbmc, label=T)


pbmc <- FindNeighbors(object = pbmc, dims = 1:7)
pbmc <- FindClusters(object = pbmc, algorithm = 1, resolution =0.13)
DimPlot(pbmc,  label=T,label.size = 8) + theme(text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))



pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)

         
         
        
        DimPlot(object = pbmc,label=T,label.size = 8)
        pbmc <- SubsetData(object = pbmc2, ident.use = 1)
        write.csv(pbmc@meta.data, "mpl12/dE07.5.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 2)
        write.csv(pbmc@meta.data, "mpl12/dE08.5.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 3)
        write.csv(pbmc@meta.data, "mpl12/dE09.5.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 4)
        write.csv(pbmc@meta.data, "mpl12/dE09.5d.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 5)
        write.csv(pbmc@meta.data, "mpl12/dE10.5.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 6)
        write.csv(pbmc@meta.data, "mpl12/dE10.5h.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 7)
        write.csv(pbmc@meta.data, "mpl12/dE11.5h.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 8)
        write.csv(pbmc@meta.data, "mpl12/dE12.5.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 9)
        write.csv(pbmc@meta.data, "mpl12/dE12.5h.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 10)
        write.csv(pbmc@meta.data, "mpl12/dE13.5.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 11)
        write.csv(pbmc@meta.data, "mpl12/dE14.5.csv")
        pbmc <- SubsetData(object = pbmc2, ident.use = 12)
        write.csv(pbmc@meta.data, "mpl12/dE14.5d.csv")
        
        
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =1)])
        write.csv(pbmc.data, "mpl12/E07.5.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =2)])
        write.csv(pbmc.data, "mpl12/E08.5.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =3)])
        write.csv(pbmc.data, "mpl12/E09.5.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =4)])
        write.csv(pbmc.data, "mpl12/E09.5d.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =5)])
        write.csv(pbmc.data, "mpl12/E10.5.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =6)])
        write.csv(pbmc.data, "mpl12/E10.5h.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =7)])
        write.csv(pbmc.data, "mpl12/E11.5h.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =8)])
        write.csv(pbmc.data, "mpl12/E12.5.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =9)])
        write.csv(pbmc.data, "mpl12/E12.5h.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =10)])
        write.csv(pbmc.data, "mpl12/E13.5.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =11)])
        write.csv(pbmc.data, "mpl12/E14.5.csv")
        pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =12)])
        write.csv(pbmc.data, "mpl12/E14.5d.csv")
        
        
        
        pbmc.data=read.csv("mpl12/E07.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE07.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        pbmc <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E07.5")
        pbmc.data=read.csv("mpl12/E08.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE08.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        E08.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E08.5")
        pbmc.data=read.csv("mpl12/E09.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE09.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        E09.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E09.5")
        pbmc.data=read.csv("mpl12/E10.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE10.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        E10.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E10.5")
        pbmc.data=read.csv("mpl12/E10.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE10.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
        E10.5h <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E10.5h")
        pbmc.data=read.csv("mpl12/E11.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE11.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
        E11.5h <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E11.5h")
        pbmc.data=read.csv("mpl12/E12.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE12.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        E12.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E12.5")
        pbmc.data=read.csv("mpl12/E12.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE12.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
        E12.5h <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E12.5h")
        pbmc.data=read.csv("mpl12/E13.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE13.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        E13.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E13.5")
        pbmc.data=read.csv("mpl12/E14.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        meta.data=read.csv("mpl12/dE14.5.csv",check.names = FALSE, header=TRUE,row.names=1)
        E14.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E14.5")
        
        pbmc <- merge(x = E07.5, y = c(E08.5,E09.5,E10.5, E10.5h, E11.5h, E12.5, E12.5h, E13.5,E14.5 ),  
                      add.cell.ids = c("E07.5","E08.5","E09.5","E10.5", "E10.5h", "E11.5h", "E12.5", "E12.5h", "E13.5","E14.5"), project = "mpl")
        table(Idents(object = pbmc))
        
         
        pbmc <- NormalizeData(object = pbmc)
        # pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures =1200, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
          pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
        length(x = VariableFeatures(object = pbmc))
        pbmc <- ScaleData(object = pbmc)
        pbmc <- RunPCA(pbmc,verbose = F)
        DimPlot(object = pbmc, reduction = "pca", label=T,label.size = 8)+ theme(text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))
        ElbowPlot(object = pbmc)
        pbmc <- RunUMAP(pbmc, reduction.use = "pca", dims = 1:7)
        DimPlot(pbmc,reduction = "umap",  label=T,label.size = 8) + theme(text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))
        
        pbmc <- RunTSNE(object = pbmc, dims = 1:6)
        DimPlot(object = pbmc, reduction = "tsne", label=T,label.size = 8)+ theme(text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))
        
        FeaturePlot(pbmc, c("Cdx2"),min.cutoff = "q01", max.cutoff = "q90") 
        FeaturePlot(pbmc, c("Hand1","Cdx2"),min.cutoff = "q01", max.cutoff = "q90") 
        FeaturePlot(pbmc, c("Crip2", "Car2"),min.cutoff = "q01", max.cutoff = "q90",blend = T)
        FeaturePlot(pbmc, c("Hand1", "Car2"),min.cutoff = "q01", max.cutoff = "q90",blend = T)
        FeaturePlot(pbmc, c("Phlda2"),min.cutoff = "q01", max.cutoff = "q90")
        
       
            