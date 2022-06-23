setwd("C:/Users/Admin/Documents/Computational Biology exercise/Computational Biology exercise")

suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

#Read the data to Seurat object 
my_data <- Read10X("WT")

WT_cells <- CreateSeuratObject(counts = my_data, project = "WT")


my_data <- Read10X("AD/")

AD_cells <- CreateSeuratObject(counts = my_data, project = "AD")


all_cells_merge <- merge(WT_cells, y = AD_cells, add.cell.ids = c("WT", "AD"))


#Visualize quality control metrics. How many cells do you have from each data set?

#calculates the percentage of mitochondrial genes
all_cells_merge[["percent.mt"]] <- PercentageFeatureSet(all_cells_merge, pattern = "^mt-")

VlnPlot(all_cells_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

table(all_cells_merge@meta.data$orig.ident)

#From WT cell data set there are 4581 cells, and from AD there are 5103 cells.

#removing unwanted cells from the dataset: filter cells that have unique feature
#counts over 2,500 or less than 200 and filter cells that have >5% mitochondrial counts.

all_cells_merge <- subset(all_cells_merge, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalize the data using “LogNormalize”.


all_cells_merge <- NormalizeData(all_cells_merge, normalization.method = "LogNormalize")

#Identify highly variable features (feature selection).

all_cells_merge <- FindVariableFeatures(all_cells_merge, selection.method = "vst", nfeatures = 2000)
#find to 10 most highly variable genes
top10 <- head(VariableFeatures(all_cells_merge), 10)

#the 10 most highly variable genes: "Ngp" ,"Chil3", "Camp", "Dlk1", "Ltf", "Lcn2", "Retnlg", "H2-Eb1", "H2-Aa", "H2-Ab1".

#Plot variable features with labels.

plot1 <- VariableFeaturePlot(all_cells_merge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scaling the data

all.genes <- rownames(all_cells_merge)
all_cells_merge <- ScaleData(all_cells_merge, features = all.genes)

#Perform linear dimensional reduction.

all_cells_merge <- RunPCA(all_cells_merge, features = VariableFeatures(object = all_cells_merge))

#Cluster the cells using a graph-based clustering approach.

all_cells_merge <- FindNeighbors(all_cells_merge, dims = 1:10)
all_cells_merge <- FindClusters(all_cells_merge, resolution = 0.5)

#Run non-linear dimensional reduction (UMAP) and plot the results.

all_cells_merge <- RunUMAP(all_cells_merge, dims = 1:10)

DimPlot(all_cells_merge, reduction = "umap")

#Plot the origin of the cells (AD or WT) over the UAMP maps.

DimPlot(all_cells_merge, reduction = "umap", group.by = "orig.ident")

#Find the most differential clusters between the AD sample and the WT sample. Plot
#for each cluster the proportion of cells originating from the AD mice and the WT mice.

clusterdiff<-rep(NA, 12)

Idents(all_cells_merge) <- all_cells_merge$seurat_clusters
DefaultAssay(all_cells_merge) <- "RNA"
for(i in 0:11){
deg.cluster <- FindMarkers(all_cells_merge, ident.1 = "WT", ident.2 = "AD", verbose = TRUE, group.by="orig.ident", subset.ident = i, logfc.threshold = log(2))
clusterdiff[i+1]<-nrow(deg.cluster)
}

#The differential clusters between the AD sample and the WT sample are cluster 7,9 and 11.

#for each cluster the proportion of cells originating from the AD mice and the WT mice.

cluster7<-filter(all_cells_merge@meta.data, seurat_clusters==7)

barplot(prop.table(table(cluster7$orig.ident))*100, ylim =  c(0,60), main= "cluster 7")

cluster9<-filter(all_cells_merge@meta.data, seurat_clusters==9)

barplot(prop.table(table(cluster9$orig.ident))*100, ylim =  c(0,60), main= "cluster 9")

cluster11<-filter(all_cells_merge@meta.data, seurat_clusters==11)

barplot(prop.table(table(cluster11$orig.ident))*100, ylim =  c(0,60), main= "cluster 11")


#The top differential marker genes for the most differential clusters.

cluster7_markers<- FindMarkers(all_cells_merge, ident.1 = "WT", ident.2 = "AD", verbose = TRUE, group.by="orig.ident", subset.ident = 7, logfc.threshold = log(2))

cluster9_markers<- FindMarkers(all_cells_merge, ident.1 = "WT", ident.2 = "AD", verbose = TRUE, group.by="orig.ident", subset.ident = 9, logfc.threshold = log(2))

cluster11_markers<- FindMarkers(all_cells_merge, ident.1 = "WT", ident.2 = "AD", verbose = TRUE, group.by="orig.ident", subset.ident = 11, logfc.threshold = log(2))


#plot the most differential marker in the most differential clusters.
VlnPlot(all_cells_merge, features = row.names(cluster7_markers)[1:3],split.by = "orig.ident",idents =7)

VlnPlot(all_cells_merge, features = row.names(cluster9_markers)[1:3],split.by = "orig.ident",idents =9)

VlnPlot(all_cells_merge, features = row.names(cluster11_markers)[1:3],split.by = "orig.ident",idents =11)
