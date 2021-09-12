suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(clustree))
suppressMessages(library(varhandle))
suppressMessages(library(DoubletFinder))
suppressMessages(library(gridExtra))

sample <- Read10X("~/samples/single_cell_samples/SRR8526547/") %>% CreateSeuratObject(., project="SRR8526547")
sample <- AddMetaData(sample, PercentageFeatureSet(sample, pattern = "^MT-"), col.name= "percent.mt")

print(VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.3))

plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)

sample <- subset(sample, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 10)

sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sample)
sample <- ScaleData(sample, features = all.genes)
sample <- RunPCA(sample, features = VariableFeatures(object = sample))
ElbowPlot(sample)+labs(title= "SRR8526547")


### quantitative approach to choose PC cut-off:
# We can calculate where the principal components start to elbow by taking the larger value of:

#     The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
#     The point where the percent change in variation between the consecutive PCs is less than 0.1%. (I changed it to 0.05%)
#     Finally, we take the minimum between the two values
pct <- sample@reductions$pca@stdev / sum(sample@reductions$pca@stdev) * 100
cum <- cumsum(pct)
co1 <- which(cum > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.05),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
co2
pcs <- min(co1, co2)
pcs
### The 'quantitative approach' say that the best PC cut-off to use is PC17, and the elbow plot says so
### this is why I stayed on PC17
### I will some higher and lower values, just to verify things

## RunUMAP and RunTSNE
sample <- RunUMAP(sample, dims=1:17)
sample <- RunTSNE(sample, dims=1:17)
sample <- FindNeighbors(sample, dims = 1:17) %>% FindClusters(., resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
Idents(sample) <- "RNA_snn_res.1"
print(DimPlot(sample, label=TRUE, reduction="umap", pt.size=0.3)+labs(title=paste0("SRR8526547", " - PCs: 1:",17,", resolution: ",1)))

