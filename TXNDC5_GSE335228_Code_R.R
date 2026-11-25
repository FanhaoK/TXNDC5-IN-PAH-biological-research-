library(Seurat)
setwd("F:/桌面/TXN返修/GSE209738_RAW")
path <- "F:/桌面/TXN返修/GSE249461_RAW"

Hypoxia  <- CreateSeuratObject(Read10X(file.path(path, "Hypoxia")),  project = "Hypoxia")
Normoixa <- CreateSeuratObject(Read10X(file.path(path, "Normoixa")), project = "Normoixa")
sce <- merge(Hypoxia, y = Normoixa, add.cell.ids = c("Hypoxia", "Normoixa"))
library(Seurat)
library(harmony)


library(harmony)


# QC
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern="^mt-")
sce <- subset(sce, subset=nFeature_RNA>500 & nFeature_RNA<7000 & percent.mt<25)

# 标准流程
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, nfeatures=2000)
sce <- ScaleData(sce, vars.to.regress="percent.mt")
sce <- RunPCA(sce, npcs=30)
ElbowPlot(sce, ndims = 30)
# Harmony
sce <- RunHarmony(sce, "orig.ident")

# 聚类+UMAP
sce <- FindNeighbors(sce, reduction="harmony", dims=1:15)
sce <- FindClusters(sce, resolution=0.5)
sce <- RunUMAP(sce, reduction="harmony", dims=1:15)

DimPlot(sce, label=T)
DimPlot(sce, group.by="orig.ident")
FeaturePlot(sce, features="Txndc5",order = T)
FeaturePlot(sce, features="Col1a2",order = T)
VlnPlot(sce, features = "Txndc5", sort = "increasing")
sce <- subset(sce, idents = c("10","11","12","13","14"), invert = TRUE)###这些群细胞比例小于2%，故判定为污染
sce <- subset(sce, idents = c("7","5"), invert = TRUE)###这些细胞是成纤维细胞污染
sce <- JoinLayers(sce)


marker <- FindAllMarkers(
  sce,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  min.diff.pct = 0.1
)


map <- c(
  "0"="Arterial EC", "1"="gCap", "2"="Arterial EC", "3"="gCap",
  "4"="aCap", "6"="Pulmonary venous EC",
  "8"="Lymphatic EC", "9"="aCap"
)

sce$celltype <- unname(map[as.character(sce$seurat_clusters)])
sce$celltype[is.na(sce$celltype)] <- "Unassigned"

Idents(sce) <- "celltype"

table(Idents(sce))




library(dplyr)
library(ggplot2)
library(patchwork)

df <- FetchData(sce, vars=c("Txndc5","celltype")) %>%
  group_by(celltype) %>%
  summarise(
    total = n(),
    positive = sum(Txndc5 > 0),
    percent = positive/total*100
  )

p1 <- ggplot(df, aes(reorder(celltype, positive), positive)) +
  geom_col(fill="#689EFE", color="black", linewidth=0.5) +
  coord_flip() +
  labs(x=NULL, y="Txndc5+ cells") +
  theme_classic()

p2 <- ggplot(df, aes(reorder(celltype, percent), percent)) +
  geom_col(fill="#F56F70", color="black", linewidth=0.5) +
  coord_flip() +
  labs(x=NULL, y="Txndc5+ cells (%)") +
  theme_classic()

p<-p1 + p2
ggsave(
  "Txndc5_positive_cells.png",
  plot = p,
  width = 6,
  height = 3,
  dpi = 600
)



library(dplyr)
library(ggplot2)

df <- FetchData(sce, c("Txndc5","celltype","orig.ident")) %>%
  group_by(celltype, orig.ident) %>%
  summarise(
    total = n(),
    positive = sum(Txndc5 > 0),
    percent = positive / total * 100,
    .groups = "drop"
  )

# 设置组别顺序
df$orig.ident <- factor(df$orig.ident, levels=c("Normoixa","Hypoxia"))

p <- ggplot(df, aes(x=celltype, y=percent, fill=orig.ident)) +
  geom_col(position=position_dodge(0.8), width=0.7,
           color="black", linewidth=0.4) +
  scale_fill_manual(values=c("Normoixa"="#689EFE",
                             "Hypoxia"="#F56F70")) +
  labs(x=NULL, y="Txndc5+ cells (%)", fill=NULL) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

p





custom_cmap = mcolors.LinearSegmentedColormap.from_list(
  "smooth_gray_red",
  [
    "lightgray",   # low expression
    "#E8D7DA",
    "#F0BFC3",
    "#F59A9F",
    "#F77E82",
    "#F56F70"      # high expression
  ]
)

library(data.table)



data <- read.csv(
  "GSE229561_Processed_data_files.csv",
  row.names = 1,
  check.names = FALSE
)

barcodes <- data[,1]
rownames(data)<-barcodes
head(barcodes)

counts <- t(as.matrix(data))

seu <- CreateSeuratObject(
  counts = counts,
  project = "GSE229561",
  min.cells = 3,
  min.features = 200
)
seu$orig.ident <- sub(
  ".*_",
  "",
  colnames(seu)
)



library(Seurat)
library(harmony)

# 标准化
seu <- NormalizeData(seu)

# 高变基因
seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 3000
)

# Scale + PCA
seu <- ScaleData(seu)
seu <- RunPCA(seu)

# Harmony去批次（按orig.ident）
seu <- RunHarmony(
  seu,
  group.by.vars = "orig.ident",
  dims.use = 1:30
)

# 聚类
seu <- FindNeighbors(
  seu,
  reduction = "harmony",
  dims = 1:30
)

seu <- FindClusters(
  seu,
  resolution = 0.5
)

# UMAP
seu <- RunUMAP(
  seu,
  reduction = "harmony",
  dims = 1:30
)

# 查看
DimPlot(seu, reduction = "umap", label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "orig.ident")
FeaturePlot(
  seu,
  features = "Pecam1",
  reduction = "umap",order = T
)












