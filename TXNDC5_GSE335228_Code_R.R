library(Seurat)
library(ggplot2)
library(scales)
######################################################
seu_anno <- readRDS("data/seu_anno.rds")
###################################################### Fig1 K

DimPlot(
  seu_anno,
  reduction = "umap.harmony",
  group.by = "celltype_fine",
  label = TRUE
)
DimPlot(
  seu_anno,
  reduction = "umap.harmony",
  group.by = "group",
  label = TRUE
)
######################################################Fig1 L
df <- FetchData(
  seu_anno,
  vars = c("TXNDC5", "celltype_fine")
)
cell_order <- names(sort(
  tapply(df$TXNDC5, df$celltype_fine, mean, na.rm = TRUE),
  decreasing = TRUE
))

df$celltype_fine <- factor(
  df$celltype_fine,
  levels = cell_order
)
my_colors <- hue_pal()(length(cell_order))
names(my_colors) <- cell_order
p <- ggplot(
  df,
  aes(
    x = celltype_fine,
    y = TXNDC5,
    fill = celltype_fine
  )
) +
  geom_violin(
    scale = "width",
    width = 0.85,
    bw = 0.1,          
    trim = FALSE,     
    color = "black",
    linewidth = 0.45
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2.5,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  scale_fill_manual(values = my_colors) +
  coord_cartesian(
    ylim = c(0, 3.5),
    expand = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(0, 3.5, by = 0.5)
  ) +
  labs(
    x = NULL,
    y = "TXNDC5"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 13,
      color = "black"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

print(p)

ggsave(
  "TXN_umap_human.png",dpi = 600,
  p,
  width =6.09,
  height = 5.33
)
######################################################Fig1 M
DefaultAssay(seu_endo) <- "RNA"


# Seurat v5
counts <- GetAssayData(
  seu_endo,
  assay = "RNA",
  layer = "counts"
)
meta <- seu_endo@meta.data[colnames(counts), , drop = FALSE]
check_group <- meta %>%
  group_by(.data[[sample_col]]) %>%
  summarise(n_group = n_distinct(.data[[group_col]]))

check_group
sample_factor <- factor(meta[[sample_col]])

mm <- sparse.model.matrix(~ 0 + sample_factor)

colnames(mm) <- levels(sample_factor)

pb_counts <- counts %*% mm
sample_info <- meta %>%
  select(
    sample = all_of(sample_col),
    group  = all_of(group_col)
  ) %>%
  distinct()

sample_info <- sample_info[
  match(colnames(pb_counts), sample_info$sample),
]

sample_info


y <- DGEList(
  counts = pb_counts,
  samples = sample_info
)

y <- calcNormFactors(y)

cpm_mat <- cpm(y, log = FALSE)
logcpm_mat <- cpm(y, log = TRUE, prior.count = 1)
txn_pseudobulk <- data.frame(
  sample = colnames(logcpm_mat),
  TXNDC5_logCPM = as.numeric(logcpm_mat["TXNDC5", ])
)

txn_pseudobulk <- txn_pseudobulk %>%
  left_join(sample_info, by = "sample")

txn_pseudobulk



library(ggplot2)
t.test(
  TXNDC5_logCPM ~ group,
  data = txn_pseudobulk
)
p <- ggplot(txn_pseudobulk,
            aes(x = group, y = TXNDC5_logCPM, fill = group, color = group)) +
  geom_boxplot(width = 0.55, linewidth = 0.8, alpha = 0.35,
               outlier.shape = NA) +
  geom_jitter(width = 0.10, size = 3.5, alpha = 0.9) +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 3.5,
               fill = "white", color = "black") +
  
  annotate("segment",
           x = 1, xend = 2,
           y = 7.2, yend = 7.2,
           linewidth = 0.7, color = "black") +
  annotate("segment",
           x = 1, xend = 1,
           y = 7.0, yend = 7.2,
           linewidth = 0.7, color = "black") +
  annotate("segment",
           x = 2, xend = 2,
           y = 7.0, yend = 7.2,
           linewidth = 0.7, color = "black") +
  
  # P-value
  annotate("text",
           x = 1.5, y = 7.45,
           label = "P = 0.010960",###pvalue=0.010960 after t.test(TXNDC5_logCPM ~ group,data = txn_pseudobulk)
           size = 5) +
  
  scale_y_continuous(
    limits = c(0, 8),
    expand = expansion(mult = c(0, 0.02))
  )+ scale_fill_manual(
    values = c("Control" = "#689EFE", "Disease" = "#F56F70")
  ) +
  scale_color_manual(
    values = c("Control" = "#689EFE", "Disease" = "#F56F70")
  ) +
  labs(
    x = NULL,
    y = expression("TXNDC5 pseudobulk expression (log"[2]*"CPM)")
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.7)
  )

p
ggsave(
  filename = "TXNDC5_pseudobulk_GSE335228_totalEC.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 1200
)

######################################################FigS9 E,H
seu_endo <- NormalizeData(seu_endo)
seu_endo <- FindVariableFeatures(seu_endo, nfeatures = 3000)
seu_endo <- ScaleData(seu_endo)
seu_endo <- RunPCA(seu_endo, npcs = 30)
seu_endo <- RunHarmony(
  seu_endo,
  group.by.vars = "orig.ident"
)
ElbowPlot(seu_endo, ndims = 30)
#i = 15, the resulting UMAP (Fig. S9H) shows a clustering pattern similar to 6 clusters structure observed in GSE273062.
#i = 10, the resulting UMAP (Fig. S9E) shows a clustering pattern similar to 3 clusters structure observed in GSE273062.
seu_endo <- FindNeighbors(seu_endo, reduction = "harmony", dims = 1:i)
seu_endo <- FindClusters(seu_endo, resolution = 0.2)
seu_endo <- RunUMAP(seu_endo, reduction = "harmony", dims = 1:i)
DimPlot(seu_endo, reduction = "umap", label = TRUE,pt.size = 1,group.by = "endo_annotation")
DimPlot(seu_endo, reduction = "umap", group.by = "group",pt.size = 2)

######################################################FigS9 F
genes <- c(
  "APLN", "CA4","EDNRB",
  "FTL", "FABP4","CTSB","TXNDC5", "BGN", "SNAI1"
)

p <- DotPlot(
  seu_endo,
  features = genes,
  group.by = "endo_annotation_a",
  dot.scale = 7
)
p$layers[[1]] <- geom_point(
  aes(size = pct.exp, fill = avg.exp.scaled),
  shape = 21,
  colour = "black",
  stroke = 0.4
)

p <- p +
  scale_fill_gradient(
    low = "#689EFE",
    high = "#F56F70",
    name = "Average Expression"
  ) +
  scale_radius(range = c(5, 12), name = "Percent Expressed") +
  guides(colour = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

p

TransEC_markers <- FindMarkers(
  object = seu_endo_filtered,
  group.by = "celltype",
  ident.1 = "Trans ECs",
  ident.2 = NULL,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

ggsave(
  "endo_dotplot.png",
  plot = p,
  width = 9.29,
  height = 2.28,
  units = "in",
  bg = "white",dpi = 600
)


######################################################FigS9 G
### ECM producing ECs
seu_ecmec <- subset(
  seu_endo,
  subset = endo_annotation_a == "ECM producing ECs")
sample_col <- "orig.ident"  
group_col  <- "group"

DefaultAssay(seu_ecmec) <- "RNA"


# Seurat v5
counts <- GetAssayData(
  seu_ecmec,
  assay = "RNA",
  layer = "counts"
)


meta <- seu_ecmec@meta.data[colnames(counts), , drop = FALSE]


check_group <- meta %>%
  group_by(.data[[sample_col]]) %>%
  summarise(n_group = n_distinct(.data[[group_col]]))

check_group
sample_factor <- factor(meta[[sample_col]])

mm <- sparse.model.matrix(~ 0 + sample_factor)

colnames(mm) <- levels(sample_factor)

pb_counts <- counts %*% mm
sample_info <- meta %>%
  select(
    sample = all_of(sample_col),
    group  = all_of(group_col)
  ) %>%
  distinct()

sample_info <- sample_info[
  match(colnames(pb_counts), sample_info$sample),
]

sample_info


y <- DGEList(
  counts = pb_counts,
  samples = sample_info
)

y <- calcNormFactors(y)

cpm_mat <- cpm(y, log = FALSE)
logcpm_mat <- cpm(y, log = TRUE, prior.count = 1)
txn_pseudobulk <- data.frame(
  sample = colnames(logcpm_mat),
  TXNDC5_logCPM = as.numeric(logcpm_mat["TXNDC5", ])
)

txn_pseudobulk <- txn_pseudobulk %>%
  left_join(sample_info, by = "sample")

txn_pseudobulk



library(ggplot2)
t.test(
  TXNDC5_logCPM ~ group,
  data = txn_pseudobulk
)
p <- ggplot(txn_pseudobulk,
            aes(x = group, y = TXNDC5_logCPM, fill = group, color = group)) +
  geom_boxplot(width = 0.55, linewidth = 0.8, alpha = 0.35,
               outlier.shape = NA) +
  geom_jitter(width = 0.10, size = 3.5, alpha = 0.9) +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 3.5,
               fill = "white", color = "black") +
  annotate("segment",
           x = 1, xend = 2,
           y = 7.2, yend = 7.2,
           linewidth = 0.7, color = "black") +
  annotate("segment",
           x = 1, xend = 1,
           y = 7.0, yend = 7.2,
           linewidth = 0.7, color = "black") +
  annotate("segment",
           x = 2, xend = 2,
           y = 7.0, yend = 7.2,
           linewidth = 0.7, color = "black") +
  
  # Pvalue
  annotate("text",
           x = 1.5, y = 7.45,
           label = "P = 0.038417",###pvalue=0.038417 after t.test(TXNDC5_logCPM ~ group,data = txn_pseudobulk)
           size = 5) +
  
  scale_y_continuous(
    limits = c(0, 8),
    expand = expansion(mult = c(0, 0.02))
  )+ scale_fill_manual(
    values = c("Control" = "#689EFE", "Disease" = "#F56F70")
  ) +
  scale_color_manual(
    values = c("Control" = "#689EFE", "Disease" = "#F56F70")
  ) +
  labs(
    x = NULL,
    y = expression("TXNDC5 pseudobulk expression (log"[2]*"CPM)")
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.7)
  )

p
ggsave(
  filename = "TXNDC5_pseudobulk_GSE335228_totalEC.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 1200
)


######################################################FigS9 I

#1. Contingency Table of Raw Cell Counts

tab <- with(
  seu_endo@meta.data,
  table(celltype, endo_annotation)
)

tab <- tab[
  rowSums(tab) > 0,
  colSums(tab) > 0,
  drop = FALSE
]

if (sum(tab) == 0) {
  stop("No valid cells have both annotations.")
}

# 2. Calculate Expected Counts Under the Assumption of Independence
expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
dimnames(expected) <- dimnames(tab)

# 3. Calculate O/E and log2(O/E)
oe <- tab / expected
log2_oe <- log2(oe)

# 4. Prepare Data for Plotting and Retain Observed and Expected Counts
df <- as.data.frame(tab)
names(df) <- c("celltype", "endo_annotation", "n")

df$expected <- as.vector(expected)
df$OE <- as.vector(oe)
df$log2_OE <- as.vector(log2_oe)
lim <- max(1, abs(df$log2_OE[is.finite(df$log2_OE)]))
df$fill_value <- pmax(-lim, pmin(lim, df$log2_OE))
df$label <- ifelse(
  df$n == 0,
  "-Inf",
  sprintf("%.2f", df$log2_OE)
)
p <- ggplot(
  df,
  aes(
    x = endo_annotation,
    y = celltype,
    fill = fill_value
  )
) +
  geom_tile(colour = "black", linewidth = 0.5) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(
    low = "#FFF8F0",
    mid = "#F2C14E",
    high = "#F56F70",
    midpoint = -1,
    limits = c(-lim, 3),
    name = expression(log[2](O/E))
  ) +
  labs(
    x = NULL,
    y = NULL,
    caption = "Values: log2(O/E); -Inf: zero observed cells."
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p


library(dplyr)
library(ggplot2)

df_raw <- df

df <- df %>%
  group_by(celltype) %>%
  mutate(
    fill_value = {
      x <- fill_value
      x[!is.finite(x)] <- NA_real_
      s <- sd(x, na.rm = TRUE)
      
      if (is.finite(s) && s > 0) {
        (x - mean(x, na.rm = TRUE)) / s
      } else {
        ifelse(is.na(x), NA_real_, 0)
      }
    },
    label = if_else(
      is.na(fill_value),
      "NA",
      sprintf("%.2f", fill_value)
    )
  ) %>%
  ungroup()

p <- ggplot(
  df,
  aes(x = endo_annotation, y = celltype, fill = fill_value)
) +
  geom_tile(colour = "black", linewidth = 0.5) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(
    low = "#FFF8F0",
    mid = "#F2C14E",
    high = "#F56F70",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    na.value = "grey90",
    name = "Row Z-score"
  ) +
  labs(
    x = NULL,
    y = NULL,
    caption = "Values: row Z-score."
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p
ggsave(
  "cor_celltype.png",
  plot = p,
  width = 8,
  height = 5,
  units = "in",
  bg = "white",dpi = 600
)


######################################################FigS9 J
plot_data <- DotPlot(
  seu_endo,
  features = c("TXNDC5", "BGN","CCN2","FN1"),
  group.by = "endo_annotation"
)$data
p<-ggplot(plot_data, aes(x = features.plot, y = id, fill = avg.exp.scaled)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("#FFF8F0", "#F2C14E", "#F56F70")) +
  labs(x = NULL, y = NULL, fill = "Scaled expression") +
  theme_classic()
ggsave(
  "endo_hu_umapplot_5cluster.png",
  plot = p,
  width = 7.38,
  height = 5.65,
  units = "in",
  bg = "white",dpi = 600
)




