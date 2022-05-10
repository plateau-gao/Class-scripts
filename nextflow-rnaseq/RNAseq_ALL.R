# ============================================================================ #
# LOADING PACKAGES
# ---------------------------------------------------------------------------- #
library(dplyr)
library(stringr)
library(DESeq2)
library(tximport)
library(tximportData)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(gprofiler2)
library(RColorBrewer)
# ============================================================================ #


# ============================================================================ #
# SETTING INPUT FILES
# ---------------------------------------------------------------------------- #
workdir = "/Users/plateau/Desktop/5632_ca2/result/GSE112593"
outdir = "All_cells"
input = "pipeline_result"

#File Path:
path_genetpm = file.path(input, "salmon.merged.gene_tpm.tsv" )
path_genecount = file.path(input,"salmon.merged.gene_counts.tsv")
path_tx2gene = file.path(input,"salmon_tx2gene.tsv")
path_metadata = file.path(input, "metadata.csv")

setwd(workdir)
# ============================================================================ #



# ============================================================================ #
# Differential expression analysis
# 1. metadata preparation
# 2. select gene with differenctial expression
# 3. Visuralization
# ---------------------------------------------------------------------------- #
# metadata preparation
metadata = read.csv(path_metadata, header = TRUE)
metadata = dplyr::mutate(metadata, 
                         path = file.path(input,
                                          "salmon",
                                          sample,
                                          "quant.sf"))
metadata$condition = factor(metadata$condition, levels = c("wt", "yap_KO"))
metadata$cell_type = factor(metadata$cell_type, levels = c("naive", "treg_unstm", "treg_stm"))
files = metadata$path
names(files) = metadata$sample
colnames(metadata) = c("sample", "genetype", "cell_type", "path")

# Select gene with differential expression
t2g = read.csv(path_tx2gene, sep = "\t", header = FALSE)
txi = tximport(files, type = "salmon", tx2gene = t2g)
dds = DESeqDataSetFromTximport(txi, metadata, ~ cell_type + genetype )
dds = dds[rowSums(DESeq2::counts(dds)) > 1, ]
dds = DESeq(dds)

res_yap = results(dds, contrast = c("condition", 'yap_KO', 'wt'))
summary(res_yap)

yap_ko = res_yap[!is.na(res_yap$padj), ]
yap_ko = yap_ko[yap_ko$padj < 0.1, ]
summary(yap_ko)

yap_log2fc = yap_ko[abs(yap_ko$log2FoldChange) > 1, ]
summary(yap_log2fc)

# Visualization
# 1. MA-plot
gene_names = read.csv(path_genecount, sep = "\t", row.names = "gene_id")
gene_names = gene_names[rownames(res_yap), ]
res_yap_name = res_yap
rownames(res_yap_name) = gene_names$gene_name

png(file.path(outdir, "MA-plot.png") ,
    width = 750, height = 600, units = "px")

ggmaplot(res_yap_name,
         main = "Different gene expression between Yap_KO and wide type",
         fdr = 0.05, fc = 2, size = 2,
         legend = "top", top = 10,
         xlab = "log2 mean expression",
         ylab = "log2 fold change",
         palette = c("#B31B21", "#1465AC", "lightgray"),
         font.label = c("bold", 15),
         font.legend = c("bold", 15),
         font.main = c("bold",15),
         ggtheme = theme_minimal())
dev.off()

# 2. Heat Map
### select gene with log2 fold change > 1
countsNormalize = DESeq2::counts(dds, normalize = TRUE)
gene_lfc1 = countsNormalize[row.names(yap_log2fc), ]
gene_names = read.csv(path_genecount, sep = "\t", row.names = "gene_id")
gene_names = gene_names[row.names(gene_lfc1), ]
row.names(gene_lfc1) = gene_names$gene_name
### draw heat map
heatmap = pheatmap(log2(gene_lfc1 + 1),
         show_rownames = TRUE,
         color = colorRampPalette(c("#1465AC", "#FFFFFF", "#B31B21"))(50),
         clustering_distance_rows = "euclidean", 
         fontsize_row = 12,
         scale = "row",
         filename = file.path(outdir, "heatmap.png"),
         width = 7,
         height = 16,
         main = "Gene expression level",
         )

# 3. Bar plot
### prepare data
yap_log2fc_name = yap_log2fc
gene_names = read.csv(path_genecount, sep = "\t", row.names = "gene_id")
gene_names = gene_names[row.names(yap_log2fc), ]
row.names(yap_log2fc_name) = gene_names$gene_name

df_yap_log2fc = as.data.frame(yap_log2fc_name)
df_yap_log2fc$gene_name = row.names(df_yap_log2fc)
df_yap_log2fc = df_yap_log2fc[order(df_yap_log2fc$log2FoldChange), ]
df_yap_log2fc$group = ifelse (df_yap_log2fc$log2FoldChange > 0, "up", "down")
for (i in seq_len(NROW(df_yap_log2fc))) {
  if (df_yap_log2fc$log2FoldChange[i] < 0) {
    df_yap_log2fc$log2FoldChange[i] = abs(df_yap_log2fc$log2FoldChange[i])
  }
}
df_yap_log2fc$gene_name = factor(df_yap_log2fc$gene_name,
                                 levels = df_yap_log2fc$gene_name)
### draw bar plot
png(file.path(outdir, "hits.png") ,
    width = 500, height = 1600, units = "px")
hist_lfc = ggplot(data = df_yap_log2fc, 
                  aes(x = gene_name,
                      y = log2FoldChange,
                      fill = group)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c( "#1465AC", "#B31B21")) +
  theme_minimal() + 
  coord_flip() +
  ggtitle("Log2 fold change (abs>1) of Yap_KO / Wide Type T cells") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 14, face = "bold"))
hist_lfc
dev.off()
# ============================================================================ #


# ============================================================================ #
# Functional Enrichment analysis
# 1. GO term analysis
# 2. KEGG analysis
# ---------------------------------------------------------------------------- #
# GO analysis
#### combine analysis
goResult_lfc = gost(query = rownames(yap_log2fc),
                organism = "mmusculus",
                significant = TRUE,
                sources = "GO")
df_goResult_lfc = goResult_lfc$result[, c(1:13)]
write.csv(df_goResult_lfc, file = file.path(outdir, "goResult_lfc.csv"))

goResult = gost(query = rownames(yap_ko),
                    organism = "mmusculus",
                    significant = TRUE,
                    sources = "GO")
df_goResult = goResult$result[, c(1:13)]
write.csv(df_goResult, file = file.path(outdir, "goResult.csv"))


#### GO analysis of up-regulated genes
goResult_up_lfc = gost(query = rownames(yap_log2fc[yap_log2fc$log2FoldChange>0, ]),
                organism = "mmusculus",
                significant = TRUE,
                sources = "GO")
df_goResult_up_lfc = goResult_up_lfc$result[, c(1:13)]
write.csv(df_goResult_up_lfc, file = file.path(outdir, "goResult_up_lfc.csv"))

goResult_up = gost(query = rownames(yap_ko[yap_ko$log2FoldChange>0, ]),
                       organism = "mmusculus",
                       significant = TRUE,
                       sources = "GO")
df_goResult_up = goResult_up$result[, c(1:13)]
write.csv(df_goResult_up, file = file.path(outdir, "goResult_up.csv"))



#### GO analysis of down-regulated genes
goResult_down_lfc = gost(query = rownames(yap_log2fc[yap_log2fc$log2FoldChange<0, ]),
                   organism = "mmusculus",
                   significant = TRUE,
                   sources = "GO")
df_goResult_down_lfc = goResult_down_lfc$result[, c(1:13)]
write.csv(df_goResult_down_lfc, file = file.path(outdir, "goResult_down_lfc.csv"))

goResult_down= gost(query = rownames(yap_ko[yap_ko$log2FoldChange<0, ]),
                         organism = "mmusculus",
                         significant = TRUE,
                         sources = "GO")
df_goResult_down = goResult_down$result[, c(1:13)]
write.csv(df_goResult_down, file = file.path(outdir, "goResult_down.csv"))



# KEGG analysis
keggResult_up_lfc = gost(query = rownames(yap_log2fc[yap_log2fc$log2FoldChange>0, ]),
                  organism = "mmusculus",
                  significant = TRUE,
                  sources = "KEGG")
keggResult_up_lfc$result
df_keggResult_up_lfc = keggResult_up_lfc$result[, c(1:13)]
write.csv(df_keggResult_up_lfc, file = file.path(outdir, "keggResult_up_lfc.csv"))

keggResult_up = gost(query = rownames(yap_ko[yap_ko$log2FoldChange>0, ]),
                         organism = "mmusculus",
                         significant = TRUE,
                         sources = "KEGG")
keggResult_up$result
df_keggResult_up = keggResult_up$result[, c(1:13)]
write.csv(df_keggResult_up, file = file.path(outdir, "keggResult_up.csv"))

