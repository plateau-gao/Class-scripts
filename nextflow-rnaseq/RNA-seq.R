# ============================================================================ #
# LOADING PACKAGES
# ---------------------------------------------------------------------------- #
##### data manipulation 
library(tidyverse)

##### Differential expression analysis
library(DESeq2)
library(tximport)
library(apeglm)

##### Visualization
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)

library(PCAtools)
library(stats)
library(ggalt)
library(glmpca)

##### GO and KEGG analysis
library(clusterProfiler)
library(enrichplot)
library(ggridges)
library(GOplot)
library(topGO)
library(ggnewscale)
library(ggupset)
library(org.Hs.eg.db)

##### STRINGdb
library(STRINGdb)
# ============================================================================ #



# ============================================================================ #
# SETTING INPUT FILES
# ---------------------------------------------------------------------------- #
workdir = "/Users/plateau/Desktop/5632_ca2/result/GSE151090"
input = "nf-core_result"
output = "r_result"

#File Path:
path_tx2gene = file.path(input,"salmon_tx2gene.tsv")
path_sraruntable = "SraRunTable-2.csv"

setwd(workdir)
# ============================================================================ #



# ============================================================================ #
# TXIMPORT RESULT OF SALMON
# 1. prepare metadata
# 2. import analysis result of salmon
# ---------------------------------------------------------------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1. Prepare Metadata
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Metadata for all data
metadata = read.csv(path_sraruntable, header = TRUE) %>%
  dplyr::select(Age, timepoint, treatment) %>%
  tidyr::separate(timepoint, c("time", "v"), sep = " ")

metadata$Age = factor(metadata$Age, levels = unique(metadata$Age))
metadata$time = factor(metadata$time, levels = unique(metadata$time))
metadata$treatment = factor(metadata$treatment, levels = unique(metadata$treatment))
metadata$treatment = relevel(metadata$treatment, ref = "DMSO")
levels(metadata$Age) = c("P1", "P2", "P3")
levels(metadata$time) = c("24h", "48h", "4d")
colnames(metadata) = c("patient", "timepoint", "v", "treatment")

metadata = metadata %>% 
  dplyr::mutate(sample = paste0(patient, "_", timepoint, "_", treatment),
                path = file.path(input,
                                "salmon",
                                 paste0(patient, "_", timepoint, "_", treatment),
                                 "quant.sf")) %>%
  dplyr::select(c("sample", "patient", "timepoint", "treatment", "path"))

# metadata based on timepoint
metadata_24h = dplyr::filter(metadata, timepoint == "24h")
metadata_48h = dplyr::filter(metadata, timepoint == "48h")
metadata_4d  = dplyr::filter(metadata, timepoint == "4d")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 2. Import data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
t2g = read.csv(path_tx2gene, sep = "\t", header = FALSE)

files = metadata$path
names(files) = metadata$sample
txi = tximport(files = files, 
               type = "salmon",
               tx2gene = t2g)

files_24h = metadata_24h$path
names(files_24h) = metadata_24h$sample
txi_24h = tximport(files = files_24h,
                   type = "salmon",
                   tx2gene = t2g)

files_48h = metadata_48h$path
names(files_48h) = metadata_48h$sample
txi_48h = tximport(files = files_48h,
                   type = "salmon",
                   tx2gene = t2g)

files_4d = metadata_4d$path
names(files_4d) = metadata_4d$sample
txi_4d = tximport(files = files_4d,
                   type = "salmon",
                   tx2gene = t2g)
# ============================================================================ #



# ============================================================================ #
# DIFFERENT EXPRESSION ANALYSIS
# ---------------------------------------------------------------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1. Design formular
##### In order to benefit from the default settings of the package
##### Put the variable of interest at the end of the formula 
##### Make sure the control level is the first level.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
design = ~ patient + timepoint + treatment
design_time = ~ patient + treatment

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 2. create DESeqDataSet from tximport
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
dds = DESeqDataSetFromTximport(txi, metadata, design = design )
dds_24h = DESeqDataSetFromTximport(txi_24h, metadata_24h, design = design_time )
dds_48h = DESeqDataSetFromTximport(txi_48h, metadata_48h, design = design_time )
dds_4d  = DESeqDataSetFromTximport(txi_4d, metadata_4d, design = design_time )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 3. Pre-filter (counts > 0)
##### All dim: 34506 18
##### 24h dim: 29806 6
##### 48h dim: 30235 6
##### 4d  dim: 30187 6
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
dds = dds[rowSums(counts(dds)) > 0, ]
dds_24h = dds_24h[rowSums(counts(dds_24h)) > 0, ]
dds_48h = dds_48h[rowSums(counts(dds_48h)) > 0, ]
dds_4d = dds_4d[rowSums(counts(dds_4d)) > 0, ]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 4. Differential expression analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
dds = DESeq(dds)
res = results(dds, name = "treatment_VP_vs_DMSO")
resSig = subset(res, padj < 0.1 )

dds_24h = DESeq(dds_24h)
res_24h = results(dds_24h, name = "treatment_VP_vs_DMSO")
resSig_24h = subset(res_24h, padj < 0.1)

dds_48h = DESeq(dds_48h)
res_48h = results(dds_48h, name = "treatment_VP_vs_DMSO")
resSig_48h = subset(res_48h, padj < 0.1)

dds_4d  = DESeq(dds_4d)
res_4d = results(dds_4d, name = "treatment_VP_vs_DMSO")
resSig_4d = subset(res_4d, padj < 0.1)
# ============================================================================ #




# ============================================================================ #
# VISUALIZATION FOR DIFFERENT EXPRESSION ANALYSIS
# 1. MA-plot
# 2. PCA-plot
# 3. Heat map
# 4. Volcano-plot
# ---------------------------------------------------------------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1. MA-plot                      
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
res_name = data.frame(geneID = rownames(res))
colnames(t2g) = c("transcriptID","geneID", "gene_name")
res_name = dplyr::left_join(res_name, t2g, by = "geneID") %>%
  dplyr::select(geneID, gene_name) %>%
  dplyr::distinct()
resShrink = lfcShrink(dds, coef = "treatment_VP_vs_DMSO", type = "apeglm")

png(filename = file.path(output, "MA-VPvsDMSO.png"),
    width = 700, height = 700, units = "px")
ggmaplot(resShrink,
         main = "Differential expression between VP and DMSO",
         fdr = 0.05, fc = 2, size = 2,
         genenames = res_name$gene_name,
         legend = "top", top = 10,
         xlab = "Log2 mean expression",
         ylab = "log2 fold change",
         palette = c("#B31B21", "#1465AC", "lightgray"),
         font.label = c("bold", 12),
         label.rectangle = TRUE,
         font.legend = c("bold", 15),
         font.main = c("bold",15),
         ggtheme = theme_minimal())
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 2. PCA plot                     
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
vsd = vst(dds, blind = FALSE) %>% assay()
gpca = glmpca(vsd, L = 3)
gpca.data = gpca$factors %>%
  dplyr::mutate(patient = dds$patient,
         timepoint = dds$timepoint,
         treatment = dds$treatment)

cols = c("24h" = "#B31B21", "48h" = "#1465AC", "4d" = "#FBB731")
sizes = c("P1" = 2, "P2" = 6, "P3" = 10)
shapes = c("VP" = 15, "DMSO" = 16)
png(file.path(output, "PCA_plot.png") ,
    width = 500, height = 400 , units = "px") 
  ggplot(gpca.data, 
         aes(x = dim1, y = dim2,
             color = timepoint,
             size = patient,
             shape = treatment))  +
    geom_point() +
    scale_size_manual(values = sizes)  +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    labs(title = "Generalized principal component analysis",
         x = "PC1",
         y = "PC2") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text = element_blank())
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 3. Gene Clustering                  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# resSig: padj < 0.1
resSig_lfc = resSig[abs(resSig$log2FoldChange) > 1, ]
gene_sig = counts(dds)[rownames(resSig_lfc), ]
hclust_counts = gene_sig %>%
  t() %>%
  scale() %>%
  t()
ann_col = colData(dds)[, 2:4]
pheatmap(hclust_counts,
         show_rownames = FALSE,
         color = colorRampPalette(c("#1465AC", "#FFFFFF", "#B31B21"))(50),
         fontsize_row = 12,
         method = "average",
         annotation_col = as.data.frame(ann_col),
         filename = file.path(output, "heatmap_counts.png"),
         width = 7,
         height = 16,
         main = "Gene expression level")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 4. Volcano Plot                
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### Define function, which can largely simplify the code
markUpDown = function(dds_result) {
  df = as.data.frame(dds_result)
  df_filter = df[!is.na(df$padj), ]
  df_marked = df_filter %>%
    dplyr::mutate(gene_type = case_when(log2FoldChange > 1 & padj < 0.05 ~ "up",
                                        log2FoldChange < -1 & padj < 0.05 ~ "down",
                                        TRUE ~ "ns"))
  upNum   = nrow(subset(df_marked, df_marked$gene_type == "up"))
  downNum = nrow(subset(df_marked, df_marked$gene_type == "down"))
  return(list(df_mark = df_marked,
              up = upNum,
              down = downNum))
}

draw_vocalno = function(time, list_mark) {
  volcano = ggplot(list_mark$df_mark,
                       aes(x = log2FoldChange,
                           y = -log10(padj),
                           color = gene_type)) +
    geom_point(shape = 4) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") +
    geom_vline(xintercept = c(log2(0.5), log2(2)),
               linetype = "dashed") +
    scale_color_manual(values = cols,
                       breaks = c("up", "down", "ns"),
                       labels = c(paste0("Upregulated Significant genes = ", list_mark$up),
                                  paste0("Downregulated Significant genes = ", list_mark$down),
                                  "Non-significant genes")) +
    scale_y_continuous(limits = c(0, 14),
                       breaks = c(seq(0, 14, 2))) + 
    labs(x = "Log2 ( Fold Change )",
         y = "-Log10 ( adjust p-value ) ",
         title = paste0(time, ", VP/DMSO")) +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = c(0.3, 0.92))
  return(volcano)
}

##### Mark the upregulated or downregulated gene
list_marked_all = markUpDown(res)
list_marked_24h = markUpDown(res_24h)
list_marked_48h = markUpDown(res_48h)
list_marked_4d  = markUpDown(res_4d)

##### Draw volcano plot
cols = c("up" = "#B31B21", "down" = "#1465AC", "ns" = "gray")
volcano_all = draw_vocalno("All samples", list_marked_all)
volcano_24h = draw_vocalno("24h", list_marked_24h)
volcano_48h = draw_vocalno("48h", list_marked_48h)
volcano_4d  = draw_vocalno("4d", list_marked_4d)

png(filename = file.path(output, "volcano.png"),
    width = 1000, height = 1000, units = "px")
ggarrange(volcano_all, volcano_24h, volcano_48h, volcano_4d,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()
# ============================================================================ #




# ============================================================================ #
# GO AND KEGG ANALYSIS
# 1. data prepare
# 2. Over Representation analysis (GO, KEGG)
# 3. Gene set enrichment analysis (GO)
# 3. Visualization
# ---------------------------------------------------------------------------- #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data Prepare
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Database
GO_database = 'org.Hs.eg.db'
KEGG_database = 'hsa'

# 2. Gene of interest: padj < 0.05 and abs(log2FoldChange) > 1
df_res_24h_sig = dplyr::filter(list_marked_24h$df_mark, gene_type != "ns")
df_res_48h_sig = dplyr::filter(list_marked_48h$df_mark, gene_type != "ns")
df_res_4d_sig  = dplyr::filter(list_marked_4d$df_mark, gene_type != "ns" )


# 3. ID translate
idtran = function(x) {
  bitr(rownames(x),
       fromType = "ENSEMBL",
       toType = "ENTREZID",
       OrgDb = GO_database)
}
gene_trans_24h = idtran(df_res_24h_sig)
gene_trans_48h = idtran(df_res_48h_sig)
gene_trans_4d = idtran(df_res_4d_sig)

# 4. Prepare gene list for GSEA
genelist_24h       = res_24h$log2FoldChange
names(genelist_24h)= rownames(res_24h)
genelist_24h       = sort(genelist_24h, decreasing = TRUE)

genelist_48h       = res_48h$log2FoldChange
names(genelist_48h)= rownames(res_48h)
genelist_48h       = sort(genelist_48h, decreasing = TRUE)

genelist_4d        = res_4d$log2FoldChange
names(genelist_4d) = rownames(res_4d)
genelist_4d        = sort(genelist_4d, decreasing = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Over Representation analysis: GO & KEGG
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# def function
go_enrich_bp = function(x) {
  clusterProfiler::setReadable(
    clusterProfiler::simplify(
      enrichGO(x,
               OrgDb = GO_database,
               keyType = "ENSEMBL",
               ont = "BP", #BP, CC, MF, ALL
               pvalueCutoff = 0.01,
               qvalueCutoff = 0.05)
    ), OrgDb = GO_database
  )
}

go_enrich_cc = function(x) {
  clusterProfiler::setReadable(
    clusterProfiler::simplify(
      enrichGO(x,
               OrgDb = GO_database,
               keyType = "ENSEMBL",
               ont = "CC",
               pvalueCutoff = 0.01,
               qvalueCutoff = 0.05)
    ), OrgDb = GO_database
  )
}

go_enrich_mf = function(x) {
  clusterProfiler::setReadable(
    clusterProfiler::simplify(
      enrichGO(x,
               OrgDb = GO_database,
               keyType = "ENSEMBL",
               ont = "MF",
               pvalueCutoff = 0.01,
               qvalueCutoff = 0.05)
    ), OrgDb = GO_database
  )
}

kegg_enrich = function(x) {
  clusterProfiler::enrichKEGG(x,
                              organism = KEGG_database,
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)
}

# use function
eGO_24h_bp   = go_enrich_bp(rownames(df_res_24h_sig))
eGO_48h_bp   = go_enrich_bp(rownames(df_res_48h_sig))
eGO_4d_bp    = go_enrich_bp(rownames(df_res_4d_sig))

eGO_24h_mf   = go_enrich_mf(rownames(df_res_24h_sig))
eGO_48h_mf   = go_enrich_mf(rownames(df_res_48h_sig))
eGO_4d_mf    = go_enrich_mf(rownames(df_res_4d_sig))

eGO_24h_cc   = go_enrich_cc(rownames(df_res_24h_sig))
eGO_48h_cc   = go_enrich_cc(rownames(df_res_48h_sig))
eGO_4d_cc    = go_enrich_cc(rownames(df_res_4d_sig))

write.csv(eGO_24h_bp, file = file.path(output, "GO_enrichment_24h_bp.csv"))
write.csv(eGO_48h_bp, file = file.path(output, "GO_enrichment_48h_bp.csv"))
write.csv(eGO_4d_bp,  file = file.path(output, "GO_enrichment_4d_bp.csv"))

write.csv(eGO_24h_mf, file = file.path(output, "GO_enrichment_24h_mf.csv"))
write.csv(eGO_48h_mf, file = file.path(output, "GO_enrichment_48h_mf.csv"))
write.csv(eGO_4d_mf,  file = file.path(output, "GO_enrichment_4d_mf.csv"))

write.csv(eGO_24h_cc, file = file.path(output, "GO_enrichment_24h_cc.csv"))
write.csv(eGO_48h_cc, file = file.path(output, "GO_enrichment_48h_cc.csv"))
write.csv(eGO_4d_cc,  file = file.path(output, "GO_enrichment_4d_cc.csv"))

eKEGG_24h = kegg_enrich(gene_trans_24h$ENTREZID)
eKEGG_48h = kegg_enrich(gene_trans_48h$ENTREZID)
eKEGG_4d  = kegg_enrich(gene_trans_4d$ENTREZID)

eKEGG_24h_read = setReadable(eKEGG_24h, OrgDb = GO_database, keyType = "ENTREZID")
eKEGG_48h_read = setReadable(eKEGG_48h, OrgDb = GO_database, keyType = "ENTREZID")
eKEGG_4d_read = setReadable(eKEGG_4d, OrgDb = GO_database, keyType = "ENTREZID")

write.csv(eKEGG_24h, file = file.path(output, "KEGG_enrichment_24h.csv"))
write.csv(eKEGG_48h, file = file.path(output, "KEGG_enrichment_48h.csv"))
write.csv(eKEGG_4d,  file = file.path(output, "KEGG_enrichment_4d.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Gene set Enrichment analysis: GO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
go_gse = function(x) {
  clusterProfiler::gseGO(geneList = x,
                          OrgDb = GO_database,
                          keyType = "ENSEMBL",
                          minGSSize = 10,
                          maxGSSize = 500,
                          eps = 0,
                          pvalueCutoff = 0.05,
                          ont = "ALL",
                          verbose = FALSE)
}

gse_24h = go_gse(genelist_24h)
gse_48h = go_gse(genelist_48h)
gse_4d  = go_gse(genelist_4d)

df_gse_24h = as.data.frame(gse_24h)
df_gse_48h = as.data.frame(gse_48h)
df_gse_4d = as.data.frame(gse_4d)

write.csv(gse_24h, file = file.path(output, "gseGO_24h.csv"))
write.csv(gse_48h, file = file.path(output, "gseGO_48h.csv"))
write.csv(gse_4d,  file = file.path(output, "gseGO_4d.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Visualization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------------------------------------------------------------------------#
# 1. GO Bar plot
#------------------------------------------------------------------------------#
display_number = c(15,15,15)
eGO_24h_bp_part = as.data.frame(eGO_24h_bp)[1:display_number[1], ]
eGO_24h_mf_part = as.data.frame(eGO_24h_mf)[1:display_number[2], ]
eGO_24h_cc_part = as.data.frame(eGO_24h_cc)[1:display_number[3], ]

eGO_48h_bp_part = as.data.frame(eGO_48h_bp)[1:display_number[1], ]
eGO_48h_mf_part = as.data.frame(eGO_48h_mf)[1:display_number[2], ]
eGO_48h_cc_part = as.data.frame(eGO_48h_cc)[1:display_number[3], ]

eGO_4d_bp_part = as.data.frame(eGO_4d_bp)[1:display_number[1], ]
eGO_4d_mf_part = as.data.frame(eGO_4d_mf)[1:display_number[2], ]
eGO_4d_cc_part = as.data.frame(eGO_4d_cc)[1:display_number[3], ]

cols = c("#B31B21", "#1465AC", "#FBB731")

go_sum_24h = rbind(eGO_24h_bp_part, eGO_24h_cc_part, eGO_24h_mf_part)
go_sum_24h$type = c(rep("Biological Process", display_number[1]),
                    rep("Cellular Component", display_number[2]),
                    rep("Molecular Function", display_number[3]))
go_sum_24h$type = factor(go_sum_24h$type, levels = c("Biological Process","Cellular Component","Molecular Function"))
go_sum_24h$type_order = factor(c(1:sum(display_number)), labels = go_sum_24h$Description)

gobar_24h = ggplot(data = go_sum_24h, aes(x = type_order, y = Count, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cols) +
  xlab("GO term") +
  ylab("Num of Genes") +
  labs(title = "The Most Enriched GO Terms (24H, VP/DMSO)") +
  theme(axis.text = element_text(face = "bold")) +
  coord_flip() + 
  theme_bw()

go_sum_48h = rbind(eGO_48h_bp_part, eGO_48h_cc_part, eGO_48h_mf_part)
go_sum_48h$type = c(rep("Biological Process", display_number[1]),
                    rep("Cellular Component", display_number[2]),
                    rep("Molecular Function", display_number[3]))
go_sum_48h$type = factor(go_sum_48h$type, levels = c("Biological Process","Cellular Component","Molecular Function"))
go_sum_48h$type_order = factor(c(1:sum(display_number)), labels = go_sum_48h$Description)

gobar_48h = ggplot(data = go_sum_48h, aes(x = type_order, y = Count, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cols) +
  xlab("GO term") +
  ylab("Num of Genes") +
  labs(title = "The Most Enriched GO Terms (48H, VP/DMSO)") +
  theme(axis.text = element_text(face = "bold")) +
  coord_flip() + 
  theme_bw()

go_sum_4d = rbind(eGO_4d_bp_part, eGO_4d_cc_part, eGO_4d_mf_part)
go_sum_4d$type = c(rep("Biological Process", display_number[1]),
                    rep("Cellular Component", display_number[2]),
                    rep("Molecular Function", display_number[3]))
go_sum_4d$type = factor(go_sum_4d$type, levels = c("Biological Process","Cellular Component","Molecular Function"))
go_sum_4d$type_order = factor(c(1:sum(display_number)), labels = go_sum_4d$Description)

gobar_4d = ggplot(data = go_sum_4d, aes(x = type_order, y = Count, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cols) +
  xlab("GO term") +
  ylab("Num of Genes") +
  labs(title = "The Most Enriched GO Terms (4D, VP/DMSO)") +
  theme(axis.text = element_text(face = "bold")) +
  coord_flip() + 
  theme_bw()

png(filename = file.path(output, "GOterm.png"),
    width = 2100, height = 1000, units = "px")
ggarrange(gobar_24h, gobar_48h, gobar_4d,
          labels = c("A", "B", "C"),
          ncol = 3)
dev.off()

#------------------------------------------------------------------------------#
# 2. KEGG Bar Plot
#------------------------------------------------------------------------------#
kk_24h_part = as.data.frame(eKEGG_24h)[1:display_number[1], ]
kk_24h_part = kk_24h_part %>%
  dplyr::mutate(log_order = factor(kk_24h_part$Description, levels = rev(kk_24h_part$Description)))

kk_48h_part = as.data.frame(eKEGG_48h)[1:display_number[2], ]
kk_48h_part = kk_48h_part %>%
  dplyr::mutate(log_order = factor(kk_24h_part$Description, levels = rev(kk_24h_part$Description)))

kk_4d_part  = as.data.frame(eKEGG_4d)[1:display_number[3], ]
kk_4d_part = kk_4d_part %>%
  dplyr::mutate(log_order = factor(kk_24h_part$Description, levels = rev(kk_24h_part$Description)))

kkbar_24h = ggplot(data = kk_24h_part, aes(x = log_order, y = -log10(pvalue)) ) +
  geom_bar(stat = "identity", width = 0.8, fill = "#FBB731") +
  xlab( "KEGG") +
  ylab("-Log10 (pvalue)") +
  theme_bw() +
  labs(title = "The Most Enriched KEGG (24H, VP/DMSO)") +
  theme(axis.text = element_text(face = "bold")) +
  coord_flip() 

kkbar_48h = ggplot(data = kk_48h_part, aes(x = log_order, y = -log10(pvalue))) +
  geom_bar(stat = "identity", width = 0.8, fill =  "#1465AC") +
  xlab( "KEGG") +
  ylab("-Log10 (pvalue)") +
  theme_bw()+ 
  labs(title = "The Most Enriched KEGG (48H, VP/DMSO)") +
  theme(axis.text = element_text(face = "bold")) +
  coord_flip() 

kkbar_4d = ggplot(data = kk_4d_part, aes(x = log_order, y = -log10(pvalue))) +
  geom_bar(stat = "identity", width = 0.8, fill = "#B31B21") +
  xlab( "KEGG") +
  ylab("-Log10 (pvalue)") +
  theme_bw() +
  labs(title = "The Most Enriched KEGG (4D, VP/DMSO)") +
  theme(axis.text = element_text(face = "bold")) +
  coord_flip() 
png(filename = file.path(output, "KEGG.png"),
    width = 1500, height = 600, units = "px")
ggarrange(kkbar_24h,kkbar_48h, kkbar_4d,
          labels = c("A", "B", "C"),
          ncol = 3)
dev.off()

#------------------------------------------------------------------------------#
# 3. Tree Plot
#------------------------------------------------------------------------------#
tree = function(x) {
  pairwise_termsim(x) %>%
  treeplot(hclust_method = "average", nCluster = 8)
}
egoBP_24h_tree = tree(eGO_24h_bp)
egoCC_24h_tree = tree(eGO_24h_cc)
egoMF_24h_tree = tree(eGO_24h_mf)

egoBP_48h_tree = tree(eGO_48h_bp)
egoCC_48h_tree = tree(eGO_48h_cc)
egoMF_48h_tree = tree(eGO_48h_mf)

egoBP_4d_tree = tree(eGO_4d_bp)
egoCC_4d_tree = tree(eGO_4d_cc)
egoMF_4d_tree = tree(eGO_4d_mf)

png(filename = file.path(output, "egoBP_tree.png"),
    width = 3000, height = 600, units = "px")
ggarrange(egoBP_24h_tree, egoBP_48h_tree, egoBP_4d_tree,
          labels = c("A", "B", "C"),
          ncol = 3)
dev.off()

png(filename = file.path(output, "egoCC_tree.png"),
    width = 3000, height = 600, units = "px")
ggarrange(egoCC_24h_tree, egoCC_48h_tree, egoCC_4d_tree,
          labels = c("A", "B", "C"),
          ncol = 3)
dev.off()

png(filename = file.path(output, "egoMF_tree.png"),
    width = 3000, height = 600, units = "px")
ggarrange(egoMF_24h_tree, egoMF_48h_tree, egoMF_4d_tree,
          labels = c("A", "B", "C"),
          ncol = 3)
dev.off()


#------------------------------------------------------------------------------#
# 4. Ridgeline Plot
#------------------------------------------------------------------------------#
gse_24h_ridge = ridgeplot(gse_24h)
gse_48h_ridge = ridgeplot(gse_48h)
gse_4d_ridge = ridgeplot(gse_4d)
png(filename = file.path(output, "gseRidge.png"),
    width = 2400, height = 1200, units = "px")
ggarrange(gse_24h_ridge, gse_48h_ridge, gse_4d_ridge,
          labels = c("A", "B", "C"),
          ncol = 3)
dev.off()

#------------------------------------------------------------------------------#
# 5. Gene enriched by GO analysis
#------------------------------------------------------------------------------#
##### Prepare data
gene_name = t2g %>%
  dplyr::select(c("geneID", "gene_name")) %>%
  dplyr::distinct()

GO_bpId_24h = c("GO:0007178", "GO:0036503", "GO:0006506", "GO:0060389", 
                "GO:0035567", "GO:0198738", "GO:0043405", "GO:0070371",
                "GO:0070372")
GO_bpId_48h = c("GO:0007178", "GO:2001233", "GO:0016055", "GO:0198738", 
                "GO:0071900", "GO:1990778", "GO:0036503", "GO:0006506",
                "GO:0043122", "GO:0034612", "GO:0043409")
GO_bpId_4d = c("GO:0071900", "GO:0016055", "GO:0198738", "GO:0072331",
               "GO:0007178", "GO:0006506", "GO:0006505", "GO:0009615",
               "GO:0051607", "GO:0050688")
GO_bpId_intersect = intersect(intersect(GO_bpId_24h, GO_bpId_48h), GO_bpId_4d)


##### Def function
# extract genes in GO:BP of interest
GOGene = function(eGobp, goid) {
  df_GO = as.data.frame(eGobp)
  df_GO_int = df_GO[goid, ]
  geneGO = str_split(df_GO_int$geneID, "/")
  names(geneGO) = goid
  return(geneGO)
}

# find differential expression of genes of interest
### all GO BP
summary_gene = function(goid, dfResSig, geneGO, number) {
  summary = list()
  for (i in goid) {
    interest = dfResSig %>%
      dplyr::mutate(geneID = rownames(dfResSig)) %>%
      left_join(gene_name, by = "geneID") %>%
      dplyr::filter(gene_name %in% geneGO[[i]]) %>%
      dplyr::slice_max(order_by = abs(log2FoldChange), n = number)
    summary[[i]] = interest
  }
  return(summary)
}

# find differential expression of genes of interest
### shared GO BP
summary_inter_gene = function(goid, dfResSig, geneGO, time) {
  summary_inter = list()
  for (i in goid) {
    interest = dfResSig %>%
      dplyr::mutate(geneID = rownames(dfResSig)) %>%
      left_join(gene_name, by = "geneID") %>%
      dplyr::filter(gene_name %in% geneGO[[i]])
    interest$timepoint = factor(time, levels = time)
    interest = interest[order(interest$log2FoldChange), ]
    interest = dplyr::mutate(interest, 
                             gene_order = factor(interest$gene_name, 
                                                 levels = interest$gene_name))
    summary_inter[[i]] = interest
  }
  return(summary_inter)
}

summary_barplot = function(goid, summary_gene) {
  summary_plot = list()
  for (i in goid) {
    bar = ggplot(data = summary_gene[[i]],
                 aes(x = gene_name, y = log2FoldChange)) +
      geom_bar(stat='identity', width = 0.6, fill = "#1465AC") +
      labs(title = i,
           x = "Gene",
           y = "Log2 Fold Change") +
      theme_classic() +
      theme(title = element_text(face = "bold"),
            axis.title = element_text(face = "bold", size = 16),
            axis.text = element_text(size = 10))
    summary_plot[[i]] = bar
  }
  return(summary_plot)
}



GOGene_24h = GOGene(eGO_24h_bp, GO_bpId_24h)
GOGene_48h = GOGene(eGO_48h_bp, GO_bpId_48h)
GOGene_4d  = GOGene(eGO_4d_bp, GO_bpId_4d)

geneGO_inter = list()
for (i in GO_bpId_intersect) {
  gene_24h = GOGene_24h[[i]]
  gene_48h = GOGene_48h[[i]]
  gene_4d  = GOGene_4d[[i]]
  geneGO_inter[[i]] = intersect(intersect(gene_24h, gene_48h), gene_4d)
}

summaryGene_24h = summary_gene(GO_bpId_24h, df_res_24h_sig, GOGene_24h, 10)
summaryGene_48h = summary_gene(GO_bpId_48h, df_res_48h_sig, GOGene_48h, 10)
summaryGene_4d = summary_gene(GO_bpId_4d, df_res_4d_sig, GOGene_4d, 10)

summary_barplot_24h = summary_barplot(GO_bpId_24h, summaryGene_24h)
summary_barplot_48h = summary_barplot(GO_bpId_48h, summaryGene_48h)
summary_barplot_4d = summary_barplot(GO_bpId_4d, summaryGene_4d)

png(file.path(output, "GOgene_24h.png") ,
    width = 1500, height = 1200 , units = "px") 
ggarrange(summary_barplot_24h[[1]], summary_barplot_24h[[2]], summary_barplot_24h[[3]], 
          summary_barplot_24h[[4]], summary_barplot_24h[[5]], summary_barplot_24h[[6]], 
          summary_barplot_24h[[7]], summary_barplot_24h[[8]], summary_barplot_24h[[9]], 
          labels = LETTERS[1:9],
          nrow = 3, ncol = 3)
dev.off()  

png(file.path(output, "GOgene_48h.png") ,
    width = 1500, height = 1600 , units = "px") 
ggarrange(summary_barplot_48h[[1]], summary_barplot_48h[[2]], summary_barplot_48h[[3]], 
          summary_barplot_48h[[4]], summary_barplot_48h[[5]], summary_barplot_48h[[6]], 
          summary_barplot_48h[[7]], summary_barplot_48h[[8]], summary_barplot_48h[[9]], 
          summary_barplot_48h[[10]], summary_barplot_48h[[11]], 
          labels = LETTERS[1:1],
          nrow = 4, ncol = 3)
dev.off()  

png(file.path(output, "GOgene_4d.png") ,
    width = 1500, height = 1600 , units = "px") 
ggarrange(summary_barplot_4d[[1]], summary_barplot_4d[[2]], summary_barplot_4d[[3]], 
          summary_barplot_4d[[4]], summary_barplot_4d[[5]], summary_barplot_4d[[6]], 
          summary_barplot_4d[[7]], summary_barplot_4d[[8]], summary_barplot_4d[[9]], 
          summary_barplot_4d[[10]], 
          labels = LETTERS[1:10],
          nrow = 4, ncol = 3)
dev.off()  



gene_inter_24h = summary_inter_gene(GO_bpId_intersect, df_res_24h_sig, geneGO_inter, "24h")
gene_inter_48h = summary_inter_gene(GO_bpId_intersect, df_res_48h_sig, geneGO_inter, "48h")
gene_inter_4d = summary_inter_gene(GO_bpId_intersect, df_res_4d_sig, geneGO_inter, "4d")

summary_lfc_time = list()
for(i in GO_bpId_intersect) {
  time_sum = rbind(gene_inter_24h[[i]], gene_inter_48h[[i]], gene_inter_4d[[i]])
  summary_lfc_time[[i]] = time_sum
}

cols = c("#FBB731","#1465AC",  "#B31B21")
GO1_bar = ggplot(data = summary_lfc_time$`GO:0007178`,
                 aes(x = factor(gene_order), y = log2FoldChange, fill = timepoint)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge(0.7)) +
  scale_fill_manual(values = cols) +
  xlab("Gene") +
  ylab("Log2 Fold Change") +
  labs(title = "GO:0007178") +
  theme(axis.text = element_text(face = "bold")) +
  theme_bw()

GO2_bar = ggplot(data = summary_lfc_time$`GO:0006506`,
                 aes(x = factor(gene_order), y = log2FoldChange, fill = timepoint)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge(0.7)) +
  scale_fill_manual(values = cols) +
  xlab("Gene") +
  ylab("Log2 Fold Change") +
  labs(title = "GO:0006506") +
  theme(axis.text = element_text(face = "bold")) +
  theme_bw()

GO3_bar = ggplot(data = summary_lfc_time$`GO:0198738`,
                 aes(x = factor(gene_order), y = log2FoldChange, fill = timepoint)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge(0.7)) +
  scale_fill_manual(values = cols) +
  xlab("Gene") +
  ylab("Log2 Fold Change") +
  labs(title = "GO:0198738") +
  theme(axis.text = element_text(face = "bold")) +
  theme_bw()


png(file.path(output, "GO_0007178.png") ,
    width = 2500, height = 1000 , units = "px") 
  GO1_bar
dev.off() 

png(file.path(output, "GO_0006506.png") ,
    width = 1500, height = 1000 , units = "px") 
  GO2_bar
dev.off()  

png(file.path(output, "GO_0198738.png") ,
    width = 2500, height = 1000 , units = "px") 
  GO3_bar
dev.off()  
# ============================================================================ #



# ============================================================================ #
# STRINGdb: protein-protein interaction
#------------------------------------------------------------------------------#
# 1. Prepare data
string_db = STRINGdb$new(version = "11", species = 9606,
                         score_threshold = 200, input_directory = "")
##### GO id: GO_bpId_intersect
##### gene list: gene_inter_24h

# 2. Def function
ppi= function(goid, genelist) {
  list_hits = list()
  for(i in goid) {
    sample_map = string_db$map(genelist[[i]], "gene_name", removeUnmappedRows = TRUE)
    hits = sample_map$STRING_id
    list_hits[[i]] = hits
  }
  return(list_hits)
}
list_hits = ppi(GO_bpId_intersect, gene_inter_24h)

png(filename = file.path(output, "ppi_GO_0007178.png"),
    width = 1000, height = 1200, units = "px")
string_db$plot_network(list_hits$`GO:0007178`)
dev.off()

png(filename = file.path(output, "ppi_GO_0006506.png"),
    width = 1000, height =1200, units = "px")
string_db$plot_network(list_hits$`GO:0006506`)
dev.off()

png(filename = file.path(output, "ppi_GO_0198738.png"),
    width = 1000, height = 1200, units = "px")
string_db$plot_network(list_hits$`GO:0198738`)
dev.off()

png(filename = file.path(output, "ppi_summary.png"),
    width = 1000, height = 1200, units = "px")
string_db$plot_network(c(list_hits$`GO:0198738`, list_hits$`GO:0007178`))
dev.off()
# ============================================================================ #
