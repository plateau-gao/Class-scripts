#The pathway where I store the CEL.files, be careful about the last "/"
pathway = "/Users/thelittleplateau/Downloads/GSE56481_RAW/"
#The pathway where I store the processed csv.files
pathway_csv = "/Users/thelittleplateau/Desktop/data/"
# The pathway where I store the plot
pathway_plot = "/Users/thelittleplateau/Desktop/plot/"

# pathway = "C:\\Users\\lihes\\Desktop\\GSE56481_RAW\\"
# pathway_csv = "C:\\Users\\lihes\\Desktop\\csv\\"
# pathway_plot = "C:\\Users\\lihes\\Desktop\\plot\\"

# packages we need to use, if you do not have some packages, please install them first.
library(oligo)
library(GEOquery) 
library(limma) 
library(tidyverse) 
library(dplyr) 
library(readr)
library(AnnotationDbi) 
library(hugene20sttranscriptcluster.db) ##Bioconductor
library(clusterProfiler) #A universal enrichment tool for interpreting omics data ##Bioconductor
library(org.Hs.eg.db) #Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.

#get phenoData from GEO
gse = getGEO("GSE56481")
gse = gse[[1]]
gse_pdata = pData(gse)

#get raw data with annotated phenoData and find variables in the experiment.
raw_gse = oligo::read.celfiles(
  paste0(pathway, strsplit(gse_pdata$supplementary_file, "/") %>% map_chr(tail,1)),
  phenoData = phenoData(gse),
  pkgname = 'pd.hugene.2.0.st'
)
varLabels(raw_gse)
raw_gse_pdata = pData(raw_gse)
check_pdata = function(x) {unique(raw_gse_pdata[[x]])}
sapply(37:39, check_pdata) # find two variables: "diagnosis:ch1" and  "facs:ch1" 

raw_gse_rma = rma(raw_gse)

# creating a variable that represents the six groups.
raw_gse_pdata = dplyr::rename(raw_gse_pdata, diagonosis = "diagnosis:ch1", facs = "facs:ch1")
raw_gse_pdata[,c("diagonosis","facs")]
raw_gse_pdata$diagonosis = as.factor(raw_gse_pdata$diagonosis)
levels(raw_gse_pdata$diagonosis) = c("GPA", "Control")
raw_gse_pdata$group = as.factor(paste(raw_gse_pdata$diagonosis,raw_gse_pdata$facs))
levels(raw_gse_pdata$group) = c("Con.CD4","Con.CD4CD8","Con.CD8","GPA.CD4","GPA.CD4CD8","GPA.CD8")
raw_gse_pdata$group

# Identifying differentially expressed genes using linear models
design = model.matrix(~0+raw_gse_pdata$group)
colnames(design) = levels(raw_gse_pdata$group) #simplify the column name of the data.frame
contrast_martix = makeContrasts(
  gene_in_CD4 = GPA.CD4 - Con.CD4,
  gene_in_CD4CD8 = GPA.CD4CD8 - Con.CD4CD8,
  gene_in_CD8 = GPA.CD8 - Con.CD8,
  levels = design
)
raw_gse_fit = lmFit(raw_gse_rma, design)
raw_gse_fit2 = contrasts.fit(raw_gse_fit, contrasts = contrast_martix)
raw_gse_ebayes = eBayes(raw_gse_fit2)
summary(decideTests(raw_gse_ebayes, lcf = 1))

# visualization - volcanoplot
png(
  filename = paste0(pathway_plot, "CD4_volcano.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
volcanoplot(raw_gse_ebayes, coef = 1, main = "Gene alternation in  CD4+ T cell of GPA vs Control")
dev.off()

png(
  filename = paste0(pathway_plot, "CD4CD8_volcano.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
volcanoplot(raw_gse_ebayes, coef = 2, main = "Gene alternation in CD4+CD8+ T cell of GPA vs Control")
dev.off()

png(
  filename = paste0(pathway_plot, "CD8_volcano.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
volcanoplot(raw_gse_ebayes, coef = 3, main = "Gene alternation in  CD8+ T cell of GPA vs Control")
dev.off()


# find different expression of genes in CD4+ cells 
sum_CD4 = topTable(raw_gse_ebayes, coef = 1, number = Inf, p.value = 0.05, lfc = 1)
up_CD4 = sum_CD4[sum_CD4$logFC>0, ]
down_CD4 = sum_CD4[sum_CD4$logFC<0, ]
id_up_CD4 = rownames(up_CD4)
id_down_CD4 = rownames(down_CD4)


sum_CD4CD8 = topTable(raw_gse_ebayes, coef = 2, number = Inf, p.value = 0.05, lfc = 1)
up_CD4CD8 = sum_CD4CD8[sum_CD4CD8$logFC>0, ]
down_CD4CD8 = sum_CD4CD8[sum_CD4CD8$logFC<0, ]
id_up_CD4CD8 = rownames(up_CD4CD8)
id_down_CD4CD8 = rownames(down_CD4CD8)

sum_CD8 = topTable(raw_gse_ebayes, coef = 3, number = Inf, p.value = 0.05, lfc = 1)
up_CD8 = sum_CD8[sum_CD8$logFC>0, ]
down_CD8 = sum_CD8[sum_CD8$logFC<0, ]
id_up_CD8 = rownames(up_CD8)
id_down_CD8 = rownames(down_CD8)

#heatmap
raw_gse_pdata["group"] #show the pData that we need
#to show gene alternation in CD4+ T cells.
CD4 = exprs(raw_gse_rma)
CD4 = CD4[rownames(sum_CD4), c("GSM1362228", "GSM1362231", "GSM1362234", "GSM1362237", "GSM1362240", "GSM1362243")]
png(
  filename = paste0(pathway_plot, "CD4.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
heatmap(
  CD4,
  labCol = c("GPA1", "GPA2", "GPA3", "Ctl1", "Ctl2", "Ctl3"), labRow = NA,
  distfun = function(x) as.dist(1-cor(t(x))),
  main = "CD4+ T cells"
)
dev.off()

#to show gene alternation in CD8+ T cells.
CD8 = exprs(raw_gse_rma)
CD8 = CD8[rownames(sum_CD8), c("GSM1362229", "GSM1362232", "GSM1362235", "GSM1362238", "GSM1362241", "GSM1362244")]
png(
  filename = paste0(pathway_plot, "CD8.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
heatmap(
  CD8,
  labCol = c("GPA1", "GPA2", "GPA3", "Ctl1", "Ctl2", "Ctl3"), labRow = NA,
  distfun = function(x) as.dist(1-cor(t(x))),
  main = "CD8+ T cells"
)
dev.off()

#to show gene alternation in CD4+CD8+ T cells.
CD4CD8 = exprs(raw_gse_rma)
CD4CD8 = CD4CD8[rownames(sum_CD4CD8), c("GSM1362227", "GSM1362230", "GSM1362233", "GSM1362236", "GSM1362239", "GSM1362242")]
png(
  filename = paste0(pathway_plot, "CD4CD8.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
heatmap(
  CD4CD8,
  labCol = c("GPA1", "GPA2", "GPA3", "Ctl1", "Ctl2", "Ctl3"), labRow = NA,
  distfun = function(x) as.dist(1-cor(t(x))),
  main = "CD4CD8+ T cells"
)
dev.off()


# find keytype we need to use 
keytype = keytypes(hugene20sttranscriptcluster.db)
show_key = function(x) {head(keys(hugene20sttranscriptcluster.db, keytype = x))}
sapply(keytype, show_key)
#compared with data of id_up_CD4, "PROBEID" is the keytype we need.

# find annotation information
df_up_CD4 = AnnotationDbi::select(hugene20sttranscriptcluster.db, id_up_CD4, c( "SYMBOL", "ENTREZID", "GENENAME"), keytype = "PROBEID")
df_up_CD4CD8 = AnnotationDbi::select(hugene20sttranscriptcluster.db, id_up_CD4CD8, c("SYMBOL", "ENTREZID", "GENENAME"), keytype = "PROBEID")
df_up_CD8 = AnnotationDbi::select(hugene20sttranscriptcluster.db, id_up_CD8, c("SYMBOL", "ENTREZID", "GENENAME"), keytype = "PROBEID")
df_down_CD4 = AnnotationDbi::select(hugene20sttranscriptcluster.db, id_down_CD4, c("SYMBOL", "ENTREZID", "GENENAME"), keytype = "PROBEID")
df_down_CD4CD8 = AnnotationDbi::select(hugene20sttranscriptcluster.db, id_down_CD4CD8, c("SYMBOL", "ENTREZID", "GENENAME"), keytype = "PROBEID")
df_down_CD8 = AnnotationDbi::select(hugene20sttranscriptcluster.db, id_down_CD8, c("SYMBOL", "ENTREZID", "GENENAME"), keytype = "PROBEID")

# we need to delete probes which cannot be annotated.
df_up_CD4 = drop_na(df_up_CD4, "ENTREZID") 
df_up_CD4CD8 = drop_na(df_up_CD4CD8, "ENTREZID")
df_up_CD8 = drop_na(df_up_CD8, "ENTREZID")
df_down_CD4 = drop_na(df_down_CD4, "ENTREZID")
df_down_CD4CD8 = drop_na(df_down_CD4CD8, "ENTREZID")
df_down_CD4 = drop_na(df_down_CD8, "ENTREZID")

#GO enrichment analysis
go_enrich = function(x) {
  clusterProfiler::setReadable(
    clusterProfiler::simplify(
      enrichGO(x, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")
    ),
    OrgDb = org.Hs.eg.db
  )
}
CD4_down_GO = go_enrich(df_down_CD4$ENTREZID)
CD4_up_GO = go_enrich(df_up_CD4$ENTREZID)
CD4CD8_down_GO = go_enrich(df_down_CD4CD8$ENTREZID)
CD4CD8_up_GO = go_enrich(df_up_CD4CD8$ENTREZID)
CD8_down_GO = go_enrich(df_down_CD8$ENTREZID)
CD8_up_GO = go_enrich(df_up_CD8$ENTREZID)

#Visualize the GO data
#CD4+down
png(
  filename = paste0(pathway_plot, "CD4_down_GO.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
  )
barplot(CD4_down_GO, title = "Biological progress down regulated in CD4+ T cells of GPA", label_format = 30)
dev.off()

#CD4+up
png(
  filename = paste0(pathway_plot, "CD4_up_GO.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD4_up_GO, title = "Biological progress up regulated in CD4+ T cells of GPA", label_format = 30)
dev.off()

#CD8+down
png(
  filename = paste0(pathway_plot, "CD8_down_GO.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD8_down_GO, title = "Biological progress down regulated in CD8+ T cells of GPA", label_format = 30)
dev.off()

#CD8+up
png(
  filename = paste0(pathway_plot, "CD8_up_GO.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD8_up_GO, title = "Biological progress up regulated in CD8+ T cells of GPA", label_format = 30)
dev.off()

#CD4CD8_down
png(
  filename = paste0(pathway_plot, "CD4CD8_down_GO.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD4CD8_down_GO, title = "Biological progress down regulated in CD4+CD8+ T cells of GPA", label_format = 20)
dev.off()

#CD4CD8_up
png(
  filename = paste0(pathway_plot, "CD4CD8_up_GO.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD4CD8_up_GO, title = "Biological progress up regulated in CD4+CD8+ T cells of GPA", label_format = 30)
dev.off()

#KEGG enrichment

CD4_down_kegg = enrichKEGG(df_down_CD4$ENTREZID, organism = "hsa")
CD4_up_kegg = enrichKEGG(df_up_CD4$ENTREZID, organism = "hsa")
CD4CD8_down_kegg = enrichKEGG(df_down_CD4CD8$ENTREZID, organism = "hsa")
CD4CD8_up_kegg = enrichKEGG(df_up_CD4CD8$ENTREZID, organism = "hsa")
CD8_down_kegg = enrichKEGG(df_down_CD8$ENTREZID, organism = "hsa")
CD8_up_kegg = enrichKEGG(df_up_CD8$ENTREZID, organism = "hsa")

#Visualize the KEEG data
#CD4+down
png(
  filename = paste0(pathway_plot, "CD4_down_kegg.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD4_down_kegg, title = "Pathway down regulated in CD4+ T cells of GPA", label_format = 30)
dev.off()

#CD4+up
png(
  filename = paste0(pathway_plot, "CD4_up_kegg.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD4_up_kegg, title = "Pathway up regulated in CD4+ T cells of GPA", label_format = 30)
dev.off()

#CD8+down
png(
  filename = paste0(pathway_plot, "CD8_down_kegg.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD8_down_kegg, title = "Pathway down regulated in CD8+ T cells of GPA", label_format = 30)
dev.off()

#CD8+up
png(
  filename = paste0(pathway_plot, "CD8_up_kegg.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD8_up_kegg, title = "Pathway up regulated in CD8+ T cells of GPA", label_format = 30)
dev.off()

#CD4CD8_down
png(
  filename = paste0(pathway_plot, "CD4CD8_down_kegg.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD4CD8_down_kegg, title = "Pathway down regulated in CD4+CD8+ T cells of GPA", label_format = 20)
dev.off()

#CD4CD8_up
png(
  filename = paste0(pathway_plot, "CD4CD8_up_kegg.png"),
  width = 1000, height = 750, units = "px", pointsize = 20, res = 100
)
barplot(CD4CD8_up_kegg, title = "Pathway up regulated in CD4+CD8+ T cells of GPA", label_format = 30)
dev.off()


#save data CSV file
df_CD4_down_GO = data.frame(CD4_down_GO)
df_CD4_up_GO = data.frame(CD4_up_GO)
df_CD8_down_GO = data.frame(CD8_down_GO)
df_CD8_up_GO = data.frame(CD8_up_GO)
df_CD4CD8_down_GO = data.frame(CD4CD8_down_GO)
df_CD4CD8_up_GO = data.frame(CD4CD8_up_GO)

df_CD4_down_kegg = data.frame(CD4_down_kegg)
df_CD4_up_kegg = data.frame(CD4_up_kegg)
df_CD8_down_kegg = data.frame(CD8_down_kegg)
df_CD8_up_kegg = data.frame(CD8_up_kegg)
df_CD4CD8_down_kegg = data.frame(CD4CD8_down_kegg)
df_CD4CD8_up_kegg = data.frame(CD4CD8_up_kegg)

write_excel_csv(df_CD4_down_GO, file = paste0(pathway_csv, "CD4_down_GO.csv"))
write_excel_csv(df_CD4_up_GO, file = paste0(pathway_csv, "CD4_up_GO.csv"))
write_excel_csv(df_CD8_down_GO, file = paste0(pathway_csv, "CD8_down_GO.csv"))
write_excel_csv(df_CD8_up_GO, file = paste0(pathway_csv, "CD8_up_GO.csv"))
write_excel_csv(df_CD4CD8_down_GO, file = paste0(pathway_csv, "CD4CD8_down_GO.csv"))
write_excel_csv(df_CD4CD8_up_GO, file = paste0(pathway_csv, "CD4CD8_up_GO.csv"))

write_excel_csv(df_CD4_down_kegg, file = paste0(pathway_csv, "CD4_down_kegg.csv"))
write_excel_csv(df_CD4_up_kegg, file = paste0(pathway_csv, "CD4_up_kegg.csv"))
write_excel_csv(df_CD8_down_kegg, file = paste0(pathway_csv, "CD8_down_kegg.csv"))
write_excel_csv(df_CD8_up_kegg, file = paste0(pathway_csv, "CD8_up_kegg.csv"))
write_excel_csv(df_CD4CD8_down_kegg, file = paste0(pathway_csv, "CD4CD8_down_kegg.csv"))
write_excel_csv(df_CD4CD8_up_kegg, file = paste0(pathway_csv, "CD4CD8_up_kegg.csv"))

