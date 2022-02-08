# path = "/Users/thelittleplateau/Downloads/BigWigs/"
path = "C:\\Users\\Carl\\Desktop\\BigWigs"
setwd(path)
getwd()

library(rtracklayer)
library(tximport)
library(DESeq2)
library(dplyr)
library(Gviz)
library(AnnotationHub)
library(rhdf5)
library(GenomicFeatures)
library(shiny)
library(shinythemes)
library(tidyverse)

#DESeq2 analysis
metadata = read.csv("tdp43_gfap_design.csv", 
                    header = TRUE, 
                    stringsAsFactors = FALSE)
metadata = dplyr::mutate(metadata, 
                         path = file.path("kallisto_gencode", 
                                          sample, 
                                          "abundance.h5"))
metadata = dplyr::select(metadata, c("sample", "condition", "gender", "path"))
metadata$condition <- factor(metadata$condition, levels = c("ctr", "cHet", "cKO"))
metadata

t2g = read.csv("tx2gene_ensembl_v104.csv")
txi = tximport(metadata$path, type = "kallisto", tx2gene = t2g)

dds = DESeqDataSetFromTximport(txi, metadata, ~ gender + condition)
dds = DESeq(dds)
resultsNames(dds)
result_cKO = results(dds, name="condition_cKO_vs_ctr")
result_cKO <- result_cKO[order(result_cKO$padj),]
result_cKO$geneid = rownames(result_cKO)
result_cKO$geneid[1:10]
# [1] "ENSMUSG00000021508" "ENSMUSG00000016194" "ENSMUSG00000042251" "ENSMUSG00000027674" "ENSMUSG00000051497" "ENSMUSG00000020334" "ENSMUSG00000006522"
# [8] "ENSMUSG00000033065" "ENSMUSG00000005089" "ENSMUSG00000039323"
######################################################################

# creat gene model from GTF file.
txdb = makeTxDbFromGFF("gencode.vM27.primary_assembly.annotation.gtf", format = "gtf")
seqlevels(txdb) #check the seq-name, such as whether it is "1" or "chr1"
tx2gene <- transcriptsBy(txdb, by = "gene")
# change the gene name to match the gene name we get from DESeq2
names(tx2gene) = strsplit(names(tx2gene), "[.]") %>% map_chr(head, 1)


# import bigwig file, the seqlevel need to be changed
bwi_1641 = import.bw("F1641_SPC.bw")
bwi_1703 = import.bw("F1703_SPC.bw")

options(ucscChromosomeNames = F)


# PlotResultForGene function
PlotResultsForGene = function(x) {
  ensemble = tx2gene[[x]]
  sta = min(start(ensemble))
  en = max(end(ensemble))
  chr = unique(as.character(seqnames(ensemble)))
  print(c("start" = sta, "end" = en, "chromosome" = chr))
  print(ensemble)
  
  gtrack = GenomeAxisTrack()
  grTrack = GeneRegionTrack(txdb,
                            chromosome = chr, start  = sta,end = en,
                            name = "transcript",
                            transcriptAnnotation = "geneid")
  
  ctr = DataTrack(bwi_1641, chromosome = chr,
                  type = "l", name = "Control", col = "brown", background.title = "brown")
  cko = DataTrack(bwi_1703, chromosome = chr,
                  type = "l", name = "cKO", col = "darkblue", background.title = "darkblue")
  ctr1 = DataTrack(bwi_1641, chromosome = chr,
                  type = "l", name = "Ctr vs cKO", col = "brown", background.title = "brown")
  
  ot = OverlayTrack(trackList = list(ctr1, cko), background.title = "black")
  
  plotTracks(
    list(gtrack, grTrack, ctr, cko, ot),
    from =sta , to = en, 
    extend.left = 0.5, extend.right = 5000)
}


# Shiny package
ui = fluidPage(
  theme = shinytheme("lumen"),
  titlePanel("Top 50 Gene"),
  selectInput(
    inputId = "gene",
    label = "Choose a Gene ID",
    choices = result_cKO$geneid[1:50],
    selected = "ENSMUSG00000042251"
  ),
  plotOutput(outputId = "lineplot", height = "400px")
)
server = function(input, output) {
  output$lineplot = renderPlot(
    {
      PlotResultsForGene(as.character(input$gene))
    }
  )
}
shinyApp(ui = ui, server = server)

