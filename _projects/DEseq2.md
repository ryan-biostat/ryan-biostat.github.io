---
layout: page
title: Differential Expression Analysis
description: for drug discovery using DESeq2 in R
img: assets/img/deseq2_maplot.png
importance: 1
category: Biostatistics
---

#### **About**

I’m part of a grant-funded project studying how different drugs change gene expression in the human genome. It’s been a great crash course in how RNA behaves and how we can actually measure those changes in real data.

On this project I’ve gotten hands-on experience working with expression data, collaborating with people from different disciplines, and communicating results clearly. One hallmark of this project was weekly presentations of findings to our group of experts, this was a slightly intimidating but really rewarding experience that helped me build a lot of confidence in my scientific communication and teamwork skills. The insights from my analysis have already helped shape the next steps of the research, and the project is still actively moving forward.

## Purpose

Tyrosine kinase inhibitors (TKIs) are a class of targeted cancer therapies that block enzymes called tyrosine kinases, which help control cell signaling, growth, and division. According to the [National Cancer Institute](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/tyrosine-kinase-inhibitor?utm_source=chatgpt.com), these enzymes can be overactive or dysregulated in some cancers, and inhibiting them can slow tumor growth and promote cancer cell death.



In this project, I used [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), an R/Bioconductor package for RNA-seq differential expression analysis, to quantify how gene expression changed across cell types, drugs, and doses. This method models count data with a negative binomial distribution and uses moderated estimates of dispersion and fold change to identify genes that are significantly up- or down-regulated between conditions. After generating lists of differentially expressed genes, I performed pathway and gene-set enrichment analysis using [ShinyGO](https://bioinformatics.sdstate.edu/go/) and the [KEGG](https://www.genome.jp/kegg/) database. ShinyGO, developed at South Dakota State University, provides graphical enrichment analysis and integrates with KEGG and other resources to map gene lists onto functional categories and pathways. Together, DESeq2-based differential expression and pathway analysis allowed me to move from raw RNA counts to an interpretation about how TKIs influence cellular pathways.

## Methods

Analysis was performed via the following steps:

1. #### **Data Accumulation**

   Lab members treated cells with several different TKIs at a range of doses to see how their gene expression would change. After treatment, the cells were collected and processed for RNA extraction using standard molecular biology workflows. The extracted RNA was then quantified and sequenced to generate count data representing gene expression levels under each treatment condition.

2. #### **Data Normalization**

   The RNA-seq data were structured as a gene-by-sample count matrix, with genes as rows, samples as columns, and raw read counts in each cell. Because sequencing depth and compositional differences vary across samples, normalization was necessary before making meaningful comparisons. To address this, I implemented DESeq2’s [*median-of-ratios*](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) normalization, which calculates sample-specific size factors by comparing each sample’s counts to a pseudo-reference derived from gene-wise geometric means. By applying these size factors through the DESeq2 pipeline, I ensured that the normalized counts were appropriately scaled and ready for differential expression analysis.

3. #### **DESeq2**

   To perform the differential expression (DE) analysis, I implemented a workflow in R using the **DESeq2** package. After importing the raw count matrix and corresponding metadata, I first stratified the samples into *resistant* and *sensitive* groups based on their annotated cell-line classifications. I then filtered out low-abundance genes by removing any gene whose median count across samples was ≤10. This step reduces noise and improves stability in the model estimation.

   For each tyrosine kinase inhibitor (TKI) and dose combination, I constructed a DESeq2 dataset containing only the relevant treated samples. In the unpaired analysis, I modeled gene expression using a simple design formula `~ Class`, which tests for differences between resistant and sensitive groups within each treatment condition.

   I also implemented an extended **paired-controls design** to account for baseline expression. In this version of the analysis, the treated samples were combined with their plate-matched DMSO controls. I then used the design formula `~ treatment + Class`, which adjusts for untreated-versus-treated differences before estimating the resistant-versus-sensitive contrast. This allowed the analysis to isolate differences driven specifically by the drug response rather than preexisting variation.

   DESeq2’s pipeline handled dispersion estimation and fitting of the negative binomial models. For each analysis, I extracted the ordered results table and exported it as a `.tsv` file for downstream interpretation and pathway enrichment analysis.

4. #### **Pathway Analysis**

   After generating the differential expression results for each TKI and dose combination, I performed pathway enrichment analysis to interpret the biological patterns underlying the gene-level changes. For each experiment, I extracted the list of significantly differentially expressed genes and uploaded it to **ShinyGO**, an online enrichment platform that integrates Gene Ontology terms, KEGG pathways, and other functional annotation resources.

   ShinyGO provided ranked lists of enriched pathways, along with statistical metrics such as fold enrichment and adjusted p-values. I downloaded these results for every TKI–dose pair and then processed them in R to create visual summaries. Using these outputs, I generated figures highlighting the top enriched pathways for each experimental condition. These visualizations helped compare how different TKIs and doses perturbed cellular signaling and functional processes, and they provided a higher-level interpretation of the transcriptomic responses identified in the DESeq2 analysis.

## Code Sample

The following code displays the core script for running differential expression analysis via DESeq2: 

{% highlight R linenos %}

# Ryan Gallagher

#   Goal: This script will:
##          1. Import and stratify by group
##          2. Remove low abundant genes
##          3. DESeq2 of Resistant vs. Sensitive for each trt+dosage
##          4. Output as .tsv files in the analysis folder
#
#             In this approach, we incorporate control samples to account for baseline expression levels before comparing 
#             resistant and sensitive groups. This helps remove genes that are already differentially 
#             expressed at baseline (before treatment) and focuses on treatment-induced differences 
#             between resistant and sensitive samples.


library(tidyverse)
library(DESeq2)
library(matrixStats)
library(readxl)

cts = read.table("../data/TKI_phaseI_raw_read_counts_913_samples_final.txt")
cts = as.matrix(cts)

meta = read_excel("../data/metadata/TKI Phase I samples to use for analysis.xlsx", sheet=2)
meta = meta %>% mutate(treatment = ifelse(Group == "DMSO","untreated", "treated"))
meta$treatment = factor(meta$treatment, levels=c("untreated", "treated"))

cell_lines = unique(sort(sub("_.*","", colnames(cts))))
KIs = unique(na.omit(str_extract(colnames(cts), "(?<=_)[A-Z]+(?=_)")))
dosage = c("1x", "3x", "10x")

## The cell lines:
sensitive.3 = c("c1120", "c1343", "c1367")
resistant.3 = c("c1254", "c1297", "c1115")

meta.filtered = meta %>%
  filter(str_detect(SampleName, paste0("^(", paste(c(sensitive.3, resistant.3), collapse = "|"), ")"))) %>%
  mutate(Class = ifelse(CellLine %in% resistant.3, 'resistant', 'sensitive'))
cts.filtered = cts[, colnames(cts) %in% meta.filtered$SampleName]

## Separate Controls and Treatments
meta.controls = meta.filtered %>% filter(Group == 'DMSO')
cts.controls = cts.filtered[, colnames(cts.filtered) %in% meta.controls$SampleName]

meta.treatments = meta.filtered %>% filter(Group != "DMSO")
cts.treatments = cts.filtered[, colnames(cts.filtered) %in% meta.treatments$SampleName]


RvS_DESeq2 = function(KI, dose, paired=F) {

  # Check inputs
  if(!is.character(KI) | !is.character(dose)) return(NULL)

  meta.instance = meta.treatments %>% filter(Condition == KI, Group == dose)
  cts.instance = cts.treatments[, match(meta.instance$SampleName, colnames(cts.treatments))]
  meta.instance$Class = factor(meta.instance$Class, levels=c('resistant', 'sensitive'))

  ## REMOVE LOW ABUNDANCE
  d1 = dim(cts.instance)
  gene.medians=rowMedians(as.matrix(cts.instance))

  ## get the indexes of the genes with less than 10 median read count
  index.lowCountGenes=which(gene.medians<=10)

  ## remove these low abundance genes
  cts.instance=cts.instance[-index.lowCountGenes,]
  if (paired == F) {
  d2 = dim(cts.instance)
  print(paste("Low Abdundance Filtering Removed:", d1[1]-d2[1], "Gene Rows.", sep=' '))
  }
  #Skip if no other classification group
  if (length(unique(meta.instance$Class)) < 2) {
    return("skipping")
  }

  if (nrow(meta.instance) == 6) {
    directory = "REDACTED"
  }

  if (nrow(meta.instance) != 6) {
    directory = "REDACTED"
  }

  dds = DESeqDataSetFromMatrix(countData = cts.instance,
                               colData = meta.instance,
                               design = ~ Class) 

  file.name = paste(directory, "/", KI, "_", dose, ".tsv", sep='')

  if (paired==T) {
    meta.controls.instance = meta.controls %>% filter(Plate == unique(meta.instance$Plate))
    cts.controls.instance = cts.controls[, colnames(cts.controls) %in% meta.controls.instance$SampleName]
    
    plate = unique(meta.instance$Plate)
    
    meta.instance = bind_rows(meta.instance, meta.controls.instance)
    cts.instance = cts.filtered[, match(meta.instance$SampleName, colnames(cts.filtered))]
    meta.instance$Class = factor(meta.instance$Class, levels=c('resistant', 'sensitive'))
    meta.instance$treatment = factor(meta.instance$treatment, levels=c('untreated', 'treated'))
    
    print(paste("# of Genes in Table", nrow(cts.instance)))
    
    if (length(unique(meta.instance$treatment)) < 2) {
      return("skipping")
    }


​    
​    ## REMOVE LOW ABUNDANCE
​    d1 = dim(cts.instance)
​    gene.medians=rowMedians(as.matrix(cts.instance))
​    index.lowCountGenes=which(gene.medians<=10)
​    cts.instance=cts.instance[-index.lowCountGenes,]
​    d2 = dim(cts.instance)
​    print(paste("Low Abdundance Filtering Removed:", d1[1]-d2[1], "Gene Rows.", sep=' '))
​      
​    if (nrow(meta.instance) == 12) {
​      directory = "REDACTED"
​    }
​    
    if (nrow(meta.instance) != 12) {
      directory = "REDACTED"
    }
    
    dds = DESeqDataSetFromMatrix(countData = cts.instance,
                                 colData = meta.instance,
                                 design = ~ treatment + Class)
    
    file.name = paste(directory, "/", KI, "_", dose, "_paired.tsv", sep='')
  }

  print("Beginning DESeq Functon")

  dds = DESeq(dds)
  res = results(dds)
  resOrdered = res[order(res$pvalue),]

  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }


  print(file.name)
  write.table(cbind(Genes = rownames(resOrdered), as.data.frame(resOrdered)), 
              file = file.name, 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
}

KIs = unique(meta.treatments$Condition)
dosings = unique(meta.treatments$Group)

for (KI in KIs) {
  for (doses in dosings) {
    RvS_DESeq2(KI, doses, paired=T)
  }
}

{% endhighlight %}

After pathway enrichment, multiple figures (as below) were generated for presentation.

{% include figure.liquid 
    path="assets/img/pathway_img.png"
    title="pathway Image"
    class="img-fluid rounded z-depth-1"
%}