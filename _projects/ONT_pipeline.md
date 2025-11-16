---
layout: page
title: Long-Read Sequencing Pipeline
description: customized for Whole Genomes on Oxford Nanopore Long-Read Platforms
img: assets/img/dna_sequencing.jpg
importance: 1
category: Bioinformatics
related_publications: true
---

## Intro

This pipeline is designed for processing whole genome samples from the Oxford Nanopore Technologies (ONT) PromethION long-read sequencing platform. In our lab, we process these samples on R10.9.1 flowcells, and this workflow interprets the resulting raw sequence data.

The pipeline is managed by a central Python script that feeds the data through a series of specialized bioinformatics software to move from raw data to an annotated list of genetic variants.

## Method

The pipeline consists of the following key steps:

#### 1. **Transfer**
Files are prepared for cross-server transfer. An `rsync` command is used to copy over raw sequencing data from the sequencing instrument to our lab's High-Performance Computing (HPC) cluster, named "maple".

#### 2. **Alignment**
Raw sequences are aligned to both the hg38 (GRCh38) and T2T (Telomere-to-Telomere) human reference genomes using `minimap2`. This process generates a BAM (.bam) file, which serves as the foundational file for subsequent analyses.

#### 3. **CNV Calling**
`ont-spectre`, a Copy Number Variant (CNV) calling tool developed by ONT, is utilized. Although still in development, we use it to create VCF (.vcf) files capable of identifying and annotating large structural variants, specifically those greater than 100kb. This step produces both .vcf and .bed files to aid in analysis.

#### 4. **SV Calling**
`sniffles2` is employed as our primary structural variant (SV) caller. It identifies a range of SVs, including deletions, insertions, translocations, and inversions. The typical size range for variants detected by `sniffles2` is between 10bp and 50kb. The output is a .vcf file.

#### 5. **SNP/Small INDEL Calling**
We use `deepvariant` from Google's DeepMind to call small genetic variants. This software utilizes a deep neural network to analyze the sequencing data. This step calls single nucleotide polymorphisms (SNPs)—which are the most common type of genetic variation, representing a change in a single DNA building block (nucleotide)—and small INsertions and DELetions (INDELs). This analysis also produces haplotyping information, allowing us to determine which genome strand belongs to which parent.

#### 6. **Phasing**
`whatshap` is used to assign haplotypes to our .bam files. This assignment allows us to determine which allele (parental copy) holds which mutations. This is critical for determining the inheritance mode for potential genetic diseases. `whatshap` achieves this by taking the information produced by `deepvariant` and annotating the .bam file, resulting in a .haplotagged.bam file.

#### 7. **Annotation**
`snpEFF` is used to annotate the small variants (SNPs and INDELs). It enriches the .vcf file by pulling information from multiple public databases, such as ClinVar and dbSNP, to provide aggregate data about known mutations.

We also utilize the **Geneyx** software. This is a proprietary annotation platform that adds richer annotations and formats the datasets in a more contestable, browser-accessible way.

## Results / Discussion
