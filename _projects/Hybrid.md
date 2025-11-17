---
layout: page
title: Hybrid Genome Assembly Pipeline
description: combining ONT Long-Read Sequencing and Bionano Optical Genome Mapping for highly contiguous whole genome assemblies.
img: assets/img/dna_background.jpg
importance: 2
category: Bioinformatics
giscus_comments: true
---

#### **About**

The MCW Advanced Genomics lab specializes in two sequencing platforms: Oxford Nanopore Long-Read Sequencing (LRS) and Bionano Optical Genome Mapping (OGM). These are two distinct technologies that accomplish two separate goals. OGM uses molecular barcodes to detect medium-to-large structural variants using long strands of DNA that are tagged and stretched out so their patterns can be read. LRS, on the other hand, sequences base-pairs for precise readings of DNA at the most specific of resolution. **What if we could harness the strengths of both technologies?**

In 2023, Bionano published a [Hybrid Scaffold](https://bionano.com/wp-content/uploads/2023/08/CG-30073-Bionano-Solve-Theory-of-Operation-Hybrid-Scaffold.pdf) workflow onto their integrated analysis platform. This tool enabled users to combine sequencing data with OGM data to build scaffolded genomic [assemblies](https://useast.ensembl.org/info/genome/genebuild/assembly.html). Although it was designed for short-read data and later deprecated as the company shifted focus, our lab saw an opportunity to extend its utility. 

## Purpose

Genomic assemblies are reconstructed versions of an organism’s genome, built by stitching sequencing reads into longer sequences so we can view DNA as continuous contigs, scaffolds, and ideally chromosome-scale structures instead of millions of isolated reads. LRS provides base-level detail and resolves repeats, while OGM contributes large-scale structural information that helps order and orient contigs into larger, more accurate scaffolds.

Assembly quality is evaluated using contiguity metrics ([N50](https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics)), completeness ([BUSCO scores](https://academic.oup.com/bioinformatics/article/31/19/3210/211866)), and structural accuracy. High-quality assemblies enable more reliable structural variant discovery, resolution of complex regions, and the creation of personalized or population-specific references that improve downstream analysis. By integrating ONT and OGM, this project aims to generate hybrid assemblies that combine precise sequence information with robust large-scale structure.

## Methods

I adapted the Hybrid Scaffold approach for use with our long-read sequencing data and optimized the workflow for our lab samples. 

1. #### **Find the Deprecated Code**

   The first step was locating Bionano’s original Hybrid Scaffold workflow within the Access/Compute software package [downloadable](https://bionano.com/software-downloads/) from their website. This codebase is large and tightly integrated, containing modules for mapping, assembly, annotation, visualization, and database management. Because Hybrid Scaffold was deprecated years ago, it is no longer documented or exposed through the graphical interface, so identifying it required a systematic search through the underlying directories.

   

2. #### **Learn and Ingest the Code**

   After finding the Hybrid Scaffold components, I had to understand how they actually worked. The workflow is implemented almost entirely in Perl, with dense subroutines and little documentation, so I systematically read through the scripts to map inputs, outputs, and file handoffs between stages. I specifically looked for where the sequence data, OGM maps, and parameters entered the pipeline, and how intermediate results were written out. Throughout this process, I leaned heavily on large language models to translate unfamiliar Perl idioms into Python-style logic, summarize complex functions, and clarify module interactions, allowing me to quickly build a workable mental model of the original workflow.

3. #### **Generate Sufficient Sequencing Data**

   The [Hybrid Scaffold user manual](https://bionano.com/wp-content/uploads/2023/01/30324-Bionano-Access-Hybrid-Scaffold-Report-Guidelines.pdf) recommends a minimum contig N50 of 150 kb for input assemblies, so the first requirement on the sequencing side was to produce data capable of meeting or exceeding this threshold. To achieve this, I assembled the ONT long-read data using the [Flye](https://github.com/mikolmogorov/Flye) assembler, which is optimized for noisy long reads and designed to maximize contig length. By running Flye on our raw ONT data and iterating on parameters as needed, I was able to boost the N50 of the resulting contigs and generate an assembly that satisfied the 150 kb N50 requirement, making it suitable for downstream hybrid scaffolding with OGM.

4. #### Generate & Evaluate the Assembly

   After preparing the data and integrating the required tools, I executed the Hybrid Scaffold workflow by stepping through each script in sequence to ensure the long-read assembly was formatted and processed correctly. Once the pipeline completed, we generated standard assembly metrics (N50, total assembly size, and structural concordance) to evaluate performance. The resulting hybrid assemblies were highly contiguous providing a reliable and representative reconstruction of each genome.

5. #### **Build a Pipeline**

   To streamline the entire process, I developed a Python-based pipeline that orchestrates every step required to prepare ONT data for hybrid assembly. This workflow automates long-read assembly, integrates a phasing step to separate maternal and paternal haplotypes, and formats all outputs into the structures expected by the legacy Hybrid Scaffold modules. Each component is executed with consistent logging and error checking. By consolidating these tasks into a single executable workflow, I created a system that reliably processes new samples and produces hybrid-ready assemblies with minimal manual intervention.

   

## Code Sample

A snippet of the Hybrid Assembly Pipeline Code.

{% highlight R linenos %}

# --- Step 5: Parallel Flye Assemblies ---
    logger.info("--- Step 5: Submitting Flye assembly jobs in parallel ---")
    
    assembly_fasta_paths = {
        "H1": os.path.join(dirs["flye_split"], "H1", "assembly.fasta"),
        "H2": os.path.join(dirs["flye_split"], "H2", "assembly.fasta"),
        "unsplit": os.path.join(dirs["flye_unsplit"], "unsplit", "assembly.fasta")
    }
    
    if all(os.path.exists(f) for f in assembly_fasta_paths.values()):
        logger.info("All Flye assembly FASTA files already exist. Skipping Flye assembly step.")
    else:
        flye_commands = [
            f"conda run -n flye flye --nano-hq {h1_fastq} --out-dir {os.path.join(dirs['flye_split'], 'H1')} --threads {flye_threads_per_job}",
            f"conda run -n flye flye --nano-hq {h2_fastq} --out-dir {os.path.join(dirs['flye_split'], 'H2')} --threads {flye_threads_per_job}",
            f"conda run -n flye flye --nano-hq {os.path.abspath(args.input_fastq)} --out-dir {os.path.join(dirs['flye_unsplit'], 'unsplit')} --threads {flye_threads_per_job}"
        ]
        with ProcessPoolExecutor(max_workers=3) as executor:
            futures = [executor.submit(run_command, cmd, logger) for cmd in flye_commands]
            for future in as_completed(futures):
                if not future.result():
                    logger.error("A Flye assembly job failed. Check logs. Exiting.")
                    sys.exit(1)
        logger.info("--- Flye assemblies completed. ---")
    
    # --- Steps 6-9 (Parallel Analysis Workflow) ---
    logger.info("--- Steps 6-9: Submitting main analysis workflows in parallel ---")
    final_hybrid_paths = []
    
    assemblies_to_process = {
        "H1": assembly_fasta_paths["H1"],
        "H2": assembly_fasta_paths["H2"],
        "unsplit": assembly_fasta_paths["unsplit"]
    }
    
    with ProcessPoolExecutor(max_workers=args.max_parallel_jobs) as executor:
        future_to_assembly = {
            executor.submit(process_single_assembly, base_name, fasta_path, args, dirs, threads_per_job, logger): base_name
            for base_name, fasta_path in assemblies_to_process.items() if os.path.exists(fasta_path)
        }
        for future in as_completed(future_to_assembly):
            assembly_name = future_to_assembly[future]
            try:
                result_path = future.result()
                if result_path:
                    logger.info(f"Successfully processed and generated hybrid scaffold for {assembly_name}")
                    final_hybrid_paths.append(result_path)
                else:
                    logger.error(f"Processing failed for {assembly_name}. See worker logs for details.")
            except Exception as exc:
                logger.error(f'{assembly_name} generated an exception: {exc}')



{% endhighlight %}



Here are some assembly statistics from our hybrid assemblies.

| File     | Seqs    | Total bp      | N50         | min      | max          | N75         | N90         | auN         |
| -------- | ------- | ------------- | ----------- | -------- | ------------ | ----------- | ----------- | ----------- |
| Raw 1    | 1837.00 | 2692825500.00 | 6452351.00  | 277.00   | 37508636.00  | 3446441.00  | 1505842.00  | 7832347.50  |
| Raw 2    | 1995.00 | 2741310110.00 | 6941736.00  | 440.00   | 32187494.00  | 3492158.00  | 1338670.00  | 8403763.53  |
| Hybrid 1 | 169.00  | 3148684920.00 | 62654798.00 | 83587.00 | 138302816.00 | 36582336.00 | 23812109.00 | 68756155.39 |
| Hybrid 2 | 163.00  | 3112696750.00 | 58418037.00 | 70082.00 | 138312544.00 | 37010218.00 | 24296380.00 | 69158919.08 |
