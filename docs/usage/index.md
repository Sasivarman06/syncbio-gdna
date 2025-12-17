

  <p>
    SyncBio-gDNA is a specialized pipeline for rare disease variant calling, supporting both solo (proband-only) and trio (proband + parents) analysis modes with comprehensive annotation and reporting.
  </p>

</div>  

## Overview

Welcome to SyncBio-gDNA's documentation! This guide is the main source of documentation for users who are getting started with our rare disease variant calling pipeline.

The **`syncbio-gdna`** pipeline is composed of several interrelated sub-commands to set up and run the pipeline across different systems. Each of the available sub-commands performs different functions: 

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">



</section>

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">



</section>

**SyncBio-gDNA** is a specialized pipeline focused on rare disease variant discovery. It employs best-in-class tools for alignment, variant calling, annotation, and filtering, with special attention to the needs of rare disease genomics. The pipeline supports both whole genome and exome sequencing data.

The workflow is orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

As input, the pipeline accepts either:
- FASTQ files (for raw read processing)
- BAM files (for pre-aligned data)
- A samplesheet for batch processing


## Quick Start

```bash
./syncbio-gdna run --fq1 sample_R1.fastq.gz --fq2 sample_R2.fastq.gz --output results/
```
## References

<sup>**1.** Van der Auwera GA, O’Connor BD (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra. O’Reilly Media.</sup>

<sup>**2.** Poplin R, et al. (2018). A universal SNP and small-indel variant caller using deep neural networks. *Nature Biotechnology*, 36(10): 983–987.</sup> 

<sup>**3.** Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. *PLoS ONE*, 12(5): e0177459.</sup> 

<sup>**4.** Anaconda Software Distribution. Anaconda Inc. (2016). [https://docs.anaconda.com/](https://docs.anaconda.com/)</sup> 

<sup>**5.** Koster J, Rahmann S (2018). Snakemake – a scalable bioinformatics workflow engine. *Bioinformatics*, 34(20): 3600–3600.</sup>


