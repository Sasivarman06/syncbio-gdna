
# <code>syncbio-gdna <b>run</b></code>

## 1. About 
The `syncbio-gdna` executable is composed of several inter-related sub commands. Please see `syncbio-gdna -h` for all available options.

This part of the documentation describes options and concepts for <code>syncbio-gdna <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running the SyncBio-gDNA pipeline for rare disease variant calling. 

Setting up the pipeline is fast and easy! In its most basic form, <code>syncbio-gdna <b>run</b></code> requires only an input and output specification.

## 2. Synopsis
```text
$ syncbio-gdna run [--help] \
      [--fq1 FQ1] [--fq2 FQ2] [--bam BAM] [--samplesheet SAMPLESHEET] \
      [--output OUTPUT] [--reference REFERENCE] [--ped PED] \
      [--intervals INTERVALS] [--call-sv] [--skip-trim] \
      [--dry-run] [--rerun-incomplete] [--cores CORES] \
      [--sif-cache SIF_CACHE] [--tmpdir TMPDIR]
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide an output directory to store results via `--output` argument and at least one input method.

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output /data/$USER/syncbio-gdna_results`

### 2.2 Input options (choose one method)

Each of the following arguments are part of a mutually exclusive group. Only one input method should be provided.

  `--fq1 FQ1`  
> **Forward reads (single-sample mode).**  
> *type: file*
> 
> Input FastQ file containing forward reads for single-sample analysis. This option must be used together with `--fq2`.
> 
> ***Example:*** `--fq1 sample1_R1.fastq.gz`

---  
  `--fq2 FQ2`  
> **Reverse reads (single-sample mode).**  
> *type: file*
> 
> Input FastQ file containing reverse reads for single-sample analysis. This option must be used together with `--fq1`.
> 
> ***Example:*** `--fq2 sample1_R2.fastq.gz`

---  
  `--bam BAM`  
> **Input BAM file (single-sample mode).**  
> *type: file*
> 
> Input BAM file for single-sample analysis. Use this option when starting from aligned data rather than FastQ files.
> 
> ***Example:*** `--bam sample1.aligned.bam`

---  
  `--samplesheet SAMPLESHEET`  
> **Path to samplesheet (multiple-sample mode).**  
> *type: CSV file*
> 
> CSV file containing sample information for multiple-sample analysis. This should be used when processing multiple samples simultaneously. The samplesheet format supports both solo and trio analysis, with trio relationships determined based on the `case_id` and `relation` fields. The samplesheet format is described in section 3.
> 
> ***Example:*** `--samplesheet sample_sheet.csv`

### 2.3 Analysis options

Each of the following arguments are optional, and do not need to be provided. 

  `--reference REFERENCE`            
> **Path to reference genome directory.**  
> *type: path*
> 
> Custom directory containing reference genome files. If not provided, the pipeline will use the default reference specified in the configuration.
>
> ***Example:*** `--reference /data/references/GRCh38`

---  
  `--ped PED`            
> **Path to pedigree file.**  
> *type: file*
> 
> Provides pedigree information for trio analysis. When provided, the pipeline will perform family-based analysis using the pedigree structure. Note: When using a samplesheet with trio information, this option is not required as relationships are determined from the samplesheet.
>
> ***Example:*** `--ped family.ped`

---  
  `--intervals INTERVALS`            
> **Genomic intervals for targeted variant calling.**  
> *type: file*
> 
> BED file or interval list defining target regions for variant calling. Useful for exome sequencing or panel data to restrict variant calling to specific genomic regions.
>
> ***Example:*** `--intervals target_regions.bed`

---  
  `--call-sv`            
> **Enable structural variant calling.**  
> *type: boolean flag*
> 
> When specified, the pipeline will perform structural variant (SV) detection in addition to SNV calling.
>
> ***Example:*** `--call-sv`

---  
  `--skip-trim`            
> **Skip read trimming step.**  
> *type: boolean flag*
> 
> Bypass the read trimming and quality control steps. Use this when input FastQ files are already pre-processed.
>
> ***Example:*** `--skip-trim`

---  
  `--rerun-incomplete`            
> **Force pipeline to re-run jobs with incomplete outputs.**  
> *type: boolean flag*
> 
> This flag forces the pipeline to re-execute jobs that have incomplete outputs, which is useful for resuming interrupted runs or fixing partial results.
>
> ***Example:*** `--rerun-incomplete`

### 2.4 Execution options

Each of the following arguments are optional, and do not need to be provided. 

  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything! Useful for testing the pipeline configuration before running.
>
> ***Example:*** `--dry-run`

---  
  `--cores CORES`  
> **CPU cores to use.**  
> *type: integer*  
> *default: 4*
> 
> Specifies the number of CPU cores to use for parallel execution. The default is 4 cores.
> 
> ***Example:*** `--cores 8`

---  
  `--sif-cache SIF_CACHE`  
> **Custom directory to store Singularity images.**  
> *type: path*  
>
> Sets a custom directory for storing Singularity container images. This can be useful for managing disk space or using a shared cache across multiple runs.
> 
> ***Example:*** `--sif-cache /shared/singularity_cache`

---  
  `--tmpdir TMPDIR`  
> **Directory to store temporary files.**  
> *type: path*  
>
> Specifies a directory for storing temporary files during pipeline execution. This can be useful when the default temporary directory has limited space.
> 
> ***Example:*** `--tmpdir /scratch/$USER/tmp`

### 2.5 Miscellaneous options  
Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Samplesheet Format

The samplesheet CSV file supports both solo and trio analysis. The pipeline automatically detects family relationships based on the `case_id` field - samples with the same `case_id` are treated as belonging to the same family.

### 3.1 Example Samplesheet
```csv
case_id,sample,relation,platform,sex,phenotype,fq1,fq2,bam
justhusky,hugelymodelbat,proband,illumina,1,2,test_data/hugelymodelbat_1_test.fastq.gz,test_data/hugelymodelbat_2_test.fastq.gz,
justhusky,earlycasualcaiman,father,illumina,1,1,test_data/earlycasualcaiman_1_test.fastq.gz,test_data/earlycasualcaiman_2_test.fastq.gz,
justhusky,slowlycivilbuck,mother,illumina,2,1,test_data/slowlycivilbuck_1_test.fastq.gz,test_data/slowlycivilbuck_2_test.fastq.gz,
```

### 3.2 Field Descriptions
- **case_id**: Unique identifier grouping family members
- **sample**: Unique sample identifier
- **relation**: Family role (proband, mother, father)
- **platform**: Sequencing technology (illumina, etc.)
- **sex**: Biological sex (1 = male, 2 = female)
- **phenotype**: Clinical status (1 = unaffected, 2 = affected)
- **fq1**: Path to forward read FastQ file (leave empty if using BAM)
- **fq2**: Path to reverse read FastQ file (leave empty if using BAM)
- **bam**: Path to BAM file (leave empty if using FastQ files)

## 4. Example Commands
```bash 
# Option 1: Dry run to test configuration
syncbio-gdna run \
    --fq1 sample4_R1.fastq.gz \
    --fq2 sample4_R2.fastq.gz \
    --output results/sample4 \
    --dry-run
    
# Option 2: Run pipeline with FastQ inputs (single sample)
syncbio-gdna run \
    --fq1 sample1_R1.fastq.gz \
    --fq2 sample1_R2.fastq.gz \
    --output results/sample1 \
    --cores 8

# Option 3: Run pipeline with BAM input (single sample)
syncbio-gdna run \
    --bam sample2.aligned.bam \
    --output results/sample2 \
    --sif-cache /shared/singularity_cache

# Option 4: Run pipeline with samplesheet (trio analysis)
syncbio-gdna run \
    --samplesheet trio_samplesheet.csv \
    --output results-trio \
    --call-sv \
    --tmpdir /scratch/$USER/tmp

# Option 5: Run pipeline with targeted regions
syncbio-gdna run \
    --samplesheet exome_samples.csv \
    --output results/exome_analysis \
    --intervals target_regions.bed \
    --cores 12

# Option 6: Run pipeline skipping trimming (pre-processed data)
syncbio-gdna run \
    --fq1 sample3_R1.fastq.gz \
    --fq2 sample3_R2.fastq.gz \
    --output results/sample3 \
    --skip-trim
```


