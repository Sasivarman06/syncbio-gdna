from pathlib import Path
import pandas as pd
import re 
import argparse
 
 
# Generate a temporary samplesheet for single-sample runs   
def generate_temp_samplesheet(fq1, fq2=None, output_dir="config", phenotype=None):
    # get full filename (no dirs)
    fname = Path(fq1).name  
    # remove extensions like .fastq.gz, .fq.gz, .fastq, .fq
    fname = re.sub(r'\.(fastq|fq)(\.gz)?$', '', fname)
    # strip read suffix (_1, _2, _R1, _R2, _R1_001, etc.)
    sample_id = re.sub(r'(_R?[12](?:_\d{3})?)$', '', fname)
    temp_csv = Path(output_dir) / "temp_samplesheet.csv"
    data = {
        "sample": [sample_id],
        "relation": ["proband"],
        "fq1": [str(fq1)],
        "fq2": [str(fq2) if fq2 else pd.NA],
        "bam": [pd.NA],
        "platform": ["ILLUMINA"],
        "case_id": [sample_id],
        "clinical_data": [str(phenotype) if phenotype else pd.NA]
    }
    pd.DataFrame(data).to_csv(temp_csv, index=False)
    return temp_csv

def generate_temp_bamsheet(bam, output_dir="config", phenotype=None):
    sample_id = Path(bam).stem.replace(".bam", "")
    temp_csv = Path(output_dir) / "temp_bamsheet.csv"
    data = {
        "sample": [sample_id],
        "relation": ["proband"],
        "fq1": [pd.NA],
        "fq2": [pd.NA],
        "bam": [str(bam)],
        "platform": ["ILLUMINA"],
        "case_id": [sample_id],
        "clinical_data": [str(phenotype) if phenotype else pd.NA]
    }
    pd.DataFrame(data).to_csv(temp_csv, index=False)
    return temp_csv


# Common helper functions shared across the entire workflow
def provided(samplelist, condition):
    """
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is not met (i.e. False) then it will not try to run that rule.
    """

    if not condition:
        # If condition is False, 
        # returns an empty list 
        # to prevent rule from 
        # running
        samplelist = []

    return samplelist


def ignore(samplelist, condition):
    """
    Filters samplelist based on condition.
    condition can be:
    - A list of samples to exclude
    - A boolean Series (True = exclude)
    - A boolean value (True = exclude all)
    """
    if isinstance(condition, (list, pd.Index)):
        exclude = set(condition)
        return [s for s in samplelist if s not in exclude]
    elif isinstance(condition, pd.Series):
        return samplelist[~condition]
    elif condition:
        # If condition is True (non-empty), return empty list
        return []
    return samplelist

def valid_bam(path: str) -> Path:
    p = Path(path)
    if not p.suffix == ".bam":
        raise argparse.ArgumentTypeError(f"{path} is not a valid BAM file (must end with .bam)")
    return p

def valid_csv(path: str) -> Path:
    p = Path(path)
    if not p.suffix == ".csv":
        raise argparse.ArgumentTypeError(f"{path} is not a valid CSV file (must end with .csv)")
    return p


def valid_fastq(path: str) -> Path:
    """
    Validate that a file is a FASTQ file (.fastq or .fastq.gz). which is passed as a command-line argument.
    """
    p = Path(path)
    if not (p.suffix == ".fastq" or (p.suffixes == [".fastq", ".gz"]) or
            p.suffix == ".fq" or (p.suffixes == [".fq", ".gz"])):
        raise argparse.ArgumentTypeError(
            f"{path} is not a valid FASTQ file (must end with .fastq or .fq or .fq.gz or .fastq.gz)"
        )
    return p
