import pandas as pd
from snakemake.utils import validate

# ---- Load samplesheet ----
configfile: "config/config.yaml"
samples = pd.read_csv(config["samples"]).set_index("sample", drop=False)

validate(samples, schema="../schemas/samples.schema.yaml")

# ---- Detect modes ----
def detect_modes(samples: pd.DataFrame):
    """
    Detect which analysis modes exist (solo/trio), list probands, 
    and map trio probands to their parents.
    Solo samples are ignored for trio mapping.
    """
    modes = set()
    probands = []  
    trio_map = {}  

    if "relation" in samples.columns and "case_id" in samples.columns:
        # Group by case_id
        for case_id, group in samples.groupby("case_id", dropna=False):
            relations = set(group["relation"])
            case_probands = group.loc[group["relation"] == "proband", "sample"].tolist()

            # Only full trios go into trio_map
            if {"proband", "mother", "father"}.issubset(relations):
                for proband in case_probands:
                    father = group.loc[group["relation"] == "father", "sample"].tolist()
                    mother = group.loc[group["relation"] == "mother", "sample"].tolist()
                    trio_map[proband] = {
                        "father": father[0],
                        "mother": mother[0],
                    }
                modes.add("trio")
            else:
                # Probands without full trio → solo
                probands.extend(case_probands)
                modes.add("solo")
    else:
        # No case_id or relation column → all solo
        probands = samples["sample"].tolist()
        modes.add("solo")

    # Separate solo and trio probands
    solo_probands = [p for p in probands if p not in trio_map]
    trio_probands = list(trio_map.keys())
    solo_trio_probands = solo_probands + trio_probands

    return  modes, solo_probands, trio_probands, trio_map, solo_trio_probands


# ---- Run detection ----
modes, solo_probands, trio_probands, trio_map, solo_trio_probands = detect_modes(samples)

solo_analysis = "solo" in modes
trio_analysis = "trio" in modes

# ---- Wildcards ----
wildcard_constraints:
    sample="|".join(samples.index),
    proband="|".join(solo_probands + trio_probands),
    mother="|".join(samples[samples["relation"] == "mother"].index) if "relation" in samples and not samples[samples["relation"] == "mother"].empty else "",
    father="|".join(samples[samples["relation"] == "father"].index) if "relation" in samples and not samples[samples["relation"] == "father"].empty else ""


def get_fastq(wildcards):
    """Return dict with paired-end FASTQ files; raise error if missing."""
    fq1 = str(samples.loc[wildcards.sample, "fq1"]).strip()
    fq2 = str(samples.loc[wildcards.sample, "fq2"]).strip()

    if not fq1 or fq1.lower() == "nan":
        raise ValueError(f"No valid FASTQ1 for sample {wildcards.sample}")

    if not fq2 or fq2.lower() == "nan":
        raise ValueError(f"No valid FASTQ2 for sample {wildcards.sample}. Paired-end required.")

    return {"r1": fq1, "r2": fq2}


def get_read_group(wildcards):
    return "@RG\tID:{sample}\tSM:{sample}\tPL:{platform}".format(
        sample=wildcards.sample,
        platform=samples.loc[(wildcards.sample), "platform"],
    )

def get_bam(sample):
    try:
        bam = str(samples.loc[sample, "bam"]).strip()
        if bam and bam.lower() != "nan":
            return bam
    except KeyError:
        pass
    return f"{workpath}/bam/{sample}.bam"

def get_map_reads_input(wildcards):
    sample = wildcards.sample
    if sample in skip_align:
        # Alignment should not run → return empty, prevent rule execution
        return None

    if sample in skip_trim:
        return {"r1": samples.loc[sample,"fq1"], "r2": samples.loc[sample,"fq2"]}
    else:
        return {"r1": f"{workpath}/trimmed/{sample}_1.fastq.gz",
                "r2": f"{workpath}/trimmed/{sample}_2.fastq.gz"}


def get_clinical_data(proband):
    try:
        if proband not in samples.index:
            return None

        clinical_file = str(samples.loc[proband, "clinical_data"]).strip()

        if clinical_file and clinical_file.lower() != "nan":
            return clinical_file

        return None
    except KeyError:
        return None

def get_ped_file(proband):
    """Generate PED file content for a trio given the proband sample name."""
    if proband not in trio_map:
        return None

    father = trio_map[proband]["father"]
    mother = trio_map[proband]["mother"]

    ped_lines = [
        f"{proband}\t{proband}_family\t{mother}\t{father}\t1\t2",
        f"{mother}\t{proband}_family\t0\t0\t2\t1",
        f"{father}\t{proband}_family\t0\t0\t1\t1",
    ]

    return "\n".join(ped_lines)