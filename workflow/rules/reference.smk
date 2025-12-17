rule get_genome:
    output:
        f"{reference_dir}/genome/hg38.fa.gz"
    log:
        "logs/reference/get-genome.log"
    params:
        url=config["genome"]["url"]
    conda:
        "../envs/snpsift.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) > {log} 2>&1
        wget {params.url} -O {output} >> {log} 2>&1
        """

rule extract_genome:
    input:
        f"{reference_dir}/genome/hg38.fa.gz"
    output:
        f"{reference_dir}/genome/hg38.fa"
    log:
        "logs/reference/extract-genome.log"
    shell:
        """
        gunzip -c {input} > {output} 2> {log}
        """

rule genome_faidx:
    input:
        f"{reference_dir}/genome/hg38.fa"
    output:
        f"{reference_dir}/genome/hg38.fa.fai"
    log:
        "logs/reference/genome-faidx.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} > {log} 2>&1"

rule genome_dict:
    input:
        f"{reference_dir}/genome/hg38.fa"
    output:
        f"{reference_dir}/genome/hg38.dict"
    log:
        "logs/reference/genome_dict.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools dict {input} > {output} 2> {log}"

rule bwa_index:
    input:
        f"{reference_dir}/genome/hg38.fa"
    output:
        multiext(f"{reference_dir}/genome/hg38.fa", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
    log:
        "logs/reference/bwa_index.log"
    conda:
        "../envs/bwa.yaml"
    threads: 8
    resources:
        mem_mb=369000
    shell:
        """
        bwa-mem2 index {input} > {log} 2>&1
        """