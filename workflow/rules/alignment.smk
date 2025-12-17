rule Trimming:
    input:
        unpack(get_fastq)
    output:
        r1 = f"{workpath}/trimmed/{{sample}}_1.fastq.gz",
        r2 = f"{workpath}/trimmed/{{sample}}_2.fastq.gz",
        json = f"{workpath}/trimmed/{{sample}}.json",
        html = f"{workpath}/trimmed/{{sample}}.html"
    log:
        "logs/alignment/trimming/{{sample}}.log"
    threads: 4
    conda:
        "../envs/trim.yaml"
    shell:
        """
        fastp -g -x -w {threads} \
            -D --dup_calc_accuracy 6 \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            -h {output.html} -j {output.json} -R {wildcards.sample} \
            > {log} 2>&1
        """

rule map_reads:
    input:
        unpack(get_map_reads_input),
        reference = f"{reference_dir}/genome/hg38.fa",
        fai = f"{reference_dir}/genome/hg38.fa.fai",
        dict = f"{reference_dir}/genome/hg38.dict",
        bwa_index = expand(f"{reference_dir}/genome/hg38.fa{{ext}}", ext=[".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"])
    output:
        f"{workpath}/sam/{{sample}}_raw.sam", 
    log:
        "logs/alignment/map_reads/{{sample}}.log"
    params:
        read_group = get_read_group,
    threads: 8
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa-mem2 mem -R '{params.read_group}' -t {threads} {input.reference} \
        {input.r1} {input.r2} > {output}  2> {log}
        """

rule mark_duplicates:
    input:
        f"{workpath}/sam/{{sample}}_raw.sam"
    output:
        raw_bam = f"{workpath}/bam/{{sample}}_raw.bam", 
        flagstat = f"{workpath}/bam/{{sample}}_raw.flagstat.txt", 
        dedup_bam = f"{workpath}/bam/{{sample}}_dedup.bam", 
        bam = f"{workpath}/bam/{{sample}}.bam", 
        bam_idx = f"{workpath}/bam/{{sample}}.bam.bai"
    log:
        "logs/alignment/mark_duplicates/{{sample}}.log"
    resources:
        mem_mb = 8000
    threads: 8
    conda:
        "../envs/sambamba.yaml"
    shell:
        """
        sambamba view -p -t {threads} -l 9 -S {input} -f bam -o {output.raw_bam} > {log} 2>&1
        sambamba flagstat -p -t {threads} {output.raw_bam} > {output.flagstat} 2>> {log}
        sambamba markdup -r -p -t {threads} -l 9 {output.raw_bam} {output.dedup_bam} >> {log} 2>&1
        sambamba sort -m {resources.mem_mb}M -p -t {threads} -l 9 {output.dedup_bam} -o {output.bam} >> {log} 2>&1
        sambamba index -p -t {threads} {output.bam} >> {log} 2>&1
        """
