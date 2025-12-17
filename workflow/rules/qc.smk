rule Fastqc:
    input:
        unpack(get_map_reads_input)
    output:
        r1_zip = f"{workpath}/qc/fastqc/{{sample}}_1_fastqc.zip",
        r2_zip = f"{workpath}/qc/fastqc/{{sample}}_2_fastqc.zip",
        r1_html = f"{workpath}/qc/fastqc/{{sample}}_1_fastqc.html",
        r2_html = f"{workpath}/qc/fastqc/{{sample}}_2_fastqc.html",
    params:
        outdir = f"{workpath}/qc/fastqc/",
    log:
        "logs/qc/Fastqc/{{sample}}.log"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} {input} -o {params.outdir} > {log} 2>&1
        """

rule samtools_stats:
    input:
        bam=lambda w: get_bam(w.sample),
    output:
        stats=f"{workpath}/qc/samtools-stats/{{sample}}.txt",
    log:
        "logs/qc/Samtools_stats/{{sample}}.log"
    threads: 2
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools stats {input.bam} > {output.stats} 2> {log}
        """

rule qualimap:
    input:
        bam=lambda w: get_bam(w.sample),
        reference=f"{reference_dir}/genome/hg38.fa",
    output:
        txt=f"{workpath}/qc/qualimap/{{sample}}/genome_results.txt",
        html=f"{workpath}/qc/qualimap/{{sample}}/qualimapReport.html",
    params:
        outdir=f"{workpath}/qc/qualimap/{{sample}}",
    log:
        "logs/qc/qualimap/{{sample}}.log"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        """
        qualimap bamqc -bam {input.bam} -outdir {params.outdir} -outformat HTML -nt {threads} --skip-duplicated > {log} 2>&1
        """

rule multiqc:
    input:
        expand("{workpath}/qc/fastqc/{sample}_1_fastqc.zip", workpath=workpath, sample=ignore(samples.index, skip_align)),
        expand("{workpath}/qc/fastqc/{sample}_2_fastqc.zip", workpath=workpath, sample=ignore(samples.index, skip_align)),
        expand("{workpath}/qc/samtools-stats/{sample}.txt", workpath=workpath, sample=samples.index),
        expand("{workpath}/qc/qualimap/{sample}/qualimapReport.html", workpath=workpath, sample=samples.index)
    output:
        report_dir = directory(f"{workpath}/qc/multiqc"),
        report = f"{workpath}/qc/multiqc/multiqc_report.html",
    log:
        "logs/qc/multiqc.log"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {input} -o {output.report_dir} > {log} 2>&1
        """
