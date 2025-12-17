rule SV_index:
    input:
        f"{reference_dir}/genome/hg38.fa"
    output:
        multiext(f"{reference_dir}/genome/hg38.fa", ".bwt", ".sa")
    log:
        "logs/call_SV/SV_index/hg38_index.log"
    conda:
        "../envs/bwa.yaml"
    threads: 4
    shell:
        """
        bwa index {input} > {log} 2>&1
        """


rule gridss_SV:
    input:
        proband_bam = lambda w: get_bam(w.proband),
        bam_idx = lambda w: get_bam(w.proband) + ".bai",
        reference = f"{reference_dir}/genome/hg38.fa",
        reference_idx = f"{reference_dir}/genome/hg38.fa.fai",
        sequence_dict = f"{reference_dir}/genome/hg38.dict",
        idx = multiext(f"{reference_dir}/genome/hg38.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")          
    output:
        vcf = f"{workpath}/SV/{{proband}}_output.vcf",
        vcf_idx = f"{workpath}/SV/{{proband}}_output.vcf.idx",
        assembly = f"{workpath}/SV/{{proband}}_assembly.bam",
        working_dir = directory(f"{workpath}/SV/{{proband}}_gridss_working"),
    threads: 4
    conda:
        "../envs/gridss.yaml"
    log:
        "logs/call_SV/gridss_SV/{{proband}}.log"
    shell:     
        """
        gridss \
          --reference {input.reference} \
          --output {output.vcf} \
          --assembly {output.assembly} \
          --workingdir {output.working_dir} \
          {input.proband_bam} \
          --threads {threads} > {log} 2>&1
        """


rule filter_vcf:
    input:
        vcf=f"{workpath}/SV/{{proband}}_output.vcf",
        vcf_idx=f"{workpath}/SV/{{proband}}_output.vcf.idx"
    output:
        filtered_vcf=f"{workpath}/SV/{{proband}}.filtered.vcf"
    log:
        "logs/call_SV/filter_vcf/{proband}.log"
    shell:
        """
        grep -E "#|PASS" {input.vcf} > {output.filtered_vcf} 2> {log}
        """


rule annotate_vep_SV:
    input:
        converted_vcf = f"{workpath}/SV/{{proband}}.filtered.vcf",
        cache_dir=f"{reference_dir}/vep"
    output:
        vep_vcf = f"{workpath}/SV/{{proband}}.filtered.vep.vcf",
    params:
        version=config["vep"]["version"]
    threads: 4
    conda:
        "../envs/vep.yaml"
    log:
        "logs/call_SV/annotate_vep_SV/{proband}.log"
    shell:
        """
        vep -i {input.converted_vcf} \
            -o {output.vep_vcf} \
            --vcf \
            --cache \
            --dir_cache {input.cache_dir} \
            --assembly GRCh38 \
            --everything \
            --fork {threads} \
            --force_overwrite \
            > {log} 2>&1
        """


rule annotate_snpsift:
    input:
        vep_vcf=f"{workpath}/SV/{{proband}}.filtered.vep.vcf",
        dbnsfp=expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz", version=config["dbNSFP"]["version"], reference_dir=reference_dir),
        index=expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz.tbi", version=config["dbNSFP"]["version"], reference_dir=reference_dir),
    output:
        dbnsfp_vcf=f"{workpath}/SV/{{proband}}.filtered.dbnsfp.vcf"
    log:
        "logs/call_SV/annotate_snpsift/{proband}.log"
    conda:
        "../envs/snpsift.yaml"
    shell:
        """
        SnpSift dbnsfp -v -db {input.dbnsfp} \
        {input.vep_vcf} > {output.dbnsfp_vcf} 2> {log}
        """