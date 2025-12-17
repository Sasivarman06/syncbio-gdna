rule get_vep_cache:
    output:
        directory(f"{reference_dir}/vep")
    params:
        version=config["vep"]["version"],
        species="homo_sapiens",
        assembly="GRCh38"
    conda:
        "../envs/vep.yaml"
    log:
        "logs/database/get_vep_cache.log"
    shell:
        """
        vep_install -a cf \
                    -s {params.species} \
                    -y {params.assembly} \
                    -c {output} > {log} 2>&1
        """

rule get_dbNSFP:
    output:
        expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz", version=config["dbNSFP"]["version"], reference_dir=reference_dir)
    params:
        database=config["dbNSFP"]["database"]
    log:
        "logs/database/get_dbNSFP.log"
    shell:
        """
        wget {params.database} -O {output} > {log} 2>&1 
        """

rule extract_dbNSFP:
    input:
        expand("{reference_dir}/dbNSFP/dbNSFP{version}.zip", version=config["dbNSFP"]["version"], reference_dir=reference_dir)
    output:
        expand("{reference_dir}/dbNSFP/dbNSFP{version}_variant.chr{chr}.gz", chr=[*range(1,23), "X", "Y", "M"], version=config["dbNSFP"]["version"], reference_dir=reference_dir),
    log:
        "logs/database/extract_dbNSFP.log"
    shell:
        """
        unzip {input} -d resources/dbNSFP/ > {log} 2>&1
        """

rule header_dbNSFP:
    input:
        chr_file=expand("{reference_dir}/dbNSFP/dbNSFP{version}_variant.chr1.gz",version=config["dbNSFP"]["version"], reference_dir=reference_dir)
    output:
        header=expand("{reference_dir}/dbNSFP/header_{version}.gz", version=config["dbNSFP"]["version"], reference_dir=reference_dir)
    conda:
        "../envs/bgzip.yaml"
    threads: 6 
    log:
        "logs/database/header_dbNSFP.log"
    shell:
        """
        set -ex
        {{ zcat {input.chr_file} | head -n 1 | bgzip --threads {threads} > {output.header}; }} || true  > {log} 2>&1
        """

rule merge_dbNSFP:
    input:
        chr_files = expand(
            f"{reference_dir}/dbNSFP/dbNSFP{{version}}_variant.chr{{chr}}.gz",
            chr=[*range(1,23), "X", "Y", "M"],
            version=config["dbNSFP"]["version"],
        ),
        header = f"{reference_dir}/dbNSFP/header_{{version}}.gz",

    output:    
        merged = f"{reference_dir}/dbNSFP/dbNSFPv{{version}}_custom.gz", 
        final  = f"{reference_dir}/dbNSFP/dbNSFPv{{version}}_custombuild.gz", 

    threads: 6
    log:
        f"logs/database/merge_dbNSFP/merge_dbNSFP_{{version}}.log"

    shell:
        """
        (
            zcat {input.chr_files} | grep -v '#chr' | bgzip --threads {threads} > {output.merged}
            cat {input.header} {output.merged} > {output.final}
        ) > {log} 2>&1
        """

rule count_unique_chromosomes:
    input:
      expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz", version=config["dbNSFP"]["version"], reference_dir=reference_dir)
    output:
      expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz.tbi", version=config["dbNSFP"]["version"], reference_dir=reference_dir)
    log:
        "logs/database/count_unique_chromosomes.log"
    shell:
        """
        zcat {input} | awk '!/^#/' | cut -f 1 | sort | uniq | wc -l > {log}
        tabix -s 1 -b 2 -e 2 {input} >> {log} 2>&1
        """

rule get_dbSNP:
    output:
        dbsnp_vcf=f"{reference_dir}/dbSNP/GCF_000001405.40.gz",
        dbsnp_index=f"{reference_dir}/dbSNP/GCF_000001405.40.gz.tbi",
        assembly_report=f"{reference_dir}/dbSNP/GCF_000001405.40_GRCh38.p14_assembly_report.txt",
        chr_map=f"{reference_dir}/dbSNP/GCF_000001405.40_GRCh38.p14_assembly_report_revised_snpsift.chrnames"
    params:
        dbSNP_url=config["dbSNP"]["ref"],
        index_url=config["dbSNP"]["index"],
        assembly_report_url=config["dbSNP"]["assembly_report"],

    conda:
        "../envs/aria2.yaml"
    log:
        "logs/database/get_dbSNP.log"
    shell:
        """
        mkdir -p $(dirname {output.dbsnp_vcf}) > {log} 2>&1
        wget {params.dbSNP_url} -O {output.dbsnp_vcf} >> {log} 2>&1 || exit 1
        wget {params.index_url} -O {output.dbsnp_index} >> {log} 2>&1 || exit 1
        wget {params.assembly_report_url} -O {output.assembly_report} >> {log} 2>&1 || exit 1
        awk '$1 !~ /^#/ {{print $7"\\t"$1}}' {output.assembly_report} > {output.chr_map}
        """

rule rename_chr:
    input:
        dbSNP_vcf=f"{reference_dir}/dbSNP/GCF_000001405.40.gz",
        chr_map=f"{reference_dir}/dbSNP/GCF_000001405.40_GRCh38.p14_assembly_report_revised_snpsift.chrnames"
    output:
        renamed_vcf=f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.vcf.gz"
    conda:
        "../envs/bcf.yaml"
    log:
        "logs/database/rename_chr.log"
    shell:
        """
        bcftools annotate \
          --rename-chrs {input.chr_map} \
          --threads $(nproc) -Oz \
          -o {output.renamed_vcf} \
          {input.dbSNP_vcf} > {log} 2>&1
        """

rule fix_forbidden_chars:
    input:
        renamed_vcf=f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.vcf.gz"
    output:
        fixed_vcf=f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf"
    conda:
        "../envs/snpsift.yaml"
    log:
        "logs/database/fix_forbidden_chars.log"
    shell:
        r"""
        SnpSift vcfCheck {input.renamed_vcf} 2>&1 | grep "INFO field" | cut -f 2 -d "'" | sort | uniq -c > {log}
        zcat {input.renamed_vcf} \
            | sed 's/\&base_change=/\&base_change%3D/g' \
            | sed 's/A=;/A%3D;/' \
            | sed 's/C=;/C%3D;/' \
            | sed 's/G=;/G%3D;/' \
            | sed 's/T=;/T%3D;/' \
            | sed 's/=,;/%3D,/' \
            >> {log} 2>&1
        """

rule compress_and_index:
    input:
        fixed_vcf=f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf"
    output:
        compressed_vcf=f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf.gz",
        vcf_index=f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf.gz.tbi"
    conda:
        "../envs/snpsift.yaml"
    log:
        "logs/database/compress_and_index.log"
    shell:
        """
        bgzip -c {input.fixed_vcf} > {output.compressed_vcf}
        tabix {output.compressed_vcf}
        SnpSift vcfCheck {output.compressed_vcf} 2>&1 | grep "INFO field" | cut -f 2 -d "'" | sort | uniq -c > {log}
        """

rule get_ClinVar:
    output:
        clinvar_vcf=expand("{reference_dir}/ClinVar/clinvar_{version}.vcf.gz", version=config["ClinVar"]["version"], reference_dir=reference_dir),
        vcf_index=expand("{reference_dir}/ClinVar/clinvar_{version}.vcf.gz.tbi", version=config["ClinVar"]["version"], reference_dir=reference_dir)
    params:
        data=config["ClinVar"]["data"]
    conda:
        "../envs/snpsift.yaml"
    log:
        "logs/database/get_ClinVar.log"
    shell:
        """
        wget {params.data} -O {output.clinvar_vcf} > {log} 2>&1
        tabix {output.clinvar_vcf} >> {log} 2>&1
        """

rule download_exomiser:
    output:
        exomiser_db=expand("{reference_dir}/exomiser/exomiser-cli-{version}-distribution.zip", version=config["exomiser"]["jar"]["version"], reference_dir=reference_dir),
        assembly=expand("{reference_dir}/exomiser/{version}_hg38.zip", version=config["exomiser"]["assembly"]["version"], reference_dir=reference_dir),
        phenotype_db=expand("{reference_dir}/exomiser/{version}_phenotype.zip", version=config["exomiser"]["phenotype_db"]["version"], reference_dir=reference_dir)
    params:
        jar=config["exomiser"]["jar"]["url"],
        assembly=config["exomiser"]["assembly"]["url"],
        phenotype_db=config["exomiser"]["phenotype_db"]["url"],
        outdir=f"{reference_dir}/exomiser"
    log:
        "logs/database/download_exomiser.log"
    shell:
        """
        wget {params.jar} -O {output.exomiser_db} > {log} 2>&1
        wget {params.assembly} -O {output.assembly} >> {log} 2>&1
        wget {params.phenotype_db} -O {output.phenotype_db} >> {log} 2>&1
        unzip -o {output.exomiser_db} -d {params.outdir} >> {log} 2>&1
        unzip -o {output.assembly} -d {params.outdir} >> {log} 2>&1
        unzip -o {output.phenotype_db} -d {params.outdir} >> {log} 2>&1
        """