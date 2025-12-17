if trio_analysis:
    rule rename_chromosomes_trio:
        input:
            annotated_vcf = f"{workpath}/vcf/{{proband}}_trio_annotated.output.vcf.gz",
            annotated_vcf_idx = f"{workpath}/vcf/{{proband}}_trio_annotated.output.vcf.gz.tbi",
            chr_map = f"{reference_dir}/dbSNP/GCF_000001405.40_GRCh38.p14_assembly_report_revised_snpsift.chrnames"
        output:
            converted_vcf = f"{workpath}/vcf/{{proband}}-converted-trio.vcf.gz",
            converted_vcf_idx = f"{workpath}/vcf/{{proband}}-converted-trio.vcf.gz.tbi"
        threads: 8
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/annotation_trio/rename_chromosomes_trio/{proband}.log"
        shell:
            r"""
            (
                bcftools annotate --rename-chrs {input.chr_map} \
                    --threads {threads} -Oz -o {output.converted_vcf} {input.annotated_vcf}
                tabix -p vcf {output.converted_vcf}
            ) > {log} 2>&1
            """


    rule annotate_vep_trio:
        input:
            converted_vcf = f"{workpath}/vcf/{{proband}}-converted-trio.vcf.gz",
            converted_vcf_idx = f"{workpath}/vcf/{{proband}}-converted-trio.vcf.gz.tbi",
            cache_dir = f"{reference_dir}/vep",
        output:
            vep_vcf = f"{workpath}/vcf/{{proband}}-vep-trio.vcf.gz",
            vep_vcf_idx = f"{workpath}/vcf/{{proband}}-vep-trio.vcf.gz.tbi"
        params:
            version = config["vep"]["version"]
        threads: 8
        conda:
            "../envs/vep.yaml"
        log:
            "logs/annotation_trio/annotate_vep_trio/{proband}.log"
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
                --compress_output bgzip \
                > {log} 2>&1
            tabix -p vcf {output.vep_vcf} >> {log} 2>&1
            """

    rule annotate_dbSNP_clinVar_trio:
        input:
            annotated_vep_vcf = f"{workpath}/vcf/{{proband}}-vep-trio.vcf.gz",
            annotated_vep_tbi = f"{workpath}/vcf/{{proband}}-vep-trio.vcf.gz.tbi",
            dbSNP = f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf.gz",
            dbSNP_idx = f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf.gz.tbi",
            ClinVar = expand("{reference_dir}/ClinVar/clinvar_{version}.vcf.gz", version=config["ClinVar"]["version"], reference_dir=reference_dir),
            ClinVar_idx = expand("{reference_dir}/ClinVar/clinvar_{version}.vcf.gz.tbi", version=config["ClinVar"]["version"], reference_dir=reference_dir)
        output:
            annotated_dbSNP_vcf = f"{workpath}/vcf/{{proband}}-vep-dbSNP-trio.vcf.gz",
            annotated_dbNSFP_vcf_idx = f"{workpath}/vcf/{{proband}}-vep-dbSNP-trio.vcf.gz.tbi",
            annotated_clinVar_vcf = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-trio.vcf.gz",
            annotated_clinVar_vcf_idx = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-trio.vcf.gz.tbi",
        conda:
            "../envs/snpsift.yaml"
        log:
            "logs/annotation_trio/annotate_dbSNP_clinVar_trio/{proband}.log"
        shell:
            r"""
            (
                SnpSift annotate -v {input.dbSNP} {input.annotated_vep_vcf} | bgzip -c > {output.annotated_dbSNP_vcf}
                tabix -p vcf {output.annotated_dbSNP_vcf}

                SnpSift annotate -v {input.ClinVar} {output.annotated_dbSNP_vcf} | bgzip -c > {output.annotated_clinVar_vcf}
                tabix -p vcf {output.annotated_clinVar_vcf}
            ) > {log} 2>&1
            """


    rule annotate_dbNSFP_trio:
        input:
            annotated_clinVar_vcf = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-trio.vcf.gz",
            annotated_clinVar_tbi = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-trio.vcf.gz.tbi",
            dbnsfp =  expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz", version=config["dbNSFP"]["version"], reference_dir=reference_dir),
            dbnsfp_idx = expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz.tbi", version=config["dbNSFP"]["version"], reference_dir=reference_dir)
        output:
            annotated_dbNSFP_vcf = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-dbNSFP_annotated-trio.vcf.gz",
            annotated_dbNSFP_vcf_idx = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-dbNSFP_annotated-trio.vcf.gz.tbi"
        conda:
            "../envs/snpsift.yaml"
        log:
            "logs/annotation_trio/annotate_dbNSFP_trio/{proband}.log"
        shell:
            r"""
            (
                SnpSift dbnsfp -v -db {input.dbnsfp} {input.annotated_clinVar_vcf} \
                | bgzip -c > {output.annotated_dbNSFP_vcf}
                tabix -p vcf {output.annotated_dbNSFP_vcf}
            ) > {log} 2>&1
            """
