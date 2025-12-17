if solo_analysis:
    rule rename_chromosomes_solo:
        input:
            filter_vcf = f"{workpath}/vcf/{{proband}}_merged.pass.vcf.gz",
            filter_vcf_tbi = f"{workpath}/vcf/{{proband}}_merged.pass.vcf.gz.tbi",
            chr_map = f"{reference_dir}/dbSNP/GCF_000001405.40_GRCh38.p14_assembly_report_revised_snpsift.chrnames"
        output:
            converted_vcf = f"{workpath}/vcf/{{proband}}-converted-solo.vcf.gz",
            converted_vcf_tbi = f"{workpath}/vcf/{{proband}}-converted-solo.vcf.gz.tbi"
        conda:
            "../envs/bcf.yaml"
        threads: 8
        log:
            "logs/annotation_solo/rename_chromosomes_solo/{proband}.log"
        shell:
            r"""
            (
                bcftools annotate --rename-chrs {input.chr_map} \
                    --threads {threads} -Oz -o {output.converted_vcf} {input.filter_vcf}
                tabix -p vcf {output.converted_vcf}
            ) > {log} 2>&1
            """


    rule annotate_vep_solo:
        input:
            converted_vcf = f"{workpath}/vcf/{{proband}}-converted-solo.vcf.gz",
            converted_vcf_tbi = f"{workpath}/vcf/{{proband}}-converted-solo.vcf.gz.tbi",
            cache_dir = f"{reference_dir}/vep"
        output:
            vep_vcf = f"{workpath}/vcf/{{proband}}-vep-solo.vcf.gz",
            vep_vcf_tbi = f"{workpath}/vcf/{{proband}}-vep-solo.vcf.gz.tbi"
        params:
            version = config["vep"]["version"]
        threads: 8
        conda:
            "../envs/vep.yaml"
        log:
            "logs/annotation_solo/annotate_vep_solo/{proband}.log"
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

    rule annotate_dbSNP_clinVar_solo:
        input:
            annotated_vep_vcf = f"{workpath}/vcf/{{proband}}-vep-solo.vcf.gz",
            annotated_vep_tbi = f"{workpath}/vcf/{{proband}}-vep-solo.vcf.gz.tbi",
            dbSNP = f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf.gz",
            dbSNP_idx = f"{reference_dir}/dbSNP/Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf.gz.tbi",
            ClinVar = expand("{reference_dir}/ClinVar/clinvar_{version}.vcf.gz", version=config["ClinVar"]["version"], reference_dir=reference_dir),
            ClinVar_idx = expand("{reference_dir}/ClinVar/clinvar_{version}.vcf.gz.tbi", version=config["ClinVar"]["version"], reference_dir=reference_dir),

        output:
            annotated_dbSNP_vcf = f"{workpath}/vcf/{{proband}}-vep-dbSNP-solo.vcf.gz",
            annotated_dbSNP_tbi = f"{workpath}/vcf/{{proband}}-vep-dbSNP-solo.vcf.gz.tbi",
            annotated_clinVar_vcf = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-solo.vcf.gz",
            annotated_clinVar_tbi = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-solo.vcf.gz.tbi"
        conda:
            "../envs/snpsift.yaml"
        log:
            "logs/annotation_solo/annotate_dbSNP_clinVar_solo/{proband}.log"
        shell:
            r"""
            (
                SnpSift annotate -v {input.dbSNP} {input.annotated_vep_vcf} | bgzip -c > {output.annotated_dbSNP_vcf}
                tabix -p vcf {output.annotated_dbSNP_vcf}

                SnpSift annotate -v {input.ClinVar} {output.annotated_dbSNP_vcf} | bgzip -c > {output.annotated_clinVar_vcf}
                tabix -p vcf {output.annotated_clinVar_vcf}
            ) > {log} 2>&1
            """

    rule annotate_dbNSFP_solo:
        input:
            annotated_clinVar_vcf = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-solo.vcf.gz",
            annotated_clinVar_tbi = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-solo.vcf.gz.tbi",
            dbnsfp = expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz", version=config["dbNSFP"]["version"], reference_dir=reference_dir),
            index = expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz.tbi", version=config["dbNSFP"]["version"], reference_dir=reference_dir),
            data_types = expand("{reference_dir}/dbNSFP/dbNSFPv{version}_custombuild.gz.data_types", version=config["dbNSFP"]["version"], reference_dir=reference_dir)    
        output:
            annotated_dbNSFP_vcf = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-dbNSFP_annotated-solo.vcf.gz",
            annotated_dbNSFP_tbi = f"{workpath}/vcf/{{proband}}-vep-dbSNP-ClinVar-dbNSFP_annotated-solo.vcf.gz.tbi", 
        conda:
            "../envs/snpsift.yaml"
        log:
            "logs/annotation_solo/annotate_dbNSFP_solo/{proband}.log"
        shell:
            r"""
            (
                SnpSift dbnsfp -v -db {input.dbnsfp} {input.annotated_clinVar_vcf} \
                | bgzip -c > {output.annotated_dbNSFP_vcf}
                tabix -p vcf {output.annotated_dbNSFP_vcf}
            ) > {log} 2>&1
            """
