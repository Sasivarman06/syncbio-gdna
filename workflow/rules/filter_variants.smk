if solo_analysis:
    # Pre-Annotation Filtration for Solo ---
    rule gatk_variant_filtration_solo:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_merged.vcf.gz",
            merged_idx = f"{workpath}/vcf/{{proband}}_merged.vcf.gz.tbi",
            reference = f"{reference_dir}/genome/hg38.fa",
            rederence_dict = f"{reference_dir}/genome/hg38.dict",
            reference_fai = f"{reference_dir}/genome/hg38.fa.fai"
        output:
            filtered_vcf = f"{workpath}/vcf/{{proband}}_merged.filtered.vcf.gz",
            filtered_vcf_tbi = f"{workpath}/vcf/{{proband}}_merged.filtered.vcf.gz.tbi"
        log:
            "logs/filter_variants/gatk_variant_filtration_solo/{{proband}}.log"
        conda:
            "../envs/gatk.yaml"
        shell:
            """
            gatk --java-options "-Xmx4g" VariantFiltration \
                -R {input.reference} \
                -V {input.vcf} \
                -O {output.filtered_vcf} \
                --filter-name "LowQual_QUAL" --filter-expression "QUAL < 30.0" \
                --filter-name "LowQual_QD"   --filter-expression "QD < 2.0" \
                --filter-name "LowQual_FS"   --filter-expression "FS > 60.0" \
                --filter-name "LowQual_MQ"   --filter-expression "MQ < 40.0" \
                --filter-name "LowQual_SOR"  --filter-expression "SOR > 3.0" \
                --filter-name "LowQual_DP"   --filter-expression "DP < 10" \
                > {log} 2>&1
            """

    rule gatk_select_pass_solo:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_merged.filtered.vcf.gz",
            vcf_idx = f"{workpath}/vcf/{{proband}}_merged.filtered.vcf.gz.tbi",
            reference = f"{reference_dir}/genome/hg38.fa",
            reference_dict = f"{reference_dir}/genome/hg38.dict",
            reference_fai = f"{reference_dir}/genome/hg38.fa.fai"
        output:
            final_vcf = f"{workpath}/vcf/{{proband}}_merged.pass.vcf.gz",
            final_vcf_tbi = f"{workpath}/vcf/{{proband}}_merged.pass.vcf.gz.tbi"
        log:
            "logs/filter_variants/gatk_select_pass_solo/{{proband}}.log"
        conda:
            "../envs/gatk.yaml"
        shell:
            """
            gatk --java-options "-Xmx4g" SelectVariants \
                -R {input.reference} \
                -V {input.vcf} \
                -O {output.final_vcf} \
                --exclude-filtered \
                > {log} 2>&1
            """

if trio_analysis:
    # Pre-Annotation Filtration for Trio ---
    rule gatk_variant_filtration_trio:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_trio_merged.vcf.gz",
            merged_idx = f"{workpath}/vcf/{{proband}}_trio_merged.vcf.gz.tbi",
            reference = f"{reference_dir}/genome/hg38.fa",
            rederence_dict = f"{reference_dir}/genome/hg38.dict",
            reference_fai = f"{reference_dir}/genome/hg38.fa.fai"
        output:
            filtered_vcf = f"{workpath}/vcf/{{proband}}_trio_merged.filtered.vcf.gz",
            filtered_vcf_tbi = f"{workpath}/vcf/{{proband}}_trio_merged.filtered.vcf.gz.tbi"
        log:
            "logs/filter_variants/gatk_variant_filtration_trio/{{proband}}.log"
        conda:
            "../envs/gatk.yaml"
        shell:
            """
            gatk --java-options "-Xmx4g" VariantFiltration \
                -R {input.reference} \
                -V {input.vcf} \
                -O {output.filtered_vcf} \
                --filter-name "LowQual_QUAL" --filter-expression "QUAL < 30.0" \
                --filter-name "LowQual_QD"   --filter-expression "QD < 2.0" \
                --filter-name "LowQual_FS"   --filter-expression "FS > 60.0" \
                --filter-name "LowQual_MQ"   --filter-expression "MQ < 40.0" \
                --filter-name "LowQual_SOR"  --filter-expression "SOR > 3.0" \
                --filter-name "LowQual_DP"   --filter-expression "DP < 10" \
                > {log} 2>&1
            """

    rule gatk_select_pass_trio:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_trio_merged.filtered.vcf.gz",
            vcf_idx = f"{workpath}/vcf/{{proband}}_trio_merged.filtered.vcf.gz.tbi",
            reference = f"{reference_dir}/genome/hg38.fa",
            reference_dict = f"{reference_dir}/genome/hg38.dict",   
            reference_fai = f"{reference_dir}/genome/hg38.fa.fai"
        output:
            final_vcf = f"{workpath}/vcf/{{proband}}_trio_merged.pass.vcf.gz",
            final_vcf_tbi = f"{workpath}/vcf/{{proband}}_trio_merged.pass.vcf.gz.tbi"
        log:
            "logs/filter_variants/gatk_select_pass_trio/{{proband}}.log"
        conda:
            "../envs/gatk.yaml"
        shell:
            """
            gatk --java-options "-Xmx4g" SelectVariants \
                -R {input.reference} \
                -V {input.vcf} \
                -O {output.final_vcf} \
                --exclude-filtered \
                > {log} 2>&1
            """