if trio_analysis:
    rule genome_sdf:
        input:
            reference=f"{reference_dir}/genome/hg38.fa",
            reference_idx=f"{reference_dir}/genome/hg38.fa.fai",
            reference_dict=f"{reference_dir}/genome/hg38.dict"
        output:
            sdf=directory(f"{reference_dir}/sdf/hg38.sdf")
        conda:
            "../envs/rtg-tools.yaml"
        log:
            "logs/mendelian_annotation/genome_sdf/genome_sdf.log"
        shell:
            """
            rtg format -o {output.sdf} {input.reference} > log 2>&1
            """

    rule vcf_stats:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_trio_merged.pass.vcf.gz",
            vcf_idx = f"{workpath}/vcf/{{proband}}_trio_merged.pass.vcf.gz.tbi"
        output:
            stats=f"{workpath}/vcf/vcf_stats/{{proband}}_trio.vcfstats.txt"
        conda:
            "../envs/rtg-tools.yaml"
        log:
            "logs/mendelian_annotation/vcf_stats/{proband}.log"
        shell:
            """
            rtg vcfstats {input.vcf} > {output.stats} 2> {log}
            """

    rule mendelian_annotate:
        input:
            vcf=f"{workpath}/vcf/{{proband}}_trio_merged.pass.vcf.gz",
            vcf_idx=f"{workpath}/vcf/{{proband}}_trio_merged.pass.vcf.gz.tbi",
            sdf="resources/sdf/hg38.sdf/",
            pedigree=config["ped"]
        output:
            annotated_vcf=f"{workpath}/vcf/{{proband}}_trio_annotated.output.vcf.gz",
            annotated_vcf_idx=f"{workpath}/vcf/{{proband}}_trio_annotated.output.vcf.gz.tbi",
            mendelian_stats=f"{workpath}/vcf/{{proband}}_trio.mendelian.txt"
        conda:
            "../envs/rtg-tools.yaml"
        log:
            "logs/mendelian_annotation/mendelian_annotate/{proband}.log"
        shell:
            """
            rtg mendelian -i {input.vcf} -o {output.annotated_vcf} \
                          --pedigree {input.pedigree} -t {input.sdf} \
                          > {output.mendelian_stats} 2> {log}
            """
