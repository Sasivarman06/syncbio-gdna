rule phenotype_mapping:
    input:
        clinical = lambda wildcards: get_clinical_data(wildcards.proband),
        ped = lambda wildcards: config["ped"] if wildcards.proband in trio_probands else [],
        vcf = lambda wildcards: ( 
            f"{workpath}/vcf/{wildcards.proband}-vep-dbSNP-ClinVar-dbNSFP_annotated-trio.vcf.gz"
            if wildcards.proband in trio_probands
            else f"{workpath}/vcf/{wildcards.proband}-vep-dbSNP-ClinVar-dbNSFP_annotated-solo.vcf.gz" 
        )
    output:
        yml = f"{workpath}/phenotype/{{proband}}_exomiser.yaml",
        hpo_list = f"{workpath}/phenotype/{{proband}}_transformed_hpo_set.txt",
        hpo_tsv = f"{workpath}/phenotype/{{proband}}_transformed_hpo_set.tsv",
        omim_tsv = f"{workpath}/phenotype/{{proband}}_transformed_omim_set.tsv",
        ddx_sim = f"{workpath}/phenotype/{{proband}}_differential_diagnosis_similarity.tsv",
        gene_sim = f"{workpath}/phenotype/{{proband}}_recommended_gene_similarity.tsv",
        disease_sim = f"{workpath}/phenotype/{{proband}}_recommended_disease_similarity.tsv",
        dendro_all = f"{workpath}/phenotype/{{proband}}_linkage_all_ddx.png",
        dendro_filt = f"{workpath}/phenotype/{{proband}}_linkage_filt_ddx.png",
        dendro_gene = f"{workpath}/phenotype/{{proband}}_linkage_gene.png",
        dendro_disease = f"{workpath}/phenotype/{{proband}}_linkage_disease.png"
    params:
        proband = "{proband}",
        threshold = 0.4,  
        differential = 10,     
        recommendation = 20,
        workpath = workpath,
        trio_analysis=trio_analysis,
        solo_analysis=solo_analysis,
        trio_probands = trio_probands,
    log:
        "logs/phenotype/phenotype_mapping/{proband}-phenotype_mapping.log"
    script:
        "../scripts/phenotype_mapping.py"

rule run_exomiser:
    input:
        exomiser_jar = expand("{reference_dir}/exomiser/exomiser-cli-{version}/exomiser-cli-{version}.jar", version=config["exomiser"]["jar"]["version"],reference_dir=reference_dir),
        yml = f"{workpath}/phenotype/{{proband}}_exomiser.yaml",
        config = expand("{reference_dir}/exomiser/exomiser-cli-{version}/application.properties",version=config["exomiser"]["jar"]["version"],reference_dir=reference_dir),
    output:
        html = f"{workpath}/exomiser/{{proband}}/{{proband}}-exomiser.html",
    params:
        outformats = "HTML,JSON,TSV_GENE,TSV_VARIANT,VCF",
    log:
        "logs/exomiser/run_exomiser/{proband}-run_exomiser.log"
    shell:
        """
        java -jar {input.exomiser_jar} \
            --analysis {input.yml} \
            --spring.config.location={input.config} \
            --output-format {params.outformats} \
            > {log} 2>&1
        """

