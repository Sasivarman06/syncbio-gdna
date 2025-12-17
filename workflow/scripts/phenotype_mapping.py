# iderare_phenomizing_snake.py
from iderare_pheno.converter import batchconvert
from iderare_pheno.simrec import hpo2omim_similarity, omim_recommendation, hpo2name, omim2name
from iderare_pheno.utils import linkage_dendrogram, list2tsv
import yaml
import os

workpath = os.getcwd()


# -----------------------------
# Default parameters function
# -----------------------------
def iderare_pheno_param(threshold=0.4, differential=10, recommendation=20, add_yml=False):
    return threshold, differential, recommendation, add_yml

# -----------------------------
# Exomiser YAML generation function
# -----------------------------
def generate_exomiser_yml(hpo_sets, vcf_file, ped_file, proband_id, output_file, is_trio=False):
    exomiser_template = {
        'analysis': {
            'genomeAssembly': 'hg38',
            'vcf': vcf_file,
            'ped': ped_file,
            'proband': proband_id,
            'hpoIds': hpo_sets,
            'inheritanceModes': {
                'AUTOSOMAL_DOMINANT': 0.1,
                'AUTOSOMAL_RECESSIVE_HOM_ALT': 0.1,
                'AUTOSOMAL_RECESSIVE_COMP_HET': 2.0,
                'X_DOMINANT': 0.1,
                'X_RECESSIVE_HOM_ALT': 0.1,
                'X_RECESSIVE_COMP_HET': 2.0,
                'MITOCHONDRIAL': 0.2
            },
            'analysisMode': 'PASS_ONLY',
            'frequencySources': [
                'THOUSAND_GENOMES', 'TOPMED', 'UK10K',
                'ESP_AFRICAN_AMERICAN','ESP_EUROPEAN_AMERICAN','ESP_ALL',
                'EXAC_AFRICAN_INC_AFRICAN_AMERICAN','EXAC_AMERICAN','EXAC_SOUTH_ASIAN','EXAC_EAST_ASIAN',
                'EXAC_FINNISH','EXAC_NON_FINNISH_EUROPEAN','EXAC_OTHER',
                'GNOMAD_E_AFR','GNOMAD_E_AMR','GNOMAD_E_EAS','GNOMAD_E_FIN','GNOMAD_E_NFE','GNOMAD_E_OTH','GNOMAD_E_SAS',
                'GNOMAD_G_AFR','GNOMAD_G_AMR','GNOMAD_G_EAS','GNOMAD_G_FIN','GNOMAD_G_NFE','GNOMAD_G_OTH','GNOMAD_G_SAS'
            ],
            'pathogenicitySources': ['REVEL', 'MVP'],
            'steps': [
                {'failedVariantFilter': {}},
                {'variantEffectFilter': {'remove': [
                    'FIVE_PRIME_UTR_EXON_VARIANT','FIVE_PRIME_UTR_INTRON_VARIANT',
                    'THREE_PRIME_UTR_EXON_VARIANT','THREE_PRIME_UTR_INTRON_VARIANT',
                    'NON_CODING_TRANSCRIPT_EXON_VARIANT','UPSTREAM_GENE_VARIANT',
                    'INTERGENIC_VARIANT','REGULATORY_REGION_VARIANT',
                    'CODING_TRANSCRIPT_INTRON_VARIANT','NON_CODING_TRANSCRIPT_INTRON_VARIANT',
                    'DOWNSTREAM_GENE_VARIANT'
                ]}},
                {'frequencyFilter': {'maxFrequency': 2.0}},
                {'pathogenicityFilter': {'keepNonPathogenic': True}},
                {'inheritanceFilter': {}},
                {'omimPrioritiser': {}},
                {'hiPhivePrioritiser': {}},
                {'phenixPrioritiser': {}}
            ]
        },
        'outputOptions': {
            'outputContributingVariantsOnly': True,
            'numGenes': 50,
            'outputDirectory': f'{workpath}/exomiser/{proband_id}',
            'outputFileName': f'{proband_id}-exomiser',
            'outputFormats': ['HTML', 'JSON', 'TSV_GENE', 'TSV_VARIANT', 'VCF']
        }
    }

    if is_trio and ped_file:
        exomiser_template['analysis']['ped'] = ped_file

    with open(output_file, 'w') as f:
        yaml.dump(exomiser_template, f, sort_keys=False)
    print(f"Exomiser YAML written to: {output_file}")


# -----------------------------
# Main Snakemake execution
# -----------------------------
if "snakemake" in globals():
    clinical_data = snakemake.input.clinical
    vcf_file = snakemake.input.vcf
    if isinstance(vcf_file, str):
        vcf_paths = [vcf_file]
    elif isinstance(vcf_file, (list, tuple)):
        vcf_paths = list(vcf_file)
    else:
        # NamedList or unexpected type
        try:
            vcf_paths = list(vcf_file)
        except Exception:
            vcf_paths = [str(vcf_file)]

    ped_file = snakemake.input.ped
    proband_id = snakemake.params.proband
    workpath = snakemake.params.workpath
    trio_analysis = snakemake.params.trio_analysis
    solo_analysis = snakemake.params.solo_analysis
    trio_probands = snakemake.params.trio_probands


    threshold = getattr(snakemake.params, "threshold", 0.4)
    differential = getattr(snakemake.params, "differential", 10)
    recommendation = getattr(snakemake.params, "recommendation", 20)

    # Read clinical data
    with open(clinical_data, 'r') as file:
        clinical_data_list = file.read().splitlines()

    # Phenotype processing
    hpo_sets, diagnosis_sets = batchconvert(clinical_data_list)
    s_sim, [lnk_all, sr_dis_name, sr_dis_id], [lnk_thr, sr_dis_name_thr, sr_dis_id_thr] = hpo2omim_similarity(
        diagnosis_sets, hpo_sets, threshold=threshold, differential=differential
    )

    # -----------------------------
    # Use snakemake.output for all files
    # -----------------------------
    # Expected snakemake.output should be a dictionary with all output paths
    # Example in your Snakefile:
    # output:
    #     yml = "results/{proband}_exomiser.yml",
    #     hpo_list = "results/{proband}_transformed_hpo_set.txt",
    #     hpo_tsv = "results/{proband}_transformed_hpo_set.tsv",
    #     omim_tsv = "results/{proband}_transformed_omim_set.tsv",
    #     ddx_sim = "results/{proband}_differential_diagnosis_similarity.tsv",
    #     gene_sim = "results/{proband}_recommended_gene_similarity.tsv",
    #     disease_sim = "results/{proband}_recommended_disease_similarity.tsv",
    #     dendro_all = "results/{proband}_linkage_all_ddx.png",
    #     dendro_filt = "results/{proband}_linkage_filt_ddx.png",
    #     dendro_gene = "results/{proband}_linkage_gene.png",
    #     dendro_disease = "results/{proband}_linkage_disease.png"

    output = snakemake.output

    # Dendrograms
    linkage_dendrogram(lnk_all, sr_dis_name, title='Linkage of All DDx', threshold=threshold,
                       path_to_save=output.dendro_all)
    linkage_dendrogram(lnk_thr, sr_dis_name_thr, title='Linkage of Filtered DDx', threshold=threshold,
                       path_to_save=output.dendro_filt)
    rg_s_sim, [rg_lnk_all, rg_sr_dis_name, rg_sr_dis_id], [rg_lnk_thr, rg_sr_dis_name_thr, rg_sr_dis_id_thr] = omim_recommendation(
        hpo_sets, type='gene', threshold=threshold, recommendation=recommendation
    )
    linkage_dendrogram(rg_lnk_thr, rg_sr_dis_name_thr, title='Linkage of Causative Gene', threshold=threshold,
                       path_to_save=output.dendro_gene)
    rd_s_sim, [rd_lnk_all, rd_sr_dis_name, rd_sr_dis_id], [rd_lnk_thr, rd_sr_dis_name_thr, rd_sr_dis_id_thr] = omim_recommendation(
        hpo_sets, type='disease', threshold=threshold, recommendation=recommendation
    )
    linkage_dendrogram(rd_lnk_thr, rd_sr_dis_name_thr, title='Linkage of Causative Disease', threshold=threshold,
                       path_to_save=output.dendro_disease)

    # TSV files
    hpo_name = hpo2name(hpo_sets)
    list2tsv(hpo_sets, hpo_name, filename=output.hpo_tsv.replace(".tsv", ""))
    with open(output.hpo_list, 'w') as f:
        f.write(str(hpo_sets))
    omim_name = omim2name(diagnosis_sets)
    list2tsv(diagnosis_sets, omim_name, filename=output.omim_tsv.replace(".tsv", ""))
    list2tsv(sr_dis_id, sr_dis_name, s_sim, filename=output.ddx_sim.replace(".tsv", ""))
    list2tsv(rg_sr_dis_id, rg_sr_dis_name, rg_s_sim, filename=output.gene_sim.replace(".tsv", ""))
    list2tsv(rd_sr_dis_id, rd_sr_dis_name, rd_s_sim, filename=output.disease_sim.replace(".tsv", ""))

    # Check if proband is trio
    is_trio = proband_id in trio_probands

    # Exomiser YAML
    generate_exomiser_yml(
        hpo_sets=hpo_sets,
        vcf_file=vcf_paths,
        ped_file=snakemake.input.ped if is_trio else None,
        proband_id=proband_id,
        output_file=output.yml,
        is_trio=is_trio
    )
