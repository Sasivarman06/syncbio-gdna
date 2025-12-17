if trio_analysis:
    rule gatk_haplotypecaller_trio:
        input:
            reference = f"{reference_dir}/genome/hg38.fa",
            reference_idx = f"{reference_dir}/genome/hg38.fa.fai",
            sequence_dict = f"{reference_dir}/genome/hg38.dict",
            child_bam = lambda w: get_bam(w.proband),
            child_bam_idx = lambda w: get_bam(w.proband) + ".bai",
            father_bam = lambda w: get_bam(trio_map[w.proband]["father"]),
            father_bam_idx = lambda w: get_bam(trio_map[w.proband]["father"]) + ".bai",
            mother_bam = lambda w: get_bam(trio_map[w.proband]["mother"]),
            mother_bam_idx = lambda w: get_bam(trio_map[w.proband]["mother"]) + ".bai"
        output:
            gvcf_child = f"{workpath}/vcf/{{proband}}_gatk.g.vcf.gz",
            gvcf_child_idx = f"{workpath}/vcf/{{proband}}_gatk.g.vcf.gz.tbi",
            gvcf_father = f"{workpath}/vcf/{{proband}}_father_gatk.g.vcf.gz",
            gvcf_father_idx = f"{workpath}/vcf/{{proband}}_father_gatk.g.vcf.gz.tbi",
            gvcf_mother = f"{workpath}/vcf/{{proband}}_mother_gatk.g.vcf.gz",
            gvcf_mother_idx = f"{workpath}/vcf/{{proband}}_mother_gatk.g.vcf.gz.tbi"
        log:
            "logs/variant-calling/gatk_haplotypecaller_trio/{{proband}}.log"
        conda:
            "../envs/gatk.yaml"
        shell:
            """
            gatk --java-options "-Xmx6g" HaplotypeCaller \
              -R {input.reference} -I {input.child_bam} -O {output.gvcf_child} -ERC GVCF > {log} 2>&1
            gatk --java-options "-Xmx6g" HaplotypeCaller \
              -R {input.reference} -I {input.father_bam} -O {output.gvcf_father} -ERC GVCF >> {log} 2>&1
            gatk --java-options "-Xmx6g" HaplotypeCaller \
              -R {input.reference} -I {input.mother_bam} -O {output.gvcf_mother} -ERC GVCF >> {log} 2>&1
            """

    rule genomicsdb_import:
        input:
            gvcf_child = f"{workpath}/vcf/{{proband}}_gatk.g.vcf.gz",
            gvcf_child_idx = f"{workpath}/vcf/{{proband}}_gatk.g.vcf.gz.tbi",
            gvcf_father = f"{workpath}/vcf/{{proband}}_father_gatk.g.vcf.gz",
            gvcf_father_idx = f"{workpath}/vcf/{{proband}}_father_gatk.g.vcf.gz.tbi",
            gvcf_mother = f"{workpath}/vcf/{{proband}}_mother_gatk.g.vcf.gz",
            gvcf_mother_idx = f"{workpath}/vcf/{{proband}}_mother_gatk.g.vcf.gz.tbi",
            intervals=f"{reference_dir}/genome/intervals.list"
        output:
            db_dir=directory(f"{workpath}/vcf/{{proband}}_genomicsdb")
        params:
            intervals=config["intervals"]
        log:
            "logs/variant-calling/gatk_genomicsdb_import/{{proband}}.log"
        conda: 
            "../envs/gatk.yaml"
        threads: 4
        shell:
            """
            gatk --java-options "-Xmx8G" GenomicsDBImport \
            -V {input.gvcf_mother} \
            -V {input.gvcf_father} \
            -V {input.gvcf_child} \
            --genomicsdb-workspace-path {output.db_dir} \
            --batch-size 50 \
            --reader-threads {threads} \
            --intervals {params.intervals} \
            --tmp-dir . > {log} 2>&1
            """

    rule genotype_gvcfs_trio:
        input:
            db_dir=f"{workpath}/vcf/{{proband}}_genomicsdb",
            reference=f"{reference_dir}/genome/hg38.fa",
            reference_idx=f"{reference_dir}/genome/hg38.fa.fai",
            reference_dict=f"{reference_dir}/genome/hg38.dict"
        output:
            vcf=f"{workpath}/vcf/{{proband}}_trio_gatk.vcf.gz",
            vcf_idx=f"{workpath}/vcf/{{proband}}_trio_gatk.vcf.gz.tbi"
        log:
            "logs/variant-calling/genotype_gvcfs_trio/{{proband}}.log"
        threads: 4
        conda:
            "../envs/gatk.yaml"
        shell:
            """
            gatk --java-options "-Xmx6G" GenotypeGVCFs \
                -R {input.reference} \
                -V gendb://{input.db_dir} \
                -O {output.vcf} > {log} 2>&1
            """

    rule deeptrio:
        input:
            reference = f"{reference_dir}/genome/hg38.fa",
            reference_idx = f"{reference_dir}/genome/hg38.fa.fai",
            reference_dict = f"{reference_dir}/genome/hg38.dict",
            child_bam = lambda w: get_bam(w.proband),
            child_bam_idx = lambda w: get_bam(w.proband) + ".bai",
            father_bam = lambda w: get_bam(trio_map[w.proband]["father"]),
            father_bam_idx = lambda w: get_bam(trio_map[w.proband]["father"]) + ".bai",
            mother_bam = lambda w: get_bam(trio_map[w.proband]["mother"]),
            mother_bam_idx = lambda w: get_bam(trio_map[w.proband]["mother"]) + ".bai"
        output:
            vcf_child   = f"{workpath}/vcf/{{proband}}_deeptrio.vcf.gz",
            vcf_child_idx = f"{workpath}/vcf/{{proband}}_deeptrio.vcf.gz.tbi",
            vcf_father  = f"{workpath}/vcf/{{proband}}_father_deeptrio.vcf.gz",
            vcf_father_idx = f"{workpath}/vcf/{{proband}}_father_deeptrio.vcf.gz.tbi",
            vcf_mother  = f"{workpath}/vcf/{{proband}}_mother_deeptrio.vcf.gz",
            vcf_mother_idx = f"{workpath}/vcf/{{proband}}_mother_deeptrio.vcf.gz.tbi",
            gvcf_child  = f"{workpath}/vcf/{{proband}}_deeptrio.g.vcf.gz",
            gvcf_child_idx = f"{workpath}/vcf/{{proband}}_deeptrio.g.vcf.gz.tbi",
            gvcf_father = f"{workpath}/vcf/{{proband}}_father_deeptrio.g.vcf.gz",
            gvcf_father_idx = f"{workpath}/vcf/{{proband}}_father_deeptrio.g.vcf.gz.tbi",
            gvcf_mother = f"{workpath}/vcf/{{proband}}_mother_deeptrio.g.vcf.gz",
            gvcf_mother_idx = f"{workpath}/vcf/{{proband}}_mother_deeptrio.g.vcf.gz.tbi"
        params:
            model_type = "WES",
            father = lambda w: trio_map[w.proband]["father"],
            mother = lambda w: trio_map[w.proband]["mother"]
        log:    
            "logs/variant-calling/deeptrio/{{proband}}.log"
        threads: 8
        singularity:
            "docker://google/deepvariant:deeptrio-1.8.0"
        shell:
            """
            /opt/deepvariant/bin/deeptrio/run_deeptrio \
            --model_type {params.model_type} \
            --ref {input.reference} \
            --reads_child {input.child_bam} \
            --reads_parent1 {input.father_bam} \
            --reads_parent2 {input.mother_bam} \
            --output_vcf_child {output.vcf_child} \
            --output_vcf_parent1 {output.vcf_father} \
            --output_vcf_parent2 {output.vcf_mother} \
            --sample_name_child '{wildcards.proband}' \
            --sample_name_parent1 '{params.father}' \
            --sample_name_parent2 '{params.mother}' \
            --num_shards {threads} \
            --output_gvcf_child {output.gvcf_child} \
            --output_gvcf_parent1 {output.gvcf_father} \
            --output_gvcf_parent2 {output.gvcf_mother} > {log} 2>&1
            """

    rule merge_deeptrio:
        input:
            gvcf_child = f"{workpath}/vcf/{{proband}}_deeptrio.g.vcf.gz",
            gvcf_child_idx = f"{workpath}/vcf/{{proband}}_deeptrio.g.vcf.gz.tbi",
            gvcf_father = f"{workpath}/vcf/{{proband}}_father_deeptrio.g.vcf.gz",
            gvcf_father_idx = f"{workpath}/vcf/{{proband}}_father_deeptrio.g.vcf.gz.tbi",
            gvcf_mother = f"{workpath}/vcf/{{proband}}_mother_deeptrio.g.vcf.gz",   
            gvcf_mother_idx = f"{workpath}/vcf/{{proband}}_mother_deeptrio.g.vcf.gz.tbi"
        output:
            vcf=f"{workpath}/vcf/{{proband}}_trio_deeptrio.vcf.gz",
            vcf_idx=f"{workpath}/vcf/{{proband}}_trio_deeptrio.vcf.gz.tbi",
            db = temp(directory(f"{workpath}/vcf/{{proband}}_GLnexus.db"))       
        log:
            "logs/variant-calling/merge_deeptrio/{{proband}}.log"
        threads: 4
        conda:
            "../envs/glnexus.yaml"
        shell:
            """
            ( glnexus_cli \
                --config DeepVariant_unfiltered \
                --dir {output.db} \
                {input.gvcf_child} \
                {input.gvcf_father} \
                {input.gvcf_mother} \
                | bcftools view -Oz -o {output.vcf}
            ) > {log} 2>&1
            tabix -p vcf {output.vcf} >> {log} 2>&1
            """
    
    rule isec_gatk_deeptrio_trio:
        input:
            gatk_trio = f"{workpath}/vcf/{{proband}}_trio_gatk.vcf.gz",
            gatk_trio_idx = f"{workpath}/vcf/{{proband}}_trio_gatk.vcf.gz.tbi",
            deeptrio_trio = f"{workpath}/vcf/{{proband}}_trio_deeptrio.vcf.gz",
            deeptrio_trio_idx = f"{workpath}/vcf/{{proband}}_trio_deeptrio.vcf.gz.tbi"
        output:
            tempdir = directory(f"{workpath}/vcf/{{proband}}_trio_isec"),
            deeptrio_unique_vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0000.vcf",
            gatk_unique_vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0001.vcf",
            common_vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0002.vcf",
            other_vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0003.vcf"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/isec_gatk_deeptrio_trio/{{proband}}.log"
        shell:
            """
            bcftools isec -p {output.tempdir} -c both  {input.gatk_trio} {input.deeptrio_trio} > {log} 2>&1
            """

    rule extract_deeptrio_calls_trio:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0000.vcf"
        output:
            anno = f"{workpath}/vcf/{{proband}}_trio_dv_anno.txt",
            anno_gz = f"{workpath}/vcf/{{proband}}_trio_dv_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_trio_dv_anno.txt.gz.tbi"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/extract_deeptrio_calls_trio/{{proband}}.log"
        shell:
            """
            bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input.vcf} | awk '{{print $0 "\\tDeepTrio"}}' > {output.anno} 2> {log}
            bgzip -c {output.anno} > {output.anno_gz}  2>> {log}
            tabix -s1 -b2 -e2 {output.anno_gz} 2>> {log}
            """

    rule extract_gatk_calls_trio:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0001.vcf"
        output:
            anno = f"{workpath}/vcf/{{proband}}_trio_gatk_anno.txt",
            anno_gz = f"{workpath}/vcf/{{proband}}_trio_gatk_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_trio_gatk_anno.txt.gz.tbi"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/extract_gatk_calls_trio/{{proband}}.log"
        shell:
            """
            bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input.vcf} | awk '{{print $0 "\\tGATK"}}' > {output.anno} 2> {log}
            bgzip -c {output.anno} > {output.anno_gz}  2>> {log}
            tabix -s1 -b2 -e2 {output.anno_gz}  2>> {log}
            """

    rule extract_common_calls_trio:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0002.vcf"
        output:
            anno = f"{workpath}/vcf/{{proband}}_trio_common_anno.txt",
            anno_gz = f"{workpath}/vcf/{{proband}}_trio_common_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_trio_common_anno.txt.gz.tbi"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/extract_common_calls_trio/{{proband}}.log"
        shell:
            """
            bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input.vcf} | awk '{{print $0 "\\tDeepTrio,GATK"}}' > {output.anno} 2> {log}
            bgzip -c {output.anno} > {output.anno_gz} 2>> {log}
            tabix -s1 -b2 -e2 {output.anno_gz} 2>> {log}
            """

    rule annotate_deeptrio_trio:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0000.vcf",
            anno_gz = f"{workpath}/vcf/{{proband}}_trio_dv_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_trio_dv_anno.txt.gz.tbi"
        output:
            vcf_annotated = f"{workpath}/vcf/{{proband}}_trio_deeptrio.unique.vcf.gz",
            vcf_annotated_idx = f"{workpath}/vcf/{{proband}}_trio_deeptrio.unique.vcf.gz.tbi",
            header = f"{workpath}/vcf/{{proband}}_trio_dv_caller_header.txt"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/annotate_deeptrio_trio/{{proband}}.log"
        shell:
            """
            echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Variant caller source">' > {output.header} 2> {log}
            bcftools annotate -h {output.header} -a {input.anno_gz} -c CHROM,POS,REF,ALT,CALLER -Oz -o {output.vcf_annotated} {input.vcf} >> {log} 2>&1
            tabix -p vcf {output.vcf_annotated} >> {log} 2>&1
            """

    rule annotate_gatk_trio:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0001.vcf",
            anno_gz = f"{workpath}/vcf/{{proband}}_trio_gatk_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_trio_gatk_anno.txt.gz.tbi"
        output:
            vcf_annotated = f"{workpath}/vcf/{{proband}}_trio_gatk.unique.vcf.gz",
            vcf_annotated_idx = f"{workpath}/vcf/{{proband}}_trio_gatk.unique.vcf.gz.tbi",
            header = f"{workpath}/vcf/{{proband}}_trio_gatk_caller_header.txt"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/annotate_gatk_trio/{{proband}}.log"
        shell:
            """
            echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Variant caller source">' > {output.header} 2> {log}
            bcftools annotate -h {output.header} -a {input.anno_gz} -c CHROM,POS,REF,ALT,CALLER -Oz -o {output.vcf_annotated} {input.vcf} >> {log} 2>&1
            tabix -p vcf {output.vcf_annotated} >> {log} 2>&1
            """

    rule annotate_common_trio:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_trio_isec/0002.vcf",
            anno_gz = f"{workpath}/vcf/{{proband}}_trio_common_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_trio_common_anno.txt.gz.tbi"
        output:
            vcf_annotated = f"{workpath}/vcf/{{proband}}_trio.common.vcf.gz",
            vcf_annotated_idx = f"{workpath}/vcf/{{proband}}_trio.common.vcf.gz.tbi",
            header = f"{workpath}/vcf/{{proband}}_trio_common_caller_header.txt"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/annotate_common_trio/{{proband}}.log"
        shell:
            """
            echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Variant caller source">' > {output.header} 2> {log}
            bcftools annotate -h {output.header} -a {input.anno_gz} -c CHROM,POS,REF,ALT,CALLER -Oz -o {output.vcf_annotated} {input.vcf} >> {log} 2>&1
            tabix -p vcf {output.vcf_annotated} >> {log} 2>&1
            """

    rule concat_annotated_vcfs_trio:
        input:
            deeptrio_unique = f"{workpath}/vcf/{{proband}}_trio_deeptrio.unique.vcf.gz",
            deeptrio_unique_idx = f"{workpath}/vcf/{{proband}}_trio_deeptrio.unique.vcf.gz.tbi",
            gatk_unique = f"{workpath}/vcf/{{proband}}_trio_gatk.unique.vcf.gz",
            gatk_unique_idx = f"{workpath}/vcf/{{proband}}_trio_gatk.unique.vcf.gz.tbi",
            common = f"{workpath}/vcf/{{proband}}_trio.common.vcf.gz",
            common_idx = f"{workpath}/vcf/{{proband}}_trio.common.vcf.gz.tbi"
        output:
            merged = f"{workpath}/vcf/{{proband}}_trio_merged.vcf.gz",
            merged_idx = f"{workpath}/vcf/{{proband}}_trio_merged.vcf.gz.tbi"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/concat_annotated_vcfs_trio/{{proband}}.log"
        shell:
            """
            bcftools concat --allow-overlaps -Oz -o {output.merged} \
                {input.deeptrio_unique} {input.gatk_unique} {input.common} > {log} 2>&1
            tabix -p vcf {output.merged} >> {log} 2>&1
            """


if solo_analysis:
    rule gatk_haplotypecaller_solo:
        input:
            reference=f"{reference_dir}/genome/hg38.fa",
            reference_idx=f"{reference_dir}/genome/hg38.fa.fai",
            sequence_dict=f"{reference_dir}/genome/hg38.dict",
            proband_bam = lambda w: get_bam(w.proband),
            proband_bam_idx = lambda w: get_bam(w.proband) + ".bai"
        output:
            vcf=f"{workpath}/vcf/{{proband}}_gatk.vcf.gz",
            vcf_idx=f"{workpath}/vcf/{{proband}}_gatk.vcf.gz.tbi"
        log:
            "logs/variant-calling/gatk_haplotypecaller_solo/{{proband}}.log"
        conda:
            "../envs/gatk.yaml"
        shell:
            """
            gatk --java-options "-Xmx8g" HaplotypeCaller \
              -R {input.reference} \
              -I {input.proband_bam} \
              -O {output.vcf} > {log} 2>&1
            """


    rule deepvariant:
        input:
            reference=f"{reference_dir}/genome/hg38.fa",
            reference_idx=f"{reference_dir}/genome/hg38.fa.fai",
            reference_dict=f"{reference_dir}/genome/hg38.dict",
            proband_bam= lambda w: get_bam(w.proband),
            proband_bam_idx= lambda w: get_bam(w.proband) + ".bai"
        output:
            vcf=f"{workpath}/vcf/{{proband}}_deepvariant.vcf.gz",
            vcf_idx=f"{workpath}/vcf/{{proband}}_deepvariant.vcf.gz.tbi",
            gvcf=f"{workpath}/vcf/{{proband}}_proband_deepvariant.g.vcf.gz",
            gvcf_idx=f"{workpath}/vcf/{{proband}}_proband_deepvariant.g.vcf.gz.tbi",
            html_report=f"{workpath}/vcf/{{proband}}_deepvariant.visual_report.html"
        params:
            model_type="WES",
        threads: 4
        singularity:
            "docker://google/deepvariant:1.9.0"
        log:
            "logs/variant-calling/deepvariant/{{proband}}.log"
        shell:
            """
            /opt/deepvariant/bin/run_deepvariant \
            --model_type={params.model_type} \
            --ref={input.reference} \
            --reads={input.proband_bam} \
            --sample_name={wildcards.proband} \
            --output_vcf={output.vcf} \
            --output_gvcf={output.gvcf} \
            --vcf_stats_report \
            --num_shards={threads} > {log} 2>&1
            """


    rule isec_gatk_deepvariant:
        input:
            gatk = f"{workpath}/vcf/{{proband}}_gatk.vcf.gz",
            gatk_idx = f"{workpath}/vcf/{{proband}}_gatk.vcf.gz.tbi",
            deepvariant = f"{workpath}/vcf/{{proband}}_deepvariant.vcf.gz",
            deepvariant_idx = f"{workpath}/vcf/{{proband}}_deepvariant.vcf.gz.tbi"
        output:
            tempdir = directory(f"{workpath}/vcf/{{proband}}_isec"),
            gatk_unique_vcf = f"{workpath}/vcf/{{proband}}_isec/0000.vcf",
            deepvariant_unique_vcf = f"{workpath}/vcf/{{proband}}_isec/0001.vcf",
            common_vcf = f"{workpath}/vcf/{{proband}}_isec/0002.vcf",
            other_vcf = f"{workpath}/vcf/{{proband}}_isec/0003.vcf"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/isec_gatk_deepvariant/{{proband}}.log"
        shell:
            """
            bcftools isec -p {output.tempdir} -c both {input.gatk} {input.deepvariant} > {log} 2>&1
            """

    rule extract_deepvariant_calls:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_isec/0000.vcf"
        output:
            anno = f"{workpath}/vcf/{{proband}}_dv_anno.txt",
            anno_gz = f"{workpath}/vcf/{{proband}}_dv_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_dv_anno.txt.gz.tbi"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/extract_deepvariant_calls/{{proband}}.log"
        shell:
            """
            bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input.vcf} | awk '{{print $0 "\\tDeepVariant"}}' > {output.anno} 2> {log}
            bgzip -c {output.anno} > {output.anno_gz} 2>> {log} 
            tabix -s1 -b2 -e2 {output.anno_gz} 2>> {log}
            """

    rule extract_gatk_calls:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_isec/0001.vcf"
        output:
            anno = f"{workpath}/vcf/{{proband}}_gatk_anno.txt",
            anno_gz = f"{workpath}/vcf/{{proband}}_gatk_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_gatk_anno.txt.gz.tbi"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/extract_gatk_calls/{{proband}}.log"
        shell:
            """
            bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input.vcf} | awk '{{print $0 "\\tGATK"}}' > {output.anno} 2> {log}
            bgzip -c {output.anno} > {output.anno_gz} 2>> {log}
            tabix -s1 -b2 -e2 {output.anno_gz} 2>> {log}
            """

    rule extract_common_calls:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_isec/0002.vcf"
        output:
            anno = f"{workpath}/vcf/{{proband}}_common_anno.txt",
            anno_gz = f"{workpath}/vcf/{{proband}}_common_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_common_anno.txt.gz.tbi"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/extract_common_calls/{{proband}}.log"
        shell:
            """
            bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input.vcf} | awk '{{print $0 "\\tDeepVariant,GATK"}}' > {output.anno} 2> {log}
            bgzip -c {output.anno} > {output.anno_gz} 2>> {log}
            tabix -s1 -b2 -e2 {output.anno_gz} 2>> {log} 
            """

    rule annotate_deepvariant:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_isec/0000.vcf",
            anno_gz = f"{workpath}/vcf/{{proband}}_dv_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_dv_anno.txt.gz.tbi"
        output:
            vcf_annotated = f"{workpath}/vcf/{{proband}}_deepvariant.unique.vcf.gz",
            vcf_annotated_idx = f"{workpath}/vcf/{{proband}}_deepvariant.unique.vcf.gz.tbi",
            header = f"{workpath}/vcf/{{proband}}_dv_caller_header.txt"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/annotate_deepvariant/{{proband}}.log"
        shell:
            """
            echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Variant caller source">' > {output.header} 2> {log}
            bcftools annotate -h {output.header} -a {input.anno_gz} -c CHROM,POS,REF,ALT,CALLER -Oz -o {output.vcf_annotated} {input.vcf} >> {log} 2>&1
            tabix -p vcf {output.vcf_annotated} >> {log} 2>&1
            """

    rule annotate_gatk:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_isec/0001.vcf",
            anno_gz = f"{workpath}/vcf/{{proband}}_gatk_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_gatk_anno.txt.gz.tbi"
        output:
            vcf_annotated = f"{workpath}/vcf/{{proband}}_gatk.unique.vcf.gz",
            vcf_annotated_idx = f"{workpath}/vcf/{{proband}}_gatk.unique.vcf.gz.tbi",
            header = f"{workpath}/vcf/{{proband}}_gatk_caller_header.txt"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/annotate_gatk/{{proband}}.log"
        shell:
            """
            echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Variant caller source">' > {output.header} 2> {log}
            bcftools annotate -h {output.header} -a {input.anno_gz} -c CHROM,POS,REF,ALT,CALLER -Oz -o {output.vcf_annotated} {input.vcf} >> {log} 2>&1
            tabix -p vcf {output.vcf_annotated} >> {log} 2>&1
            """

    rule annotate_common:
        input:
            vcf = f"{workpath}/vcf/{{proband}}_isec/0002.vcf",
            anno_gz = f"{workpath}/vcf/{{proband}}_common_anno.txt.gz",
            anno_tbi = f"{workpath}/vcf/{{proband}}_common_anno.txt.gz.tbi"
        output:
            vcf_annotated = f"{workpath}/vcf/{{proband}}.common.vcf.gz",
            vcf_annotated_idx = f"{workpath}/vcf/{{proband}}.common.vcf.gz.tbi",
            header = f"{workpath}/vcf/{{proband}}_common_caller_header.txt"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/annotate_common/{{proband}}.log"
        shell:
            """
            echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Variant caller source">' > {output.header} 2> {log}
            bcftools annotate -h {output.header} -a {input.anno_gz} -c CHROM,POS,REF,ALT,CALLER -Oz -o {output.vcf_annotated} {input.vcf} >> {log} 2>&1
            tabix -p vcf {output.vcf_annotated} >> {log} 2>&1
            """

    rule concat_annotated_vcfs:
        input:
            dv_unique = f"{workpath}/vcf/{{proband}}_deepvariant.unique.vcf.gz",
            dv_unique_idx = f"{workpath}/vcf/{{proband}}_deepvariant.unique.vcf.gz.tbi",
            gatk_unique = f"{workpath}/vcf/{{proband}}_gatk.unique.vcf.gz",
            gatk_unique_idx = f"{workpath}/vcf/{{proband}}_gatk.unique.vcf.gz.tbi",
            common = f"{workpath}/vcf/{{proband}}.common.vcf.gz",
            common_idx = f"{workpath}/vcf/{{proband}}.common.vcf.gz.tbi"
        output:
            merged = f"{workpath}/vcf/{{proband}}_merged.vcf.gz",
            merged_idx = f"{workpath}/vcf/{{proband}}_merged.vcf.gz.tbi"
        conda:
            "../envs/bcf.yaml"
        log:
            "logs/variant-calling/concat_annotated_vcfs/{{proband}}.log"
        shell:
            """
            bcftools concat --allow-overlaps -Oz -o {output.merged} \
                {input.dv_unique} {input.gatk_unique} {input.common} > {log} 2>&1
            tabix -p vcf {output.merged} >> {log} 2>&1
            """