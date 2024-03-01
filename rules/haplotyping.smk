rule copy_snps:
    input:
        snp_df['file_path'].to_list()
    output:
        snp_df['out_path'].to_list(),
    run:
        from shutil import copyfile
        for i, fpath in enumerate(input):
    
            outPath = output[i]
            copyfile(fpath, outPath)


rule index_vcf:
    input:
        OUTPUT + 'vcf/{sid}.vcf.gz'
    output: 
        OUTPUT + 'vcf/{sid}.vcf.gz.tbi'
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    shell:
        "tabix -fp vcf {input}"


rule unzip_vcf:
    input:
        OUTPUT + 'vcf/{sid}.vcf.gz'
    output:
        OUTPUT + 'vcf/{sid}.vcf'
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    shell:
        """gzip -dk {input}"""


rule to_variants_table:
    input:
        vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
    output:
        OUTPUT + 'vcf/{sid}.snps.tsv',
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    shell:
        """bcftools query --print-header -f \
        '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' {input.vcf} > {output}"""


rule phase_vcf_longread:
    input:
        bam=expand(OUTPUT + 'merged_bam/{cid}.GRCm39.bam', cid=cell_ids,),
        bamindex=expand(OUTPUT + 'merged_bam/{cid}.GRCm39.bam.bai', cid=cell_ids,),
        ref=OUTPUT + 'references/GRCm39.fa',
        vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
        vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
    output:
        vcf=OUTPUT + "vcf/{sid}.phased.vcf",
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    threads:
        config['threads'] 
    params:
        longread=config['longread_path'],
        prefix=lambda wildcards: OUTPUT + "vcf/" + wildcards.sid + ".phased"
    shell:
        """{params.longread} phase -s {input.vcf} -b {input.bam} -r {input.ref} -t {threads} -o {params.prefix} --ont"""


rule longread_haplotag:
    input:
        bam=OUTPUT + 'merged_bam/{cid}.GRCm39.bam',
        bamindex=OUTPUT + 'merged_bam/{cid}.GRCm39.bam.bai',
        ref=OUTPUT + 'references/GRCm39.fa',
        vcf=OUTPUT + 'vcf/{sid}.phased.vcf',
    output:
        bam=OUTPUT + "longread/{cid}.{sid}.phased.bam",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    threads:
        config['threads'] // 4
    params:
        longread=config['longread_path'],
        prefix=lambda wildcards: OUTPUT + "longread/" + wildcards.cid + "." + wildcards.sid + ".phased"
    shell:
        """{params.longread} haplotag -r {input.ref} -s {input.vcf} -b {input.bam} -t {threads} -o {params.prefix} --log"""


rule gatk_VariantsToTable:
    input:
        vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
        vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
    output:
        OUTPUT + 'vcf/{sid}.gatk.table.tsv'
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    shell:
        """gatk VariantsToTable \
         -V {input.vcf} \
         -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -GF PL \
         -O {output}"""


rule whatshap_stats_CAST_EiJ:
    input:
        vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
        vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
    output:
        OUTPUT + "reports/whatshap/{sid}.CAST_EiJ.stats.txt"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    shell:
        """whatshap stats \
        --sample CAST_EiJ \
        {input.vcf} > {output}"""


rule whatshap_stats_129S1_SvImJ:
    input:
        vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
        vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
    output:
        OUTPUT + "reports/whatshap/{sid}.129S1_SvImJ.stats.txt"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    shell:
        """whatshap stats \
        --sample 129S1_SvImJ \
        {input.vcf} > {output}"""


rule whatshap_stats_phased_CAST_EiJ:
    input:
        vcf=OUTPUT + 'vcf/{sid}.phased.vcf',
    output:
        OUTPUT + "reports/whatshap/{sid}.CAST_EiJ.phased.stats.txt"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    shell:
        """whatshap stats \
        --sample CAST_EiJ \
        {input.vcf} > {output}"""


rule whatshap_stats_phased_129S1_SvImJ:
    input:
        vcf=OUTPUT + 'vcf/{sid}.phased.vcf',
    output:
        OUTPUT + "reports/whatshap/{sid}.129S1_SvImJ.phased.stats.txt"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    shell:
        """whatshap stats \
        --sample 129S1_SvImJ \
        {input.vcf} > {output}"""


rule whatshap_haplotyping:
    input:
        bams=expand(OUTPUT + 'merged_bam/{cid}.GRCm39.bam', cid=cell_ids),
        vcf=OUTPUT + 'vcf/{sid}.phased.vcf.gz',
        ref=OUTPUT + 'references/GRCm39.fa'
    output:
        OUTPUT + 'whatshap/{cid}.{sid}.haplotaged.txt'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    threads:
        config['threads'] // 2
    shell:
         """whatshap haplotag {input.vcf} {input.bam} \
         --skip-missing-contigs --tag-supplementary --ignore-read-groups -r {input.ref} \
         --output-threads={threads} --output-haplotag-list  {output}  -o /dev/null"""


rule nanopolish_phase_reads:
    input:
        bam=OUTPUT + 'minimap2/{cid}.GRCm39.raw.bam',
        bamindex=OUTPUT + 'minimap2/{cid}.GRCm39.raw.bam.bai',
        fastq=OUTPUT + 'fastq/{cid}.raw.fastq',
        ref=OUTPUT + 'references/GRCm39.fa',
        vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
        vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
    output:
        OUTPUT + 'nanopolish/{cid}.{sid}.raw.phased.bam'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """nanopolish phase-reads --threads {threads} --reads {input.fastq} --bam {input.bam} --genome {input.ref} {input.vcf} > {output}"""




# rule phase_cast:
#     input:
#         vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
#         vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
#         ref=OUTPUT + 'references/GRCm39.fa',
#         bams=expand(OUTPUT + 'merged_bam/{cid}.GRCm39.bam', cid=cell_ids,),
#         bamindex=expand(OUTPUT + 'merged_bam/{cid}.GRCm39.bam.bai', cid=cell_ids,)
#     output:
#         vcf=OUTPUT + 'whatshap/{sid}.GRCm39.phased.vcf.gz',
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#         sid='|'.join([re.escape(x) for x in set(snp_ids)]),
#         rid='|'.join([re.escape(x) for x in set(ref_ids) if not x == 'GRCm39']),
#     shell:
#         """whatshap phase --ignore-read-groups \
#         -o {output} --reference={input.ref} \
#         --sample CAST_EiJ \
#         {input.vcf} {input.bams} """
# 
# 

# rule whatshap_phase:
#     input:
#         vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
#         vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
#         ref=OUTPUT + 'references/{rid}.fa',
#         bams=expand(OUTPUT + 'merged_bam/{cid}.{{rid}}.bam', cid=cell_ids,),
#         bamindex=expand(OUTPUT + 'merged_bam/{cid}.{{rid}}.bam.bai', cid=cell_ids,)
#     output:
#         vcf=OUTPUT + 'whatshap/{sid}.{rid}.phased.vcf.gz',
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#         sid='|'.join([re.escape(x) for x in set(snp_ids)]),
#         rid='|'.join([re.escape(x) for x in set(ref_ids)]),
#     params:
#         sample=lambda wildcards: wildcards.rid
#     shell:
#         """whatshap phase --ignore-read-groups \
#         -o {output} --reference={input.ref} \
#         --sample {params.sample} \
#         {input.vcf} {input.bams} """
# 

#rule whatshap_stats:
#    input:
#        OUTPUT + 'whatshap/{sid}.phased.vcf.gz',
#    output:
#        OUTPUT + "reports/whatshap_phase/{sid}.stats.txt"
#    wildcard_constraints:
#        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
#    shell:
#        """whatshap stats {input} > {output}"""
#
#
#PHASE_AGAINST = config['phase_against']
#

#
#
# rule beagle_phase:
#     input:
#         vcf=OUTPUT + 'snps/{sid}.vcf.gz',
#     output:
#         touch(OUTPUT + "snps/{sid}.done")
#     wildcard_constraints:
#         sid='|'.join([re.escape(x) for x in set(snp_ids)]),
#     params:
#         prefix=lambda wildcards: OUTPUT + "snps/" + wildcards.sid + ".",
#     threads:
#         config['threads'] // 2
#     shell:
#         """beagle gt={input.vcf} nthreads={threads} out={params.prefix} """
#     
#
# rule phase_vcf_py_popgen:
#     input:
#         vcf=OUTPUT + 'snps/{sid}.vcf.gz',
#     output:
#         OUTPUT + 'snps/{sid}.phased.vcf.gz',
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#     params:
#         jar=config['beagle_path']
#     shell:
#         """vcf_phase.py --vcf {input.vcf} \
#         --phase-algorithm beagle \
#         --beagle-path {params.jar} \
#         --out {output} """
# 
# 
# rule extractHAIRS:
#     input:
#         bam=OUTPUT + 'merged_bam/{cid}.{rid}.bam',
#         bamindex=OUTPUT + 'merged_bam/{cid}.{rid}.bam.bai',
#         ref=OUTPUT + 'references/{rid}.fa',
#         vcf=OUTPUT + 'vcf/{sid}.vcf',
#         vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
#     output:
#         OUTPUT + 'hapcut/{cid}.{rid}.{sid}.fragments'
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#         sid='|'.join([re.escape(x) for x in set(snp_ids)]),
#         rid='|'.join([re.escape(x) for x in set(ref_ids)]),
#     shell:
#         """extractHAIRS  --ONT 1 \
#         --ref {input.ref} \
#         --bam {input.bam} \
#         --VCF {input.vcf} \
#         --out {output} """
# 
# 
# rule merge_fragments:
#     input:
#         bams=expand(OUTPUT + 'hapcut/{cid}.{{rid}}.{{sid}}.fragments', cid=cell_ids)
#     output:
#         OUTPUT + 'hapcut/{rid}.{sid}.merged.fragments'
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#         sid='|'.join([re.escape(x) for x in set(snp_ids)]),
#         rid='|'.join([re.escape(x) for x in set(ref_ids)]),
#     shell:
#         """cat {input.bams} > {output}"""
# 
# 
# rule HAPCUT2:
#     input:
#         frag=OUTPUT + 'hapcut/{rid}.{sid}.merged.fragments',
#         vcf=OUTPUT + 'vcf/{sid}.vcf',
#         vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
#     output:
#         OUTPUT + 'hapcut/{rid}.{sid}.phased.VCF',
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#         sid='|'.join([re.escape(x) for x in set(snp_ids)]),
#         rid='|'.join([re.escape(x) for x in set(ref_ids)]),
#     params:
#         prefix=lambda wildcards: OUTPUT + "hapcut/" + wildcards.rid + "." + wildcards.sid
#     shell:
#         """HAPCUT2 --fragments {input.frag} \
#         --VCF {input.vcf} \
#         --output {params.prefix} """