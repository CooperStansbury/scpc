rule copy_references:
    input:
        ref_df['file_path'].to_list()
    output:
        ref_df['out_path'].to_list(),
    run:
        from shutil import copyfile
        for i, fpath in enumerate(input):
    
            outPath = output[i]
            copyfile(fpath, outPath)


rule make_chromsizes:
    input:
        OUTPUT + 'references/{rid}.fa',
    output:
        OUTPUT + 'reports/chromsizes/{rid}.chrom.sizes'
    wildcard_constraints:
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        "faidx {input} -i chromsizes > {output}"

 
rule minimap2_index:
    input:
        refgenome=OUTPUT + 'references/{rid}.fa'
    output:
        OUTPUT + 'references/{rid}.mmi'
    threads:
        config['threads'] // 2
    wildcard_constraints:
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        "minimap2 -t {threads} -d {output} {input.refgenome}"


# rule bwa_index:
#     input:
#         refgenome=OUTPUT + 'references/{rid}.fa'
#     output:
#         OUTPUT + 'references/{rid}.bwt'
#     shell:
#         "bwa index {input}"


rule copy_gtfs:
    input:
        gtf_df['file_path'].to_list()
    output:
        gtf_df['out_path'].to_list(),
    run:
        from shutil import copyfile
        for i, fpath in enumerate(input):
    
            outPath = output[i]
            copyfile(fpath, outPath)


# # custom rule for adding a masked references
# rule mask_reference:
#     input:
#         ref=OUTPUT + 'references/GRCm39.fa',
#         vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
#         vcfindex=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
#     output:
#         OUTPUT + 'references/GRCm39masked.fa'    
#     wildcard_constraints:
#         sid='|'.join([re.escape(x) for x in set(snp_ids)]),
#     shell:
#         """bedtools maskfasta -fi {input.ref} \
#         -bed {input.vcf} \
#         -fo {output}
#         """
# 
# ref_ids.append("GRCm39masked")