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
        OUTPUT + 'references/{rid}.chrom.sizes'
    wildcard_constraints:
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        "samtools faidx {input} | cut -f1,2 > {output}"
        
        
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
        OUTPUT + 'snps/{sid}.vcf.gz'
    output: 
        OUTPUT + 'snps/{sid}.vcf.gz.tbi'
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
    shell:
        "tabix -fp vcf {input}"


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