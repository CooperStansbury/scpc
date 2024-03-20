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