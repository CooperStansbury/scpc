rule copy_fastq:
    input:
        fastq_df['file_path'].to_list()
    output:
        protected(fastq_df['out_path'].to_list()),
    run:
        from shutil import copyfile
        for i, fpath in enumerate(input):
    
            outPath = output[i]
            copyfile(fpath, outPath)


rule digest_fastq:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        fastq=OUTPUT + "digest/{cid}.digested.fastq",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        cutter=config['enzyme'],
    shell:
        """python scripts/digest.py {input} \
           {params.cutter} {output.fastq} """


rule fastq_raw_report:
    input:
        expand(OUTPUT + "fastq/{cid}.raw.fastq", cid=cell_ids),
    output:
        OUTPUT + "reports/seqkit/raw.fastq.report.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""


rule fastq_digested_report:
    input:
        expand(OUTPUT + "digest/{cid}.digested.fastq", cid=cell_ids),
    output:
        OUTPUT + "reports/seqkit/digested.fastq.report.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""


rule get_barcode_file:
    input:
        config['barcode_path'],
    output:
        OUTPUT + "resources/barcodes.txt"
    shell:
        """cp {input} {output}"""


rule get_barcode_fasta:
    input:
        OUTPUT + "resources/barcodes.txt"
    output:
        OUTPUT + "resources/barcodes.fasta"
    shell:
        """python scripts/barcode_fasta.py {input} {output}"""
    

# rule locate_barcodes:
#     input:
#         fastq=OUTPUT + "fastq/{cid}.raw.fastq",
#         codes=OUTPUT + "resources/barcodes.fasta"
#     output:
#         OUTPUT + 'barcode_locations/{cid}.csv'
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#     threads:
#         config['threads'] // 4
#     params:
#         m=config['n_mismatches']
#     shell:
#         """seqkit locate -m {params.m} -f {input.codes} {input.fastq} -j {threads} > {output} """
# 
# 
# rule locate_enzyme:
#     input:
#         fastq=OUTPUT + "fastq/{cid}.raw.fastq",
#     output:
#         OUTPUT + 'enzyme_locations/{cid}.csv'
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#     threads:
#         config['threads'] // 4
#     shell:
#         """seqkit locate -p CATG {input.fastq} -j {threads} > {output} """
# 
# 
#
# 
# 
# 
# rule get_read_lengths:
#     input:
#         fastq=OUTPUT + "fastq/{cid}.raw.fastq",
#     output:
#         OUTPUT + "read_stats/{cid}.read_lengths.parquet",
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#     shell:
#         """python scripts/get_read_lengths.py {input} {output}"""
# 
# 
# rule get_restriction_counts:
#     input:
#         fastq=OUTPUT + "fastq/{cid}.raw.fastq",
#     output:
#         OUTPUT + "read_stats/{cid}.restriction_counts.parquet",
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#     params:
#         cutter=config['enzyme'],
#     shell:
#         """python scripts/get_restriction_count.py {input} {params.cutter} {output}"""