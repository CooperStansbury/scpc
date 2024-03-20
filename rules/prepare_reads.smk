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
        fastq=OUTPUT + "fastq/{cid}.digested.fastq",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        cutter=config['enzyme'],
    shell:
        """python scripts/digest.py {input} \
           {params.cutter} {output.fastq} """


rule fastq_report:
    input:
        expand(OUTPUT + "fastq/{cid}.{{cond}}.fastq", cid=cell_ids),
    output:
        OUTPUT + "reports/seqkit/{cond}.fastq.report.txt",
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
    

rule locate_barcodes:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
        codes=OUTPUT + "resources/barcodes.fasta"
    output:
        OUTPUT + 'barcode_locations/{cid}.csv'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    params:
        m=config['n_mismatches']
    shell:
        """seqkit locate -m {params.m} -f {input.codes} {input.fastq} -j {threads} > {output} """


rule locate_enzyme:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        OUTPUT + 'enzyme_locations/{cid}.csv'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """seqkit locate -p CATG {input.fastq} -j {threads} > {output} """
    


rule run_sequence_report:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
        bc=OUTPUT + "resources/barcodes.txt",
    output:
        OUTPUT + "reports/sequence_reports/{cid}.report.pq",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        enzyme=config['enzyme'],
    shell:
        """python scripts/sequence_report.py {input.fastq} {input.bc} {params.enzyme} {output} """


rule fastqc_report:
    input:
        OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        html=OUTPUT + "reports/fastqc/{cid}.raw_fastqc.html",
        zip=OUTPUT + "reports/fastqc/{cid}.raw_fastqc.zip" 
    params: f"--quiet --contaminants {config['fastqc_contaminants']} --adapters {config['fastqc_adapters']}"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads: 
        config['threads'] // 4
    wrapper:
        "v1.14.1/bio/fastqc"


rule get_read_lengths:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        OUTPUT + "read_stats/{cid}.read_lengths.parquet",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """python scripts/get_read_lengths.py {input} {output}"""



rule get_restriction_counts:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        OUTPUT + "read_stats/{cid}.restriction_counts.parquet",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        cutter=config['enzyme'],
    shell:
        """python scripts/get_restriction_count.py {input} {params.cutter} {output}"""