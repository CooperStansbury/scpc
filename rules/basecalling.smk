rule copy_pod5:
    input:
        pod5_df['file_path'].to_list()
    output:
        protected(pod5_df['out_path'].to_list()),
    run:
        from shutil import copyfile
        for i, fpath in enumerate(input):
    
            outPath = output[i]
            copyfile(fpath, outPath)


rule make_fast5:
    input:
        OUTPUT + "pod5/{cid}.pod5"
    output:
        directory=directory(OUTPUT + "fast5/{cid}.fast5"),
        flag=touch(OUTPUT + "flags/{cid}.merge.done"),
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """pod5 convert to_fast5 {input} \
        -t {threads} \
        --output {output.directory}"""


rule basecall:
    input:
        pod5=OUTPUT + "pod5/{cid}.pod5",
        flag=OUTPUT + "flags/{cid}.merge.done",
    output:
        protected(OUTPUT + "fastq/{cid}.raw.fastq")
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        dorado=config['dorado_path'],
        model=config['dorado_model'],
        qscore=config['dorado_min_qscore'],
    shell:
        """{params.dorado} basecaller {params.model} --emit-fastq --min-qscore {params.qscore} --no-trim {input.pod5} > {output}"""


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


rule f5c_index:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
        fast5=OUTPUT + "fast5/{cid}.fast5",
    output:
        OUTPUT + "fastq/{cid}.raw.fastq.index",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 2
    shell:
        """f5c index -t {threads} --iop 10 -d {input.fast5} {input.fastq}"""


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
        OUTPUT + 'NlaIII_locations/{cid}.csv'
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
        OUTPUT + "read_stats/{cid}.read_lengths.pq",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """python scripts/get_read_lengths.py {input} {output}"""



rule get_restriction_counts:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        OUTPUT + "read_stats/{cid}.restriction_counts.pq",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        cutter=config['enzyme'],
    shell:
        """python scripts/get_restriction_count.py {input} {params.cutter} {output}"""




