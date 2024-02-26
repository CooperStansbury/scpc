rule copy_pod5:
    input:
        pod5_df['file_path'].to_list()
    output:
        pod5_df['out_path'].to_list(),
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
        OUTPUT + "pod5/{cid}.pod5"
    output:
        OUTPUT + "fastq/{cid}.raw.fastq"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        dorado=config['dorado_path'],
        model=config['dorado_model'],
        qscore=config['dorado_min_qscore'],
    shell:
        """{params.dorado} basecaller {params.model} --emit-fastq --min-qscore {params.qscore} --no-trim {input} > {output}"""


# rule nanopolish_index:
#     input:
#         fastq=OUTPUT + "fastq/{cid}.raw.fastq",
#         fast5=OUTPUT + "fast5/{cid}.fast5",
#     output:
#         OUTPUT + "fastq/{cid}.fastq.index",
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#     shell:
#         """nanopolish index -d {input.fast5} {input.fastq}"""
#     

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
        OUTPUT + "fastq/{cid}.fastq.index",
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


rule report_cut_sites:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        OUTPUT + "restriction_sites/{cid}.sites.pq",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        cutter=config['enzyme'],
    shell:
        """python scripts/get_cut_sites.py {input} \
           {params.cutter} {output} """

