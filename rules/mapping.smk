rule copy_fastq:
    input:
        fastq_df['file_path'].to_list()
    output:
        fastq_df['out_path'].to_list(),
    run:
        from shutil import copyfile
        for i, fpath in enumerate(input):
    
            outPath = output[i]
            copyfile(fpath, outPath)


rule raw_fastq_report:
    input:
        expand(f"{OUTPUT}fastq/{{cid}}.fastq", cid=cell_ids),
    output:
        OUTPUT + "reports/seqkit/raw.fastq.report.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""


rule digest_fastq:
    input:
        fastq=OUTPUT + "fastq/{cid}.fastq",
    output:
        fastq=OUTPUT + "digested_fastq/{cid}.fastq",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        cutter=config['enzyme'],
    shell:
        """python scripts/digest.py {input} \
           {params.cutter} {output.fastq} """


rule report_cut_sites:
    input:
        fastq=OUTPUT + "fastq/{cid}.fastq",
    output:
        OUTPUT + "restriction_sites/{cid}.sites.pq",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        cutter=config['enzyme'],
    shell:
        """python scripts/get_cut_sites.py {input} \
           {params.cutter} {output} """


rule digested_fastq_report:
    input:
        expand(f"{OUTPUT}digested_fastq/{{cid}}.fastq", cid=cell_ids),
    output:
        OUTPUT + "reports/seqkit/digested.fastq.report.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""


rule minimap_align_full:
    input:
        fastq=OUTPUT + "fastq/{cid}.fastq",
        ref=OUTPUT + 'references/{rid}.fa',
    output:
        OUTPUT + 'bam/{cid}.{rid}.full.bam'
    threads:
        config['threads'] // 2
    params:
        args=config['minimap_args']
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        """minimap2 {params.args} -t {threads} {input.ref} {input.fastq} \
        | samtools sort -O bam -o {output}"""


rule minimap_align_digested:
    input:
        fastq=OUTPUT + "digested_fastq/{cid}.fastq",
        ref=OUTPUT + 'references/{rid}.fa',
    output:
        OUTPUT + 'bam/{cid}.{rid}.digested.bam'
    threads:
        config['threads'] // 2
    params:
        args=config['minimap_args']
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        """minimap2 {params.args} -t {threads} {input.ref} {input.fastq} \
        | samtools sort -O bam -o {output}"""


rule bwa_align_full:
    input:
        fastq=OUTPUT + "fastq/{cid}.fastq",
        ref=OUTPUT + 'references/{rid}.fa',
    output:
    threads:
        config['threads'] // 2
    params:
        args=config['bwa_args']
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        """bwa {params.args} -t {threads} \
        | samtools sort -O bam -o {output}"""



rule merge_bams:
    input:
        full=OUTPUT + 'bam/{cid}.{rid}.full.bam',
        digested=OUTPUT + 'bam/{cid}.{rid}.digested.bam',
    output:    
        OUTPUT + 'merged_bam/{cid}.{rid}.bam'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """samtools merge -@ {threads} -o {output} {input.full} {input.digested} """


rule samtools_index:
    input:
        OUTPUT + 'merged_bam/{cid}.{rid}.bam'
    output:
        OUTPUT + 'merged_bam/{cid}.{rid}.bam.bai'
    threads:
        config['threads'] // 2
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        "samtools index -@ {threads} {input}"


rule mark_dupes:
    input:
        OUTPUT + 'merged_bam/{cid}.{rid}.bam'
    output:
        bam=OUTPUT + 'duplicates/{cid}.{rid}.bam',
        report=OUTPUT + 'reports/duplicates/{cid}.{rid}.txt',
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        """gatk MarkDuplicates I={input} O={output.bam} M={output.report} """


rule samtool_flagstat:
    input:
        bam=OUTPUT + 'merged_bam/{cid}.{rid}.bam'
    output:
        OUTPUT + 'reports/flagstat/{cid}.{rid}.flagstat.tsv'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    threads:
        config['threads'] // 4 
    shell:
        """samtools flagstat -@ {threads} -O 'tsv' {input.bam} > {output}"""


rule samtools_coverage_by_digest:
    input:
        bam=OUTPUT + 'bam/{cid}.{rid}.{digest}.bam'
    output:
        OUTPUT + "reports/coverage/{cid}.{rid}.{digest}.samtools.coverage.txt"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
        digest='full|digested',
    shell:
        """samtools coverage {input} > {output}"""


rule coverage_by_cell:
    input:
        bam=expand(OUTPUT + 'bam/{{cid}}.{{rid}}.{digest}.bam', digest=['full', 'digested'])
    output:
        OUTPUT + "reports/coverage_by_cell/{cid}.{rid}.samtools.coverage.txt"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        """samtools coverage {input} > {output}"""


rule coverage_by_reference:
    input:
        expand(OUTPUT + 'merged_bam/{cid}.{{rid}}.bam', cid=cell_ids)
    output:
        OUTPUT + "reports/coverage_by_reference/{rid}.samtools.coverage.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:   
        """samtools coverage {input} > {output}"""