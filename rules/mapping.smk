rule minimap_align:
    input:
        fastq=OUTPUT + "fastq/{cid}.{cond}.fastq",
        ref=OUTPUT + 'references/{rid}.fa',
    output:
        bam=OUTPUT + 'minimap2/{cid}.{rid}.{cond}.bam'
    threads:
        config['threads'] // 2
    params:
        args=config['minimap_args']
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
        cond='raw|digested',
    shell:
        """minimap2 {params.args} -t {threads} {input.ref} {input.fastq} \
        | samtools sort -O bam -o {output.bam} """


rule samtools_index:
   input:
       OUTPUT + 'minimap2/{cid}.{rid}.{cond}.bam'
   output:
       OUTPUT + 'minimap2/{cid}.{rid}.{cond}.bam.bai'
   threads:
       config['threads'] // 2
   wildcard_constraints:
       cid='|'.join([re.escape(x) for x in set(cell_ids)]),
       rid='|'.join([re.escape(x) for x in set(ref_ids)]),
       cond='raw|digested',
   shell:
       "samtools index -@ {threads} {input}"


rule merge_bams:
    input:
        full=OUTPUT + 'minimap2/{cid}.{rid}.raw.bam',
        digested=OUTPUT + 'minimap2/{cid}.{rid}.digested.bam',
    output:    
        OUTPUT + 'merged_bam/{cid}.{rid}.bam'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """samtools merge -@ {threads} -o {output} {input.full} {input.digested} """


rule samtools_index_merged:
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
        bam=OUTPUT + 'minimap2/{cid}.{rid}.{cond}.bam'
    output:
        OUTPUT + "reports/coverage/{cid}.{rid}.{cond}.samtools.coverage.txt"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
        digest='raw|digested',
    shell:
        """samtools coverage {input} > {output}"""


rule coverage_by_cell:
    input:
        bam=expand(OUTPUT + 'minimap2/{{cid}}.{{rid}}.{digest}.bam', digest=['raw', 'digested'])
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


rule samtools_stats:
    input:
        bam=OUTPUT + 'merged_bam/{cid}.{rid}.bam'
    output:
        OUTPUT + 'reports/stats/{cid}.{rid}.samtools.stats.txt'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    threads:
        config['threads'] // 4 
    shell:
        """samtools stats -@ {threads} {input.bam} > {output}"""