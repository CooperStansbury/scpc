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


rule nanostat_aligned:
   input:
       OUTPUT + 'minimap2/{cid}.{rid}.{cond}.bam'
   output:
       OUTPUT + "reports/nanostat/{cid}.{rid}.{cond}.tsv"
   wildcard_constraints:
       cid='|'.join([re.escape(x) for x in set(cell_ids)]),
       rid='|'.join([re.escape(x) for x in set(ref_ids)]),
        cond='raw|digested',
   threads:
       config['threads'] // 4
   shell:
       """NanoStat -t {threads} --bam {input} -n {output}"""


rule samtool_flagstat:
    input:
        bam=OUTPUT + 'minimap2/{cid}.{rid}.{cond}.bam'
    output:
        OUTPUT + 'reports/flagstat/{cid}.{rid}.{cond}.tsv'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
        cond='raw|digested',
    threads:
        config['threads'] // 4 
    shell:
        """samtools flagstat -@ {threads} -O 'tsv' {input.bam} > {output}"""


rule samtools_stats:
    input:
        bam=OUTPUT + 'minimap2/{cid}.{rid}.{cond}.bam'
    output:
        OUTPUT + 'reports/stats/{cid}.{rid}.{cond}.txt'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
        cond='raw|digested',
    threads:
        config['threads'] // 4 
    shell:
        """samtools stats -@ {threads} {input.bam} > {output}"""


rule samtools_coverage:
    input:
        bam=OUTPUT + 'minimap2/{cid}.{rid}.{cond}.bam'
    output:
        OUTPUT + "reports/coverage/{cid}.{rid}.{cond}.txt"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
        digest='raw|digested',
    shell:
        """samtools coverage {input} > {output}"""


rule mark_duplicates:
    input:
         OUTPUT + 'minimap2/{cid}.{rid}.{cond}.bam'
    output:
         bam=OUTPUT + 'duplicates/{cid}.{rid}.{cond}.bam',
         report=OUTPUT + 'reports/duplicates/{cid}.{rid}.{cond}.txt',
    wildcard_constraints:
         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
         rid='|'.join([re.escape(x) for x in set(ref_ids)]),
         digest='raw|digested',
    shell:
         """gatk MarkDuplicates I={input} O={output.bam} M={output.report} """

