rule epi2me_align:
    input:
        fastq=OUTPUT + 'fastq/{cid}.raw.fastq',
        bam=OUTPUT + 'minimap2/{cid}.{rid}.raw.bam',
        ref=OUTPUT + 'references/{rid}.fa',
    output:
        OUTPUT + 'epi2me_digest/{cid}.{rid}.bam'
    threads:
        config['threads'] // 4
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    params:
        enzyme=config['enzyme']
    shell:
        """scripts/epi2me-align-digested.sh -f {input.fastq} \
        -b {input.bam} -r {input.ref} -e {params.enzyme} \
        -t {threads} -o {output}"""


rule epi2me_coverage_by_cell:
    input:
        bam=OUTPUT + 'epi2me_digest/{cid}.{rid}.bam',
    output:
        OUTPUT + "reports/epi2me_coverage_by_cell/{cid}.{rid}.samtools.coverage.txt"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        """samtools sort {input} | samtools coverage - > {output}"""


rule epi2me_coverage_by_reference:
    input:
        expand(OUTPUT + 'epi2me_digest/{cid}.{{rid}}.bam', cid=cell_ids)
    output:
        OUTPUT + "reports/epi2me_coverage_by_reference/{rid}.samtools.coverage.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:   
        """samtools merge - {input} | samtools sort -| samtools coverage - > {output}"""


rule epi2me_annotate:
    input:
        bam=OUTPUT + 'epi2me_digest/{cid}.{rid}.bam'
    output:
        bam=OUTPUT + 'epi2me/{cid}.{rid}.ns.bam',
        pq=OUTPUT + 'epi2me/{cid}.{rid}.chromunity.parquet',
        sum=OUTPUT + 'epi2me/{cid}.{rid}.summary.json',
    threads:
        config['threads'] // 4
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    params:
        enzyme=config['enzyme'],
        prefix=lambda wildcards: OUTPUT + "epi2me/" + wildcards.cid + "." + wildcards.rid
    shell:
        """scripts/epi2me-annotate.sh -b {input} -t {threads} -o {params.prefix}"""