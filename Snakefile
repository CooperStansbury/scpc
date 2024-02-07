import pandas as pd
import yaml
from pathlib import Path
import re
import os
import sys
from utils import utils
from snakemake.utils import Paramspace
from tabulate import tabulate

BASE_DIR = Path(workflow.basedir)
configfile: str(BASE_DIR) + "/config/config.yaml"
OUTPUT = config['output_path']

# load the basecalls
fastq_paths = os.path.abspath(config['fastq_paths'])
fastq_df = utils.load_fastq_df(fastq_paths, OUTPUT)
cell_ids = fastq_df['cell_id'].to_list() # wildcard constraints

print(f"\n======== INPUT FILES ========")
print(tabulate(fastq_df[['cell_id', 'basename']], 
      headers='keys',
      showindex=False,
      tablefmt='psql'))

# build the references
ref_paths = os.path.abspath(config['ref_paths'])
ref_df = utils.load_ref_df(ref_paths, OUTPUT)
ref_ids = ref_df['ref_id'].to_list() # wildcard constraints

print(f"\n======== REFERENCES ========")
print(tabulate(ref_df[['ref_id', 'basename']], 
      headers='keys',
      showindex=False,
      tablefmt='psql'))


# build the snps
snps_paths = os.path.abspath(config['snps_paths'])
snp_df = utils.load_snp_df(snps_paths, OUTPUT)
snp_ids = snp_df['snp_id'].to_list() # wildcard constraints

print(f"\n======== SNPS ========")
print(tabulate(snp_df[['snp_id', 'basename']], 
      headers='keys',
      showindex=False,
      tablefmt='psql'))


# build the annotations
gtf_paths = os.path.abspath(config['gtf_paths'])
gtf_df = utils.load_gtf_df(gtf_paths, OUTPUT)
gtf_ids = gtf_df['gtf_id'].to_list() # wildcard constraints

print(f"\n======== GTF ========")
print(tabulate(gtf_df[['gtf_id', 'basename']], 
      headers='keys',
      showindex=False,
      tablefmt='psql'))



################ RULE FILES ################
include: "rules/references.smk"

  
################ ALL RULES ################
rule all:
    input:
        expand(f"{OUTPUT}fastq/{{cid}}.fastq", cid=cell_ids),
        OUTPUT + 'reports/seqkit/raw.fastq.report.txt',
        expand(f"{OUTPUT}references/{{rid}}.mmi", rid=ref_ids),
        expand(f"{OUTPUT}references/{{rid}}.chrom.sizes", rid=ref_ids),
        expand(f"{OUTPUT}snps/{{sid}}.vcf.gz.tbi", sid=snp_ids),
        expand(f"{OUTPUT}gtf/{{gid}}.gtf.gz", gid=gtf_ids),
        expand(f"{OUTPUT}digested_fastq/{{cid}}.fastq", cid=cell_ids),

        # OUTPUT + 'reports/samtools.alignment.coverage.txt',
        # expand(f"{OUTPUT}reports/flagstat/{{cid}}.flagstat.tsv", cid=cell_ids),
        # expand(f"{OUTPUT}sites/{{cid}}.cut.sites.pq", cid=cell_ids),
        # expand(f"{OUTPUT}bam/{{cid}}.bam.bai", cid=cell_ids),
        # expand(f"{OUTPUT}reports/samtools_coverage/{{cid}}.txt", cid=cell_ids),
        # expand(f"{OUTPUT}align_table/{{cid}}.alignments.pq", cid=cell_ids),
        # expand(f"{OUTPUT}concatemers/{{cid}}.concatemers.pq", cid=cell_ids),
        # expand(f"{OUTPUT}haplotypes/{{cid}}.txt", cid=cell_ids),
        # OUTPUT + "reports/seqkit.filtered.fastq.report.txt",
        # expand(f"{OUTPUT}filtered_fastq/{{cid}}.fastq", cid=cell_ids),


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
        report=OUTPUT + "reports/digest/{cid}.sites.pq",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        cutter=config['enzyme'],
    shell:
        """python scripts/digest.py {input} {params.cutter} \
        {output.fastq} {output.report}"""


rule find_sites:
    input:
        fastq=OUTPUT + "fastq/{cid}.fastq"
    output:
        sites=OUTPUT + "sites/{cid}.cut.sites.pq",
        ligated_reads=OUTPUT + "sites/{cid}.ligated.reads.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        cutter=config['cutter'],
    shell:
        """python scripts/get_cut_sites.py {input.fastq} \
           {params.cutter} {output.sites} {output.ligated_reads}"""


rule filter_fastq:
    input:
        fastq=OUTPUT + "fastq/{cid}.fastq",
        ids=OUTPUT + "sites/{cid}.ligated.reads.txt",
    output:
        OUTPUT + "filtered_fastq/{cid}.fastq"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """seqkit grep -f {input.ids} -j {threads} -o {output} {input.fastq}"""


rule filtered_fastq_report:
    input:
        expand(f"{OUTPUT}filtered_fastq/{{cid}}.fastq", cid=cell_ids),
    output:
        OUTPUT + "reports/seqkit.filtered.fastq.report.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""


rule minimap_align:
    input:
        fastq=OUTPUT + "fastq/{cid}.fastq",
        ref=OUTPUT + 'references/reference.fa',
    output:
        OUTPUT + 'bam/{cid}.bam'
    threads:
        config['threads'] // 2
    params:
        args=config['minimap_args']
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """minimap2 {params.args} -t {threads} {input.ref} {input.fastq} \
        | samtools sort -O bam -o {output}"""


rule samtools_index:
    input:
        OUTPUT + 'bam/{cid}.bam'
    output:
        OUTPUT + 'bam/{cid}.bam.bai'
    threads:
        config['threads'] // 2
    shell:
        "samtools index -@ {threads} {input}"


rule samtool_flagstat:
    input:
        bam=OUTPUT + 'bam/{cid}.bam',
    output:
        OUTPUT + 'reports/flagstat/{cid}.flagstat.tsv'
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4 
    shell:
        """samtools flagstat -@ {threads} -O 'tsv' {input.bam} > {output}"""


rule samtools_coverage:
    input:
        bam=OUTPUT + 'bam/{cid}.bam',
    output:
        OUTPUT + "reports/samtools_coverage/{cid}.txt"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """samtools coverage {input} > {output}"""


rule all_coverage:
    input:
        expand(f"{OUTPUT}bam/{{cid}}.bam", cid=cell_ids),
    output:
        OUTPUT + "reports/samtools.alignment.coverage.txt"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """samtools coverage {input} > {output}"""


rule make_alignment_table:
    input:
        bam=OUTPUT + 'bam/{cid}.bam',
    output:
        OUTPUT + "align_table/{cid}.alignments.pq"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """python scripts/get_alignment_table.py {input} {output}"""


rule get_concatemers:
    input:
        OUTPUT + "align_table/{cid}.alignments.pq"
    output:
        OUTPUT + "concatemers/{cid}.concatemers.pq"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        mapq=config['mapq_threshold']
    shell:
        """python scripts/resolve_concatemers.py {input} {params.mapq} {output}"""
    


# 
# rule get_cell_ids:
#     input:
#         OUTPUT + 'references/snps.vcf.gz'
#     output:
#         OUTPUT + 'references/cell_ids.txt'
#     shell:
#         """bcftools query -l {input} > {output}"""
# 

# rule tag_reads:
#     input:
#         vcf=OUTPUT + 'references/snps.vcf.gz',
#         ref=OUTPUT + 'references/reference.fa',
#         bam=OUTPUT + 'bam/{cid}.bam',
#         cell_ids=OUTPUT + 'references/vcf.cell_ids.txt',
#     output:
#         tab=OUTPUT + 'haplotypes/{cid}.txt',
#         bam=OUTPUT + 'haplotypes/{cid}.bam',
#     wildcard_constraints:
#         cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#     threads:
#         config['threads'] // 4
#     params:
#         cell_ids=", ".join([x.strip() for x in open(OUTPUT + 'references/vcf.cell_ids.txt')])
#     shell:
#         """whatshap haplotag {input.vcf} {input.bam} \
#         --skip-missing-contigs \
#         --ignore-read-groups \
#         --sample 'CAST_EiJ' --output-threads={threads} --output-haplotag-list {output.tab} --output {output.bam}"""
#     
