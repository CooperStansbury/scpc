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
include: "rules/mapping.smk"
include: "rules/haplotyping.smk"

  
################ ALL RULES ################
rule all:
    input:
        expand(f"{OUTPUT}fastq/{{cid}}.fastq", cid=cell_ids),
        OUTPUT + 'reports/seqkit/raw.fastq.report.txt',
        OUTPUT + 'reports/seqkit/digested.fastq.report.txt',
        expand(f"{OUTPUT}references/{{rid}}.mmi", rid=ref_ids),
        expand(f"{OUTPUT}references/{{rid}}.bwt", rid=ref_ids),
        expand(f"{OUTPUT}reports/chromsizes/{{rid}}.chrom.sizes", rid=ref_ids),
        expand(f"{OUTPUT}gtf/{{gid}}.gtf.gz", gid=gtf_ids),
        expand(f"{OUTPUT}digested_fastq/{{cid}}.fastq", cid=cell_ids),
        expand(f"{OUTPUT}restriction_sites/{{cid}}.sites.pq", cid=cell_ids),
        expand(f"{OUTPUT}merged_bam/{{cid}}.{{rid}}.bam.bai", cid=cell_ids, rid=ref_ids),
        expand(f"{OUTPUT}reports/coverage/{{cid}}.{{rid}}.{{dig}}.samtools.coverage.txt", cid=cell_ids, rid=ref_ids, dig=['full', 'digested']),
        expand(f"{OUTPUT}reports/coverage_by_cell/{{cid}}.{{rid}}.samtools.coverage.txt", cid=cell_ids, rid=ref_ids),
        expand(f"{OUTPUT}reports/flagstat/{{cid}}.{{rid}}.flagstat.tsv", cid=cell_ids, rid=ref_ids),
        expand(f"{OUTPUT}reports/coverage_by_reference/{{rid}}.samtools.coverage.txt", rid=ref_ids),
        expand(f"{OUTPUT}align_table/{{cid}}.{{rid}}.alignments.pq", cid=cell_ids, rid=ref_ids),
        expand(f"{OUTPUT}vcf/{{sid}}.vcf.gz.tbi", sid=snp_ids),
        expand(f"{OUTPUT}duplicates/{{cid}}.{{rid}}.bam", cid=cell_ids, rid=ref_ids),
        expand(f"{OUTPUT}vcf/{{sid}}.snps.tsv", sid=snp_ids),
        expand(f"{OUTPUT}whatshap/{{sid}}.GRCm39.phased.vcf.gz", sid=snp_ids),
        expand(f"{OUTPUT}reports/whatshap/{{sid}}.vcf.stats.txt", sid=snp_ids),

    
    


rule haplotagging:
    input:
        expand(f"{OUTPUT}whatshap/{{sid}}.GRCm39.phased.vcf.gz", sid=snp_ids),
        expand(f"{OUTPUT}vcf/{{sid}}.snps.tsv", sid=snp_ids),
        expand(f"{OUTPUT}vcf/{{sid}}.gatk.table.tsv", sid=snp_ids),
        expand(f"{OUTPUT}snp_positions/{{cid}}.{{rid}}.{{sid}}.snps.pq", cid=cell_ids, rid=ref_ids, sid=snp_ids),
        expand(f"{OUTPUT}hapcut/GRCm39.{{sid}}.merged.fragments", sid=snp_ids),
        expand(f"{OUTPUT}whatshap/{{sid}}.phased.vcf.gz", sid=snp_ids),
        expand(f"{OUTPUT}reports/whatshap_phase/{{sid}}.phase_stats.txt", sid=snp_ids),
        expand(f"{OUTPUT}whatshap/{{sid}}.{{rid}}.{{sid}}.haplotaged.txt", cid=cell_ids, rid=ref_ids, sid=snp_ids),
        expand(f"{OUTPUT}hapcut/{{cid}}.GRCm39.{{sid}}.fragments", cid=cell_ids, sid=snp_ids),
        expand(f"{OUTPUT}hapcut/GRCm39.{{sid}}.phased.VCF", sid=snp_ids),



rule make_alignment_table:
    input:
        bam=OUTPUT + 'merged_bam/{cid}.{rid}.bam',
        frags=OUTPUT + 'restriction_sites/{cid}.sites.pq',
    output:
        OUTPUT + "align_table/{cid}.{rid}.alignments.pq"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        """python scripts/get_alignment_table.py {input.bam} {input.frags} {output}"""


rule report_snps:
    input:
        bam=OUTPUT + 'merged_bam/{cid}.{rid}.bam',
        vcf=OUTPUT + 'vcf/{sid}.vcf.gz',
        vcf_idx=OUTPUT + 'vcf/{sid}.vcf.gz.tbi',
        sites=OUTPUT + "restriction_sites/{cid}.sites.pq",
    output:
        OUTPUT + "snp_positions/{cid}.{rid}.{sid}.snps.pq"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        sid='|'.join([re.escape(x) for x in set(snp_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
    shell:
        """python scripts/get_snp_positions.py {input.vcf} {input.bam} {input.sites} {output}"""


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
