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

# load the fastq files 
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
include: "rules/prepare_reads.smk"
include: "rules/references.smk"
include: "rules/mapping.smk"
include: "rules/epi2me.smk"
include: "rules/haplotyping.smk"

  
################ ALL RULES ################
rule all:
    input:
        expand(f"{OUTPUT}gtf/{{gid}}.gtf.gz", gid=gtf_ids),
        expand(f"{OUTPUT}references/{{rid}}.mmi", rid=ref_ids),
        expand(f"{OUTPUT}fastq/{{cid}}.raw.fastq", cid=cell_ids),
        expand(f"{OUTPUT}fastq/{{cid}}.digested.fastq", cid=cell_ids),
        # expand(f"{OUTPUT}reports/chromsizes/{{rid}}.chrom.sizes", rid=ref_ids),        
        # expand(f"{OUTPUT}minimap2/{{cid}}.{{rid}}.{{cond}}.bam", cid=cell_ids, rid=ref_ids, cond=['raw', 'digested']),
        # expand(f"{OUTPUT}reports/seqkit/{{cond}}.fastq.report.txt", cond=['raw', 'digested']),
        # expand(f"{OUTPUT}reports/fastqc/{{cid}}.raw_fastqc.html", cid=cell_ids),
        # expand(f"{OUTPUT}reports/coverage/{{cid}}.{{rid}}.{{dig}}.txt", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),
        # expand(f"{OUTPUT}reports/flagstat/{{cid}}.{{rid}}.{{dig}}.tsv", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),
        # expand(f"{OUTPUT}reports/stats/{{cid}}.{{rid}}.{{dig}}.txt", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),
        # expand(f"{OUTPUT}read_stats/{{cid}}.read_lengths.parquet", cid=cell_ids,),
        # expand(f"{OUTPUT}read_stats/{{cid}}.restriction_counts.parquet", cid=cell_ids,),
        # expand(f"{OUTPUT}barcode_locations/{{cid}}.csv", cid=cell_ids),
        # expand(f"{OUTPUT}enzyme_locations/{{cid}}.csv", cid=cell_ids),
        # expand(f"{OUTPUT}reports/epi2me_coverage/{{cid}}.{{rid}}.txt", cid=cell_ids, rid=ref_ids),
        # expand(f"{OUTPUT}epi2me_digest/{{cid}}.{{rid}}.bam", cid=cell_ids, rid=ref_ids),
        # expand(f"{OUTPUT}epi2me/{{cid}}.{{rid}}.ns.bam", cid=cell_ids, rid=ref_ids),
        # expand(f"{OUTPUT}duplicates/{{cid}}.{{rid}}.{{dig}}.bam", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),
        # expand(f"{OUTPUT}align_table/{{cid}}.{{rid}}.{{dig}}.parquet", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),
    

rule archive:
    input:
        expand(f"{OUTPUT}reports/sequence_reports/{{cid}}.report.pq", cid=cell_ids),
        expand(f"{OUTPUT}reports/nanostat/{{cid}}.{{rid}}.{{dig}}.tsv", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),


rule bam2table:
    input:
        bam=OUTPUT + "minimap2/{cid}.{rid}.{cond}.bam"
    output:
        OUTPUT + "align_table/{cid}.{rid}.{cond}.parquet"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
        rid='|'.join([re.escape(x) for x in set(ref_ids)]),
        cond='raw|digested',
    params:
        sep="_",
    shell:
        """python scripts/bam2table.py {input.bam} {params.sep} {output}"""
