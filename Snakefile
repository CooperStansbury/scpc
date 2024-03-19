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
pod5_paths = os.path.abspath(config['pod5_paths'])
pod5_df = utils.load_pod5_df(pod5_paths, OUTPUT)
cell_ids = pod5_df['cell_id'].to_list() # wildcard constraints

print(f"\n======== INPUT FILES ========")
print(tabulate(pod5_df[['cell_id', 'basename']], 
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
include: "rules/basecalling.smk"
include: "rules/references.smk"
include: "rules/mapping.smk"
include: "rules/epi2me.smk"
include: "rules/haplotyping.smk"

  
################ ALL RULES ################
rule all:
    input:
        expand(f"{OUTPUT}reports/chromsizes/{{rid}}.chrom.sizes", rid=ref_ids),        
        expand(f"{OUTPUT}minimap2/{{cid}}.{{rid}}.{{cond}}.bam", cid=cell_ids, rid=ref_ids, cond=['raw', 'digested']),
        expand(f"{OUTPUT}reports/seqkit/{{cond}}.fastq.report.txt", cond=['raw', 'digested']),
        expand(f"{OUTPUT}reports/fastqc/{{cid}}.raw_fastqc.html", cid=cell_ids),
        expand(f"{OUTPUT}reports/coverage/{{cid}}.{{rid}}.{{dig}}.txt", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),
        expand(f"{OUTPUT}reports/flagstat/{{cid}}.{{rid}}.{{dig}}.tsv", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),
        expand(f"{OUTPUT}reports/stats/{{cid}}.{{rid}}.{{dig}}.txt", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),
        expand(f"{OUTPUT}read_stats/{{cid}}.read_lengths.pq", cid=cell_ids,),
        expand(f"{OUTPUT}read_stats/{{cid}}.restriction_counts.pq", cid=cell_ids,),
        expand(f"{OUTPUT}barcode_locations/{{cid}}.csv", cid=cell_ids),
        expand(f"{OUTPUT}NlaIII_locations/{{cid}}.csv", cid=cell_ids),
        # expand(f"{OUTPUT}reports/epi2me_coverage_by_cell/{{cid}}.{{rid}}.samtools.coverage.txt", cid=cell_ids, rid=ref_ids),
        # expand(f"{OUTPUT}duplicates/{{cid}}.{{rid}}.bam", cid=cell_ids, rid=ref_ids),
        # expand(f"{OUTPUT}epi2me_digest/{{cid}}.{{rid}}.bam", cid=cell_ids, rid=ref_ids),
        # expand(f"{OUTPUT}epi2me/{{cid}}.{{rid}}.ns.bam", cid=cell_ids, rid=ref_ids),
        # expand(f"{OUTPUT}align_table/{{cid}}.{{rid}}.{{cond}}.parquet", cid=cell_ids, rid=ref_ids, cond=['raw', 'digested']),
    

rule basecalling:
    input:
        expand(f"{OUTPUT}gtf/{{gid}}.gtf.gz", gid=gtf_ids),
        expand(f"{OUTPUT}references/{{rid}}.mmi", rid=ref_ids),
        expand(f"{OUTPUT}pod5/{{cid}}.pod5", cid=cell_ids),
        expand(f"{OUTPUT}fast5/{{cid}}.fast5", cid=cell_ids),
        expand(f"{OUTPUT}fastq/{{cid}}.raw.fastq", cid=cell_ids),
        expand(f"{OUTPUT}fastq/{{cid}}.raw.fastq.index", cid=cell_ids),
        expand(f"{OUTPUT}fastq/{{cid}}.digested.fastq", cid=cell_ids),


rule archive:
    input:
        expand(f"{OUTPUT}reports/sequence_reports/{{cid}}.report.pq", cid=cell_ids),
        expand(f"{OUTPUT}reports/nanostat/{{cid}}.{{rid}}.{{dig}}.tsv", cid=cell_ids, rid=ref_ids, dig=['raw', 'digested']),


rule haplotagging:
    input:
        expand(f"{OUTPUT}vcf/{{sid}}.vcf.gz.tbi", sid=snp_ids),
        expand(f"{OUTPUT}vcf/{{sid}}.snps.tsv", sid=snp_ids),
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
        expand(f"{OUTPUT}whatshap/{{sid}}.GRCm39.phased.vcf.gz", sid=snp_ids),
        expand(f"{OUTPUT}reports/whatshap/{{sid}}.vcf.stats.txt", sid=snp_ids),
        expand(f"{OUTPUT}nanopolish/{{cid}}.{{sid}}.raw.phased.bam", cid=cell_ids, sid=snp_ids),
        expand(f"{OUTPUT}longread/{{cid}}.{{sid}}.phased.bam", cid=cell_ids, sid=snp_ids),
        expand(f"{OUTPUT}vcf/{{sid}}.phased.vcf", sid=snp_ids),
        expand(f"{OUTPUT}reports/whatshap/{{sid}}.CAST_EiJ.stats.txt", sid=snp_ids),
        expand(f"{OUTPUT}reports/whatshap/{{sid}}.129S1_SvImJ.stats.txt", sid=snp_ids),
        expand(f"{OUTPUT}reports/whatshap/{{sid}}.CAST_EiJ.phased.stats.txt", sid=snp_ids),
        expand(f"{OUTPUT}reports/whatshap/{{sid}}.129S1_SvImJ.phased.stats.txt", sid=snp_ids),



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
    