import sys
import re
import pandas as pd
import pysam


if __name__ == "__main__":
    vcf_path = sys.argv[1]  
    bam_path = sys.argv[2]  
    restriction_sites = sys.argv[3]  
    outpath = sys.argv[4]

    # load the snps and aligments
    vcf_file = pysam.VariantFile(vcf_path)
    bam_file = pysam.AlignmentFile(bam_path)

    # load the restriction sites
    sites = pd.read_parquet(restriction_sites)
    
    sites['query_name'] = sites['read_name'] + "_" + sites['fragment_idx'].astype(str)
    
    digest_map = sites.copy()
    digest_map = digest_map.set_index('query_name')
    digest_map = digest_map[['read_start', 'read_end', 'seq']]
    
    digest_map = digest_map.to_dict(orient='index')
    del sites


    # loop through all reads
    align_dict = {}
    res = []
    samples = []
    for align in bam_file:
            
        # resolve digested reads
        if "_" in align.query_name:
            read_name = align.query_name.split("_")[0]
            align_type = 'digested'
        else:
            read_name = align.query_name
            align_type = 'full'

        # extract the reference positions of the read
        read_start = align.reference_start
        read_end = align.reference_end
        chrom = align.reference_name 

        # assign an id to each alignment
        if read_name in align_dict:
            align_dict[read_name] += 1
        else:
            align_dict[read_name] = 0
    
        align_idx = align_dict[read_name]

        # if the chrom not in the VCF, skip the read
        if chrom is None:
            continue
        
        if not vcf_file.is_valid_reference_name(chrom):
            continue

        # get the sequence, skip if None
        seq = align.get_forward_sequence()
        if seq is None:
            continue

        # attempt to extract the basecalling qualities,
        # assume zero if not assigned
        try:
            quals = align.get_forward_qualities()
        except:
            quals = [0] * len(seq)

        # extract the read mapping quality
        mapq = align.mapping_quality
        if not mapq > 0:
            continue
        
        # loop through all snps
        for vcf_rec in vcf_file.fetch(chrom, read_start, read_end):    
            
            idx = vcf_rec.pos - read_start
            if idx > len(seq) - 1:
                continue
            if idx < 0:
                continue
    
            read_base = seq[idx]
            read_qual = quals[idx]
        
            snp_row = {
                'read_name' : read_name,
                'read_length_on_ref' : read_end - read_start,
                'align_idx' : align_idx,
                'reference_base' : vcf_rec.ref,
                'read_base' : read_base,
                'read_base_quality' : read_qual,
                'read_mapping_quality' : mapq,
            }
            
            lookup = dict(vcf_rec.samples)
    
            for k, v in lookup.items():
                snp_row[k] = v.alleles[0]

                # for assigning metadata columns later 
                if not k in samples:
                    samples.append(k)
    
            res.append(snp_row)
            
    # compile 
    res = pd.DataFrame(res)

    # add metadata columns 
    res['read_match_ref'] = res['read_base'] == res['reference_base']
    for sample in samples:
        res[f'read_match_{sample}'] = res['read_base'] == res[sample]
    
    # save the file
    res.to_parquet(outpath, index=False)

