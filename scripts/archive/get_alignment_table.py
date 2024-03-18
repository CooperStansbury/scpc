import sys
import re
import pandas as pd
import pysam


if __name__ == "__main__":
    bam_path = sys.argv[1]  
    restriction_sites = sys.argv[2]  
    outpath = sys.argv[3]

    # load files
    bamfile = pysam.AlignmentFile(bam_path)
    sites = pd.read_parquet(restriction_sites)

    sites['query_name'] = sites['read_name'] + "_" + sites['fragment_idx'].astype(str)

    digest_map = sites.copy()
    digest_map = digest_map.set_index('query_name')
    digest_map = digest_map[['read_start', 'read_end', 'seq']]

    digest_map = digest_map.to_dict(orient='index')
    del sites
    
    align_dict = {}
    alignment_records = []
    
    for align in bamfile:
        # resolve digested reads
        if "_" in align.query_name:
            read_name = align.query_name.split("_")[0]
            align_type = 'digested'
    
            # get the start and end position the the digested 
            # fragment
            digest_query = digest_map[align.query_name]
    
            # offset the digested fragment based on the mappable 
            # region
            as_start = align.query_alignment_start
            as_end = align.query_alignment_end
            read_start = digest_query['read_start'] + as_start
            read_end =  digest_query['read_start'] + as_end
        else:
            read_name = align.query_name
            align_type = 'full'
            read_start = align.query_alignment_start
            read_end = align.query_alignment_end
    
        # assign an id to each alignment
        if read_name in align_dict:
            align_dict[read_name] += 1
        else:
            align_dict[read_name] = 0
    
        align_idx = align_dict[read_name]
        
        # check for the alignment tag
        if align.has_tag('AS'):
            as_tag = align.get_tag('AS')
        else:
            as_tag = 0
    
        align_row = {
                'read_name' : read_name,
                'align_idx' : align_idx,
                'align_type' : align_type,
                'read_start' : read_start,
                'read_end' : read_end,
                'chrom' : align.reference_name,
                'ref_start' : align.reference_start,
                'ref_end' : align.reference_end,
                'mapping_quality' : align.mapping_quality,
                'is_reverse' : align.is_reverse,
                'is_duplicate' : align.is_duplicate,
                'alignment_score' : as_tag,
            }
        
        alignment_records.append(align_row)
        # end for
        
    align = pd.DataFrame(alignment_records)
    align = align.sort_values(by=['read_name', 'read_start'])
    align.to_parquet(outpath, index=False)
    

    
   