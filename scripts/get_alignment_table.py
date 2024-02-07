import sys
import re
import pandas as pd
import pysam
from multiprocessing import Pool


if __name__ == "__main__":
    bam_path = sys.argv[1]  
    outpath = sys.argv[2]
    
    bamfile = pysam.AlignmentFile(bam_path)
    
    alignment_index = {}
    alignment_records = []
    for align in bamfile:
    
        # add an identifier for each alignment    
        if align.query_name in alignment_index:
            alignment_index[align.query_name] += 1
        else:
            alignment_index[align.query_name] = 0

        # get alignment score where possible
        if align.has_tag('AS'):
            as_tag = align.get_tag('AS')
        else:
            as_tag = 0

        # structure the alignments as a table
        align_row = {
            'read_name' : align.query_name,
            'alignment_idx' : alignment_index[align.query_name],
            'read_start' : align.query_alignment_start,
            'read_end' : align.query_alignment_end,
            'chrom' : align.reference_name,
            'ref_start' : align.reference_start,
            'ref_end' : align.reference_end,
            'mapping_quality' : align.mapping_quality,
            'is_reverse' : align.is_reverse,
            'is_duplicate' : align.is_duplicate,
            'alignment_score' : as_tag,
        }
        alignment_records.append(align_row)
    
    align_df = pd.DataFrame(alignment_records)
    align_df = align_df.sort_values(by=['read_name', 'alignment_idx'])
    align_df.to_parquet(outpath, index=False)
    

   
    

    

    
   