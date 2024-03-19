import sys
import re
import pandas as pd
import pysam
import numpy as np


print(f"{pysam.__version__=}")
print(f"{pysam.__file__=}")

def get_align_type(align):
    """A function to return one of: 
        'primary', 'secondary', 'supplementary'
    based on read
    """
    align_type = 'primary'
    
    if align.is_secondary:
        align_type = 'secondary'
        
    if align.is_supplementary:
        align_type = 'supplementary'
    return align_type
    

def parse_read_name(align, sep):
    """A function to parse the read name """
    read_name = align.qname
    offset = 0
    is_digested = False
    fragment_index = -1

    if sep in read_name:
        read_name = read_name.split(sep)[0]
        offset = int(align.qname.split(sep)[2])
        fragment_index = int(align.qname.split(sep)[1])
        is_digested = True

    return read_name, fragment_index, is_digested, offset


def get_mean_alignment_base_qualities(align):
    """A function to get the alignment-span
    base qualities """
    if align.query_alignment_qualities is None:
        base_qualities = [-1]
    else:
        base_qualities = align.query_alignment_qualities
    
    return np.mean(base_qualities)
    

def bam_to_df(bampath, sep="_"):
    """A function to parse a bam file into 
    a structured table """
    res = []
    bam = pysam.AlignmentFile(bampath)
    for align in bam:

        # parse read name
        read_name, fragment_index, is_digested, offset = parse_read_name(align, sep)

        # check for the alignment tag
        if align.has_tag('AS'):
            as_tag = align.get_tag('AS')
        else:
            as_tag = 0        

        # get the alignment type
        align_type = get_align_type(align)

        # get base qualities
        mean_qual = get_mean_alignment_base_qualities(align)
    
        row = {
            'read_name' : read_name,
            'fragment_index' : fragment_index,
            'is_digested_fragment' : is_digested,
            'align_type' : align_type,
            'is_forward' : align.is_forward,
            'is_mapped' : align.is_mapped,
            'mean_align_base_quality' : mean_qual,
            'read_length' : align.query_length,
            'read_start' : align.query_alignment_start + offset,
            'read_end' : align.query_alignment_end + offset, 
            'chrom' : align.reference_name,
            'reference_start' : align.reference_start,
            'reference_end' : align.reference_end,
            'mapping_quality' : align.mapping_quality,
            'alignment_score' : as_tag,
        }

        res.append(row)
    return pd.DataFrame(res)


if __name__ == "__main__":
    bam_path = sys.argv[1]  
    sep = sys.argv[2]  
    outpath = sys.argv[3]

    # load the reads in
    df = bam_to_df(bam_path, sep=sep)
    df.to_parquet(outpath, index=False)
    
    
    

    
    


   

    
   