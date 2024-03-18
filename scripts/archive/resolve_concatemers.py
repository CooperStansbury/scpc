import sys
import re
import pandas as pd
import networkx as nx
from itertools import combinations
    
def get_concatemers(df):
    """A function to resolve concatemers 
    from multiple alignments """
    concatemers = []

    for read_name, align_group in df.groupby('read_name'):
        column_sort = [
            'read_start', 
            'length_on_read', 
            'mapping_quality'
        ]

        sort_order = [
            True, 
            False, 
            False
        ]
        
        # sort nodes by the starting position on the read
        align_group = align_group.sort_values(by=column_sort, 
                                              ascending=sort_order)
        # get the longest, highest quality alignments
        align_group = align_group.drop_duplicates(subset='read_start')

        # find overlapping alignments
        intervals = []
        for rs, re in align_group[['read_start', 'read_end']].values:
            interval = pd.Interval(rs, re)
            intervals.append(interval)

        # drop the lower quality overlapping alignments
        to_drop = []
        for i, align_1 in enumerate(intervals):
            for j, align_2 in enumerate(intervals):
                if i == j:
                    continue
                    
                if align_1.overlaps(align_2):
                    overlaps = align_group.take([i, j])
                    overlaps = overlaps.sort_values(by=column_sort, 
                                                    ascending=sort_order)
                    to_drop.extend(overlaps['alignment_idx'].values[1:])

        # drop unique only
        to_drop = list(set(to_drop))
        align_group = align_group[~align_group['alignment_idx'].isin(to_drop)]

        if len(align_group) > 1:
            concatemers.append(align_group)

    if len(concatemers) == 0:
        concatemers = pd.DataFrame(columns = df.columns)
    else:
        concatemers = pd.concat(concatemers)
        
    return concatemers


if __name__ == "__main__":
    align_path = sys.argv[1]  
    mapq_threshold = sys.argv[2]
    outpath = sys.argv[3]

    # define params
    mapq_threshold = float(mapq_threshold)

    # load data
    align = pd.read_parquet(align_path)

    # filter low-quality alignments
    align = align[align['mapping_quality'] > mapq_threshold]

    # only consider reads with multiple alignments
    align['n_alignments'] = align.groupby('read_name')['alignment_idx'].transform('nunique')
    align = align[align['n_alignments'] > 1]

    # add alignment length
    align['length_on_read'] = align['read_end'] - align['read_start']

    concatemers = get_concatemers(align)
    concatemers['length_on_ref'] = concatemers['ref_end'] - concatemers['ref_start']
    concatemers['midpoint'] = concatemers['ref_start'] + (concatemers['length_on_ref'] / 2 ) 
    concatemers['midpoint'].astype(int)

    print(concatemers.dtypes)
    concatemers.to_parquet(outpath, index=False)


    
    
   
    

    

    
   