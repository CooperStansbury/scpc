import sys
import re
import pandas as pd
import networkx as nx
from itertools import combinations
    
def clique_expand(df):
    """A function to clique-expand reads """
    contacts = []

    for read_name, read_df in df.groupby("read_name"):
        if len(read_df) <= 1:
            continue

        rows = list(
            read_df.sort_values(["read_start"], ascending=True)
            .itertuples()
        )

        for (align_1, align_2) in combinations(rows, 2):
            contact_row = {
                'read_name' : align_1.read_name,
                'align1_idx' : align_1.alignment_idx,
                'align1_chrom' : align_1.chrom,
                'align1_ref_start' : align_1.ref_start,
                'align1_ref_end' : align_1.ref_end,
                'align1_read_start' : align_1.read_start,
                'align1_ref_midpoint' : align_1.midpoint,
                'align1_read_end' : align_1.read_end,
                'align1_mapping_quality' : align_1.mapping_quality,
                'align2_chrom' : align_2.chrom,
                'align2_ref_start' : align_2.ref_start,
                'align2_ref_end' : align_2.ref_end,
                'align2_ref_midpoint' : align_2.midpoint,
                'align2_read_start' : align_2.read_start,
                'align2_read_end' : align_2.read_end,
                'align2_mapping_quality' : align_2.mapping_quality,
            }

            contacts.append(contact_row)
    
    contacts = pd.DataFrame(contacts)
    return contacts
    


  

if __name__ == "__main__":
    concatemer_path = sys.argv[1]  
    outpath = sys.argv[2]

    # load data
    concatemer = pd.read_parquet(concatemer_path)
    print(concatemer.head())
    print()
    contacts = clique_expand(concatemer)
    print(contacts)
    



    
    
   
    

    

    
   