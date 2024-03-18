import sys
import re
import pandas as pd
from Bio import Restriction
from Bio.Seq import Seq
from Bio import SeqIO


if __name__ == "__main__":
    fastq_path = sys.argv[1]  
    enzyme = str(sys.argv[2])
    outpath = sys.argv[3]

    handle = open(fastq_path, "r")
    records = list(SeqIO.parse(handle, "fastq"))
    handle.close()
    
    site_df = []
    
    for read in records:
        rb = Restriction.RestrictionBatch([enzyme])
        search_results = rb.search(read.seq)
        sites = list(search_results.values())[0]
        n_sites = len(sites)
        sites.insert(0, 0) # insert the beginning of the read
        sites.append(len(read.seq)) # append the end of the read

        for i, site in enumerate(sites[:-1]):
            fragment_row = {
                'read_name' : read.name,
                'n_sites' : n_sites,
                'fragment_idx' : i,
                'read_start' : site,
                'read_end' : sites[i+1],
                'length' : sites[i+1] - site,
                'seq' : str(read.seq[site:sites[i+1]]),
            }

            site_df.append(fragment_row)
        
    site_df = pd.DataFrame(site_df)
    site_df.to_parquet(outpath, index=False)

    

    

    
   