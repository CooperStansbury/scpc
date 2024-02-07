import sys
import re
import pandas as pd
import pysam


if __name__ == "__main__":
    fastq_path = sys.argv[1]  
    cutsite = sys.argv[2]
    site_outpath = sys.argv[3]
    ligated_outpath = sys.argv[4]

    fastq = pysam.FastxFile(fastq_path)

    site_df = []
    for read in fastq:

        sites = [m.start() for m in re.finditer(cutsite, read.sequence)]

        if len(sites) == 0:
            site_str = ""
        else:
            site_str = ";".join([str(x) for x in sites])

        
        read_row = {
            'read_name' : read.name,
            'n_sites' : len(sites),
            'sites' : site_str,
        }

        site_df.append(read_row)

    site_df = pd.DataFrame(site_df)
    site_df.to_parquet(site_outpath, index=False)

    # split reads with more than 0 cut sites
    ligated_reads = site_df[site_df['n_sites'] > 0].copy()
    ligated_reads = ligated_reads['read_name'].drop_duplicates()
    ligated_reads.to_csv(ligated_outpath, index=False, header=False)
    
    

    

    
   