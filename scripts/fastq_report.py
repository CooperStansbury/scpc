import sys
import re
import os
import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq
from Bio import Restriction
from Bio import SeqIO

def parse_fastq(fpath, rb, barcode, barcode_rc):
    """a function to parse a fastq file """
    res = []
    for read in pysam.FastxFile(fpath):
        read_seq = read.sequence
        
        # look for restriction sites
        search_results = rb.search(Seq(read_seq))
        sites = list(search_results.values())[0]
        n_sites = len(sites)
        if n_sites == 0:
            sites = [-1]

        # barcode searching
        n_barcode_forward = read_seq.count(barcode)
        n_barcode_reverse_comp = read_seq.count(barcode_rc)

        forward_sites = [-1]
        if n_barcode_forward > 0:
            forward_sites = [m.start() for m in re.finditer(barcode, read_seq)]

        reverse_comp_sites = [-1]
        if n_barcode_reverse_comp > 0:
            reverse_comp_sites = [m.start() for m in re.finditer(barcode_rc, read_seq)]

        # get the base call qualities
        quals = read.get_quality_array()
        row = {
            'read_name' : read.name,
            'seq_length' : len(read_seq),
            'n_enzymes' : n_sites,
            'enzyme_sites' : ";".join([str(x) for x in sites]),
            'n_barcode_forward' : n_barcode_forward,
            'forward_sites' : ";".join([str(x) for x in forward_sites]),
            'n_barcode_reverse_comp' : n_barcode_reverse_comp,
            'reverse_comp_sites' : ";".join([str(x) for x in reverse_comp_sites]),
            'mean_base_quality' : int(np.mean(quals)),
            'median_base_quality' : int(np.median(quals)),
            'min_base_quality' : np.min(quals),
            'max_base_quality' : np.max(quals),
        }
        res.append(row)
    return pd.DataFrame(res)


if __name__ == "__main__":
    fastq_path = sys.argv[1]  
    barcode_path = sys.argv[2]  
    enzyme = str(sys.argv[3])
    outpath = sys.argv[4]  
    
    # set up restriction enzyme
    rb = Restriction.RestrictionBatch([enzyme])
    
    # set up get barcodes
    barcode_id = os.path.basename(fastq_path).split(".")[0]
    code_df = pd.read_csv(barcode_path)
    barcode = code_df[code_df['cell_id'] == barcode_id]['barcode'].values[0]
    barcode_rc = str(Seq(barcode).reverse_complement())
    
    df = parse_fastq(fastq_path, rb, barcode, barcode_rc)
    df.to_parquet(outpath, index=False)
    

    

    

    
   