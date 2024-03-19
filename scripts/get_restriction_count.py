from threading import Thread
import queue
import sys
import pandas as pd
import pysam
from Bio import Restriction
from Bio.Seq import Seq


def find_restriction_sites(read_seq, rb):
    """A function to find all restriction sites 
    in a sequence """
    search_results = rb.search(Seq(read_seq))
    sites = list(search_results.values())[0]
    n_sites = len(sites)
    return n_sites


if __name__ == "__main__":
    fastq_path = sys.argv[1]  
    enzyme = str(sys.argv[2])
    outpath = sys.argv[3]

    # set up restriction enzyme
    rb = Restriction.RestrictionBatch([enzyme])

    records = []
    with pysam.FastxFile(fastq_path, "r") as f:
        for read in f:
            record = {
                'read_name' : read.name,
                'n_enzymes' : find_restriction_sites(read.sequence, rb),
            }

            records.append(record)
          
    df = pd.DataFrame(records)
    df.to_parquet(outpath, index=False)
    


    