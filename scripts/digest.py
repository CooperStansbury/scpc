import pandas as pd
import sys
from Bio import Restriction
from Bio.Seq import Seq
from Bio import SeqIO


def digestFastq(fastqPath, enzyme, outPath):
    """A function to digest a read based on a restriction site """
    output_handle = open(outPath, "w")
    
    handle = open(fastqPath, "r")
    records = list(SeqIO.parse(handle, "fastq"))
    handle.close()
    
    for record in records:
        rb = Restriction.RestrictionBatch([enzyme])
        search_results = rb.search(record.seq)
        cut_sites = list(search_results.values())[0]
        n_sites = len(cut_sites)

        cut_sites.insert(0, 0) # insert the beginning of the read
        cut_sites.append(len(record.seq)) # append the end of the read
        
        # only write records for reads with enzyme sites            
        if n_sites > 0:

            # loop through digested fragments
            for i, s_idx in enumerate(cut_sites[:-1]):
                subread = record[s_idx : cut_sites[i+1]]

                subread.id += "_" + str(i) + "_" + str(s_idx)
                subread.name += "_" + str(i) + "_" + str(s_idx) # add the start site to the read name
                subread.description += "_" + str(i) + "_" + str(s_idx)
                output_handle.write(subread.format("fastq"))
                
            output_handle.write(subread.format("fastq"))

    output_handle.close()


if __name__ == "__main__":
    fastq_path = sys.argv[1]
    enzyme = str(sys.argv[2])
    fastq_outpath = sys.argv[3]

    # extract reads 
    digestFastq(fastq_path, enzyme, fastq_outpath)
