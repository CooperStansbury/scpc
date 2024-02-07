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

    sites = []
    
    for record in records:
        rb = Restriction.RestrictionBatch([enzyme])
        search_results = rb.search(record.seq)
        cut_sites = list(search_results.values())[0]

        # store site information for later
        site_row = {
            'read_name' : record.name,
            'n_sites' : len(cut_sites),
            'sites' : ";".join(cut_sites),
        }

        sites.append(site_row)

        # only write records for reads with enzyme sites            
        if len(cut_sites) > 0:
            
            subread_counter = 0
            for i in range(0, len(cut_sites)):
                if i == 0:
                    subread = record[0 : cut_sites[i]]
                elif i < len(cut_sites) - 1:
                    subread = record[cut_sites[i] : cut_sites[i+1]]
                else:
                    subread = record[cut_sites[i] : ]

                if len(subread.seq) <= 4:
                    continue

                subread.id += "_" + str(subread_counter)
                subread.name += "_" + str(subread_counter)
                subread.description += "_" + str(subread_counter)
                subread_counter += 1
                output_handle.write(subread.format("fastq"))

            subread = record[cut_sites[i]-1 :]
            subread.id += "_" + str(subread_counter)
            subread.name += "_" + str(subread_counter)
            subread.description += "_" + str(subread_counter)

            # handle tiny artifacts
            if len(subread.seq) <= 4:
                continue
                
            output_handle.write(subread.format("fastq"))

    output_handle.close()

    sites = pd.DataFrame(sites)
    return sites
    


if __name__ == "__main__":
    fastq_path = sys.argv[1]
    enzyme = str(sys.argv[2])
    fastq_outpath = sys.argv[3]
    report_outpath = sys.argv[4]

    # extract reads and metadata
    sites = digestFastq(fastq_path, enzyme, fastq_outpath)
    sites.to_parquet(report_outpath, index=False)