import sys
import re
import os
import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq
from Bio import Restriction
from Bio import SeqIO
from rapidfuzz import fuzz
import rapidfuzz

def find_restriction_sites(read_seq, rb):
    """A function to find all restriction sites 
    in a sequence """
    search_results = rb.search(Seq(read_seq))
    sites = list(search_results.values())[0]
    n_sites = len(sites)
    if n_sites == 0:
        sites = [-1]

    site_str = ";".join([str(x) for x in sites])
    return site_str, n_sites


def find_fuzzy(seq, tgt, min_sim):
  """
  Finds all fuzzy matches of the target sequence within the main sequence recursively.

  Args:
      seq: The main sequence to search.
      tgt: The target sequence to find fuzzy matches for.
      min_sim: Minimum similarity threshold (0-100).

  Returns:
      A list of tuples containing the matching character, matched substring, start and end positions, and similarity score.
  """
  matches = []
  best_match = None

  # Find the best alignment using a sliding window and update best_match
  for i in range(len(seq) - len(tgt) + 1):
    sub = seq[i:i + len(tgt)]
    score = fuzz.partial_ratio(tgt, sub)
    if score >= min_sim and (not best_match or score > best_match[4]):
      best_match = (seq[i], sub, i, i + len(tgt) - 1, score)

  # If a match is found, record it and continue recursively
  if best_match:
    matches.append(best_match)
    rem_seq = seq[:best_match[2]] + seq[best_match[3] + 1:]
    matches.extend(find_fuzzy(rem_seq, tgt, min_sim))

  return matches


def find_barcodes(read_seq, barcode, min_sim=92):
    """A function to find barcode sequences in reads """
    barcode_len = len(barcode)

    matches = find_fuzzy(read_seq, barcode, min_sim)
    n_matches = len(matches)
    
    if len(matches) == 0:
        pos = [-1]
        scores = [-1]
    else:
        pos = []
        scores = []
        for i, match in enumerate(matches):
            _, _, start, _, score = match
            # correct alignment positions for striped barcodes
            offset = barcode_len * i            
            start = start + offset
            pos.append(start)
            scores.append(score)

    pos_str = ";".join([str(x) for x in pos])
    score_str = ";".join([str(round(x, 2)) for x in scores])
    return pos_str, score_str, n_matches
    

def get_sequence_report(fpath, rb, barcode, barcode_rc):
    """A function to return read-level information 
    from fastq files """
    report_df = []
    
    count = -1
    stop = 500
    for read in pysam.FastxFile(fpath):
        count += 1
        if count == stop:
            break

        # get read information
        read_seq = read.sequence
        read_seq_length = len(read_seq)
        quals = read.get_quality_array()

        # find the restrictin enzyme sites
        enzyme_pos, enzyme_matches = find_restriction_sites(read_seq, rb)

        # find forward barcodes
        fw_pos_str, fw_score_str, fw_n_matches = find_barcodes(read_seq, barcode)

        # find reverse complement barcodes
        rv_pos_str, rv_score_str, rv_n_matches = find_barcodes(read_seq, barcode_rc)

        # total BC
        total_bc_matches = fw_n_matches + rv_n_matches

        row = {
            'read_name' : read.name,
            'read_seq_length' : read_seq_length,
            'enzyme_matches' : enzyme_matches,
            'enzyme_pos' : enzyme_pos,
            'total_bc_matches' : total_bc_matches,
            'forward_bc_pos' : fw_pos_str,
            'forward_bc_score' : fw_score_str,
            'forward_bc_macthes' : fw_n_matches,
            'reverse_bc_pos' : rv_pos_str,
            'reverse_bc_score' : rv_score_str,
            'reverse_bc_macthes' : rv_n_matches,
            'mean_base_quality' : int(np.mean(quals)),
            'median_base_quality' : int(np.median(quals)),
            'min_base_quality' : np.min(quals),
            'max_base_quality' : np.max(quals),
        }
        report_df.append(row)
    return pd.DataFrame(report_df)


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
    
    df = get_sequence_report(fastq_path, rb, barcode, barcode_rc)
    df.to_parquet(outpath, index=False)
    

    

    

    
   