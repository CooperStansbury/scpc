import sys
import re
import os
import numpy as np
import pandas as pd
import pysam


if __name__ == "__main__":
    barcode_path = sys.argv[1]  
    outpath = sys.argv[2]  
    
    df = pd.read_csv(barcode_path)

    # remove if exists
    if os.path.exists(outpath):
      os.remove(outpath)

    with open(outpath, 'w') as outfile:
        for idx, row in df.iterrows():
            print(f">{row['cell_id']}", file=outfile)
            print(f"{row['barcode']}", file=outfile)

            
        
    

    

    

    
   