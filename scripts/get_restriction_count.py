from threading import Thread
import queue
import sys
import pandas as pd
import pysam
from Bio import Restriction
from Bio.Seq import Seq


class FastqProcessor(Thread):
  def __init__(self, fastq_file, queue):
    super().__init__()
    self.fastq_file = fastq_file
    self.queue = queue

  def run(self):
    with pysam.FastxFile(self.fastq_file, "r") as f:
      for read in f:
        self.queue.put(read)


def find_restriction_sites(read_seq, rb):
    """A function to find all restriction sites 
    in a sequence """
    search_results = rb.search(Seq(read_seq))
    sites = list(search_results.values())[0]
    n_sites = len(sites)
    return n_sites


def process_read(read, rb):
    """A function to return sequence stats """
    record = {
        'read_name' : read.name,
        'n_enzyme' : find_restriction_sites(read.sequence, rb),
    }
    return record


if __name__ == "__main__":
    fastq_path = sys.argv[1]  
    num_threads = int(sys.argv[2])
    enzyme = str(sys.argv[3])
    outpath = sys.argv[4]

    # set up restriction enzyme
    rb = Restriction.RestrictionBatch([enzyme])

    # build the queue
    queue = queue.Queue()
    
    processor = FastqProcessor(fastq_path, queue)
    processor.start()

    # Adjust num_threads based on your hardware and file size
    threads = []
    records = []
    for _ in range(num_threads):
        thread = Thread(target=lambda: records.append(process_read(queue.get(), rb)))
        thread.start()
        threads.append(thread)
    
    # Wait for threads to finish
    for thread in threads:
        thread.join()
    
    processor.join()  # Wait for processor thread to finish

    df = pd.DataFrame(records)
    df.to_parquet(outpath, index=False)
    


    