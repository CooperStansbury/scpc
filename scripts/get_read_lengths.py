from threading import Thread
import queue
import sys
import pandas as pd
import pysam


class FastqProcessor(Thread):
  def __init__(self, fastq_file, queue):
    super().__init__()
    self.fastq_file = fastq_file
    self.queue = queue

  def run(self):
    with pysam.FastxFile(self.fastq_file, "r") as f:
      for read in f:
        self.queue.put(read)


def process_read(read):
    """A function to return sequence stats """
    record = {
        'read_name' : read.name,
        'read_length' : len(read.sequence),
    }
    return record



if __name__ == "__main__":
    fastq_path = sys.argv[1]  
    num_threads = int(sys.argv[2])
    outpath = sys.argv[3]

    # build the queue
    queue = queue.Queue()
    
    processor = FastqProcessor(fastq_path, queue)
    processor.start()

  # Adjust num_threads based on your hardware and file size
    threads = []
    records = []
    for _ in range(num_threads):
        thread = Thread(target=lambda: records.append(process_read(queue.get())))
        thread.start()
        threads.append(thread)
    
    # Wait for threads to finish
    for thread in threads:
        thread.join()
    
    processor.join()  # Wait for processor thread to finish

    df = pd.DataFrame(records)
    df.to_parquet(outpath, index=False)
    


    