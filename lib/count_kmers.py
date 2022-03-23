
import numpy as np
import pandas as pd
import gzip
import sys
import csv
import tqdm
from collections import defaultdict

def hash_reads(infile):
    if infile is None:
        return None
    reads_oi = {}
    with open(infile, "tr") as f:
        for line in csv.reader(f, delimiter = "\t"):
            read = "@" + line[0].replace("_"," ")
            reads_oi[read] = 1
    return reads_oi

def count(infile, k = 1, reads_oi = None):
    counts = defaultdict(int)
    with gzip.open(infile, "rb") as f:
        for i, line in tqdm.tqdm(enumerate(f)):
            if i % 4 == 0:
                name = line.decode("utf-8").rstrip("\n")
            if reads_oi is not None:
                if name not in reads_oi:
                    continue
            if i % 4 == 1:
                line = line.decode("utf-8")
                for j in range(0, len(line)-k):
                    counts[line[j:j+k]] += 1
    return counts

if __name__ == "__main__":
    # Fastq.gz sequencing file
    infile = sys.argv[1]
    # CSV file you want to read
    outfile = sys.argv[2]
    # If you want to only count GC in a subset of the reads. Currently only work on the metaphlan 
    if len(sys.argv) > 3:
        infile_reads_oi = sys.argv[3]
    else:
        infile_reads_oi = None
    reads_oi = hash_reads(infile_reads_oi)

    counts = count(infile, reads_oi=reads_oi)
    df_out = pd.DataFrame(counts, index=["Counts"]).T
    df_out["percentage"] = df_out["Counts"] / df_out["Counts"].sum() * 100.
    df_out.to_csv(outfile)

