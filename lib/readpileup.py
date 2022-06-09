import sys
import csv
from collections import defaultdict
import numpy as np
import pandas as pd
import logging
from pathlib import Path

csv.field_size_limit(sys.maxsize)


def read_pileup(infile):

    vals = defaultdict(list)
    genome = defaultdict(list)

    with open(infile) as f:
        for line in csv.reader(f, delimiter = "\t"):
            vals[line[0]].append(int(line[3]))
            genome[line[0]].append(line[2])
    
    for sp in vals:
        vals[sp] = np.array(vals[sp])
        genome[sp] = np.array(genome[sp])
    return vals, genome

def quantify_qc_windows(vals, genome, window_size=1000, step_size = 50):
    results = defaultdict(list)
    for sp in genome:
        assert len(genome[sp]) == len(vals[sp])
        i = 0
        while i + window_size < len(genome[sp]):
            gc_content = np.isin(genome[sp][i:i+window_size], ["G","C"]).sum() / window_size
            covg = np.mean(vals[sp][i:i+window_size])
            results[sp].append((i, gc_content, covg))
            i += step_size
    return results


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    infiles = sys.argv[1:-1]
    outfile = sys.argv[-1]
    
    df_full = pd.DataFrame()
    for infile in infiles:
        logging.info("currently processing {}".format(infile))
        vals, genome = read_pileup(infile)
        results = quantify_qc_windows(vals, genome)

        df_intermediate = pd.DataFrame()
        for sp in results:
            df = pd.DataFrame(results[sp], columns = ["pos","gc","covg"])
            df["species"] = sp
            df_intermediate = df_intermediate.append(df)
        df_intermediate["infile"] = Path(infile).parts[0]
        df_intermediate.to_csv(outfile, index=False)
        
        print(infile)
        for key in vals:
            print("key: %s,  covg of bins: %.2f std of bin covg %.2f" %(key,  np.mean(vals[key]), np.std(vals[key])))

        df_full = df_full.append(df_intermediate)
    df_full.to_csv(outfile, index=False)

