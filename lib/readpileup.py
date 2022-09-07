import sys
import csv
from collections import defaultdict
import numpy as np
import pandas as pd
import logging
import argparse
from tqdm import tqdm
from pathlib import Path

csv.field_size_limit(sys.maxsize)


def read_pileup(infile):

    vals = defaultdict(list)
    genome = defaultdict(list)
    positions = defaultdict(list)

    with open(infile) as f:
        for line in csv.reader(f, delimiter = "\t"):
            vals[line[0]].append(int(line[3]))
            genome[line[0]].append(line[2])
            positions[line[0]].append(line[1])
    
    for sp in vals:
        vals[sp] = np.array(vals[sp])
        genome[sp] = np.array(genome[sp])
        positions[sp] = np.array(positions[sp])
    return vals, genome, positions

def quantify_qc_windows(vals, genome, positions, window_sizes=[1000], step_size = 50):
    results = []
    for sp in genome:
        for window in tqdm(window_sizes, desc="going through window sizes"):
            assert len(genome[sp]) == len(vals[sp])
            i = 0
            while i + window < len(genome[sp]):
                gc_content = np.isin(genome[sp][i:i+window], ["G","C"]).sum() / window
                covg = np.mean(vals[sp][i:i+window])
                # if covg > 1000:
                #    print(i, positions[sp][i], i+window, sp, covg)
                results.append([sp,positions[sp][i], gc_content, covg, window])
                i += step_size
    return results


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument('--infiles', default=[], action="append")
    parser.add_argument('--outfile')
    parser.add_argument('--window_sizes', default=[], type=int, action="append")
    parser.add_argument('--step_size', default=50, type=int)
    args = parser.parse_args()
    
    df_full = pd.DataFrame()
    for infile in args.infiles:
        logging.info("currently processing {}".format(infile))
        vals, genome, positions = read_pileup(infile)
        results = quantify_qc_windows(vals, genome, positions, args.window_sizes, args.step_size)

        df_intermediate = pd.DataFrame()
        df = pd.DataFrame(results, columns = ["species","pos","gc","covg","window_size"])
        df_intermediate = pd.concat([df_intermediate, df])
        df_intermediate["infile"] = Path(infile).name
        for key in vals:
            print("key: %s,  covg of bins: %.2f std of bin covg %.2f" %(key,  np.mean(vals[key]), np.std(vals[key])))

        df_full = df_full.append(df_intermediate)
    
    df_full.to_csv(args.outfile, index=False)

