import sys
import csv
from collections import defaultdict
import numpy as np
import pandas as pd
import logging
import argparse
from tqdm import tqdm
from pathlib import Path
from Bio import SeqIO

csv.field_size_limit(sys.maxsize)


def read_lengths(infile):
    lens = {}
    with open(infile) as f:
        for seq in SeqIO.parse(f, "fasta"):
            lens[seq.id] = len(seq)
    return lens

def read_pileup(infile, len_dict):
    rows = []
    with open(infile) as f:
        pos_last = 0
        chrom_last = ""
        for line in tqdm(csv.reader(f, delimiter = "\t"),desc="reading pileup"):
            chrom = line[0]
            pos = int(line[1])
            ref = line[2]
            if (chrom != chrom_last) & (chrom_last != ""):
                # Fill zeros until the end of the chromosome
                while pos_last < len_dict[chrom_last]:
                    rows.append([chrom_last, ref, pos_last, 0])
                    pos_last += 1
            if chrom != chrom_last:
                pos_last = 0
                chrom_last = chrom
            while pos != pos_last + 1:
                # Add zeros
                pos_last += 1
                rows.append([chrom, ref, pos_last, 0])
            rows.append([chrom, ref, pos, int(line[3])])
            pos_last = pos
    return pd.DataFrame(rows, columns=["chromosome","ref","pos","covg"])


def quantify_qc_windows(df_covg, len_dict, window_sizes=[1000], step_size = 50):
    results = []
    for sp in df_covg["chromosome"].unique():
        df_covg_sub = df_covg.query('chromosome == "{}"'.format(sp))
        df_covg_sub.index = df_covg_sub["pos"].values
        for window in tqdm(window_sizes, desc="going through window sizes"):
            i = 1
            with tqdm(total=len_dict[sp], desc="quantifying {} window size {}".format(sp, window)) as pbar:
                while i + window < len_dict[sp]:
                    pos = i + window/2
                    #df_covg_sub_sub = df_covg_sub.query('(pos >= {}) & (pos < {})'.format(i, i+window))
                    df_covg_sub_sub = df_covg_sub.loc[i:i+window,:]
                    gc = df_covg_sub_sub["ref"].isin(["G","C"]).sum() / window
                    covg = df_covg_sub_sub["covg"].mean()
                    results.append([sp, pos, gc, covg, window])
                    i += step_size
                    pbar.update(step_size)
    return results


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument('--infiles', default=[], action="append")
    parser.add_argument('--outfile')
    parser.add_argument('--infile_seqs')
    parser.add_argument('--window_sizes', default=[], type=int, action="append")
    parser.add_argument('--step_size', default=50, type=int)
    args = parser.parse_args()
   
    len_dict = read_lengths(args.infile_seqs)
 
    df_full = pd.DataFrame()
    for infile in args.infiles:
        logging.info("currently processing {}".format(infile))
        #vals, genome, positions = read_pileup(infile)
        df_covg = read_pileup(infile, len_dict)
        print(df_covg)
        results = quantify_qc_windows(df_covg, len_dict, window_sizes=args.window_sizes)

        df_intermediate = pd.DataFrame()
        df = pd.DataFrame(results, columns = ["species","pos","gc","covg","window_size"])
        df_intermediate = pd.concat([df_intermediate, df])
        df_intermediate["infile"] = Path(infile).name
        #for key in vals:
        #    print("key: %s,  covg of bins: %.2f std of bin covg %.2f" %(key,  np.mean(vals[key]), np.std(vals[key])))

        df_full = df_full.append(df_intermediate)
    
    df_full.to_csv(args.outfile, index=False)

