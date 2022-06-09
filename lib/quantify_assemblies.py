
from Bio import SeqIO
import numpy as np
import pandas as pd
import sys
from pathlib import Path

def count_sizes(infile):
    sizes = []
    with open(infile) as f:
        for seq in SeqIO.parse(f, "fasta"):
            sizes.append(len(seq))
    return sizes

if __name__ == "__main__":
    infiles = sys.argv[1:-1]
    outdir = Path(sys.argv[-1])
    df_full = pd.DataFrame()
    for infile in infiles:
        sizes = count_sizes(infile)
        total = np.sum(sizes)
        sizes_cumsum = np.cumsum(sizes)
        upper_half = sizes[:np.where(sizes_cumsum >= total/2)[0][0]]
        n50 = upper_half[-1]
        l50 = len(upper_half)
        df_out = pd.DataFrame([["n50",n50],["l50",l50]], columns = ["metric","value"])
        df_out["infile"] = Path(infile).parts[-3]
        df_full = df_full.append(df_out)
    df_full.to_csv(outdir / "assembly_stats.csv", index=False)
