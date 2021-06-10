
from Bio import SeqIO
import gzip
import pandas as pd
from collections import defaultdict

def count_kmers(seq, k=1):
    counts = defaultdict(int)
    for j in range(0, len(seq)-k):
        if "N" in seq[j:j+k]:
            continue
        counts[seq[j:j+k]] += 1
    return counts

def load_seq_cat(infile):
    seq_cat = ""
    with gzip.open(infile,"rt") as f:
        for seq in SeqIO.parse(f, "fasta"):
            seq_cat += seq.seq
    return seq_cat

def convert_to_df(counts):
    df_out = pd.DataFrame(counts, index=["Counts"]).T
    df_out["percentage"] = df_out["Counts"] / df_out["Counts"].sum() * 100.
    return df_out


if __name__ == "__main__":
    infile_seqs = ["Listeria_monocytogenes_ATCC_19114.fasta.gz", "GCF_010374945.1_ASM1037494v1_genomic.fna.gz"]
    names = ["listeria", "Ecoli"]
    props = [2/3, 1/3]
    
    seqs = [load_seq_cat(infile) for infile in infile_seqs]
    counts = [count_kmers(s) for s in seqs]
    dfs = [convert_to_df(c) for c in counts]
    print(names[0])
    print(dfs[0])
    print(dfs[0].loc[["C","G"],:]["percentage"].sum())
    print(names[1])
    print(dfs[1])
    print(dfs[1].loc[["C","G"],:]["percentage"].sum())
    props = dfs[0]["percentage"] * props[0] + dfs[1]["percentage"] * props[1]
    print(props)
