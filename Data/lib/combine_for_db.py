
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import re
import pandas as pd
from collections import defaultdict
import argparse

def load_seq_cat(infile, mode="cat"):
    seq_cat = ""
    with open(infile,"rt") as f:
        for seq in SeqIO.parse(f, "fasta"):
            seq_cat += seq.seq
            if args.mode == "first":
                print("Using only the first sequence")
                break
    assert len(seq_cat) > 0
    return seq_cat


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-i', '--infile_genomes', action='append', default=[])
    parser.add_argument('-n', '--names', action='append', default=[])
    parser.add_argument('--mode',default="cat")
    args = parser.parse_args()
        
    seqs = [load_seq_cat(infile, args.mode) for infile in args.infile_genomes]
    seqs = [SeqRecord(s, id=args.names[i], name="", description="") for i, s in enumerate(seqs)]
    with open(args.outfile, "w") as f:
        SeqIO.write(seqs, f, "fasta")


