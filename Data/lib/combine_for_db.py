
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import pandas as pd
from collections import defaultdict
import argparse

def load_seq_cat(infile):
    seq_cat = ""
    with gzip.open(infile,"rt") as f:
        for seq in SeqIO.parse(f, "fasta"):
            seq_cat += seq.seq
    return seq_cat


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-i', '--infile_genomes', action='append', default=[])
    parser.add_argument('-n', '--names', action='append', default=[])
    args = parser.parse_args()
        
    #infile_seqs = ["Listeria_monocytogenes_ATCC_19114.fasta.gz", "GCF_010374945.1_ASM1037494v1_genomic.fna.gz"]
    #names = ["listeria", "Ecoli"]
    #props = [2/3, 1/3]
    
    seqs = [load_seq_cat(infile) for infile in args.infile_genomes]
    seqs = [SeqRecord(s, id=args.names[i], name="", description="") for i, s in enumerate(seqs)]
    with open(args.outfile, "w") as f:
        SeqIO.write(seqs, f, "fasta")


