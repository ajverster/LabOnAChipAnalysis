
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import pandas as pd
from collections import defaultdict

def load_seq_cat(infile):
    seq_cat = ""
    with gzip.open(infile,"rt") as f:
        for seq in SeqIO.parse(f, "fasta"):
            seq_cat += seq.seq
    return seq_cat


if __name__ == "__main__":
    outfile = "combined_ecoli_listeria.fasta"
    infile_seqs = ["Listeria_monocytogenes_ATCC_19114.fasta.gz", "GCF_010374945.1_ASM1037494v1_genomic.fna.gz"]
    names = ["listeria", "Ecoli"]
    props = [2/3, 1/3]
    
    seqs = [load_seq_cat(infile) for infile in infile_seqs]
    seqs = [SeqRecord(s, id=names[i], name="", description="") for i, s in enumerate(seqs)]
    with open(outfile, "w") as f:
        SeqIO.write(seqs, f, "fasta")


