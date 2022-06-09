

from pathlib import Path
from Bio import SeqIO
import sys
import gzip
import re
import pandas as pd
import argparse
from collections import defaultdict

def hash_strain_names(indir):

    name_dict = defaultdict(str)
    for infile in indir.glob("*.fna.gz"):
        with gzip.open(infile,"rt") as f:
            seq = next(SeqIO.parse(f, "fasta"))
            strain_name = seq.description.replace(seq.id,"").replace("Listeria monocytogenes","").strip()
            strain_name = re.sub(" genome assembly.*$","", strain_name)
            strain_name = re.sub(" chromosome.*$","", strain_name)
            strain_name = re.sub(", complete genome.*$","", strain_name)
            name_dict[infile.name] = strain_name
    return name_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--indir_genomes')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-m', '--mode', default="metadata")
    args = parser.parse_args()
    indir = Path(args.indir_genomes)

    name_dict = hash_strain_names(indir)
    print(name_dict)
    if args.mode == "fastANI":
        # tasks
        # fastANI output renaming
        df_fastani = pd.read_csv(args.infile, header=None, sep="\t")
        df_fastani["strain_name"] = [name_dict[Path(f).name] for f in df_fastani[1]]
        df_fastani.to_csv(args.outfile, sep = "\t", index=False, header=False)
    elif args.mode == "metadata":
        # adds the names to the ncbi genomes metadata
        df_metadata = pd.read_csv(args.infile, sep=",")
        print(df_metadata)
        df_metadata["strain_from_genome"] = [name_dict[Path(x).name + "_genomic.fna.gz"] for x in df_metadata["GenBank FTP"]]
        df_metadata.to_csv(args.outfile, sep = ",", index=False, header=True)

    # TODO: phylogenetic tree renaming


