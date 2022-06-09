
from Bio import SeqIO
from pathlib import Path
import argparse
import shutil
import subprocess

def get_species_seq(infile, species_oi):
    seqs_oi = []
    with open(infile) as f:
        for seq in SeqIO.parse(f, "fasta"):
            if seq.id == species_oi:
                seq.id = infile.name.split("_")[0]
                seqs_oi.append(seq)
    assert len(seqs_oi) == 1, "Couldn't find your species_oi in the consensus file"
    return seqs_oi[0]

def create_genome_links(indir_genomes, outdir):
    for infile in Path(indir_genomes).glob("*.fna.gz"):
        outfile = outdir / infile.name
        outfile.symlink_to(infile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--infile_consensus', action='append', default=[])
    parser.add_argument('-s', '--species')
    parser.add_argument('-g', '--indir_genomes')
    parser.add_argument('-o', '--outdir_use')
    args = parser.parse_args()

    indir_genomes = Path(args.indir_genomes)
    outdir_use = Path(args.outdir_use)
    shutil.rmtree(outdir_use)
    outdir_use.mkdir(exist_ok=True)


    create_genome_links(indir_genomes, outdir_use)

    for infile in args.infile_consensus:
        infile = Path(infile)
        seq_consensus = get_species_seq(infile, args.species)
        outfile = outdir_use / infile.name
        outfile = str(outfile).replace(".fasta",".fna")
        with open(outfile,"w") as f:
            SeqIO.write(seq_consensus, f, "fasta")

        # Gzip output
        subprocess.call(["gzip",outfile])
