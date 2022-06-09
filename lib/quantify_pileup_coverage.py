import numpy as np
import pandas as pd
import csv
import sys
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path


csv.field_size_limit(sys.maxsize)

def quantify_sequence_length(infile_seqs):
    len_dict = {}
    with open(infile_seqs) as f:
        for seq in SeqIO.parse(f, "fasta"):
            len_dict[seq.id] = len(seq)
    return len_dict


def quantify_pileup(infile, len_dict):
    data = defaultdict(list)
    i = 1
    seq_current = None
    with open(infile) as f:
        for line in csv.reader(f, delimiter="\t"):
            seq = line[0]
            pos = int(line[1])
            covg = int(line[3])
            # At the end, check to see if we've covered the whole genome, if not, add some zeros
            if (seq != seq_current) & (seq_current is not None):
                if i < len_dict[seq_current]:
                    for j in range(i, len_dict[seq_current]):
                        data[seq].append(0)
                i = 1

            if pos != i:
                # Need to add some zeros
                for j in range(i, pos):
                    data[seq].append(0)
            data[seq].append(covg)

            i = pos + 1
            seq_current = seq
    return data


if __name__ == "__main__":
    infile_seqs = sys.argv[1]
    infile_pileups = sys.argv[2:-1]
    outfile = sys.argv[-1]
    covg_cutoff_test = [10,30,100]
    len_dict = quantify_sequence_length(infile_seqs)
    lines = []
    for infile in infile_pileups:
        data = quantify_pileup(infile, len_dict)
        for seq in data:
            r = [Path(infile).name, seq, np.mean(data[seq])]
            for covg_cutoff in covg_cutoff_test:
                r.append(np.sum(np.array(data[seq]) >= covg_cutoff) / len(data[seq]))
            lines.append(r)
    df_out = pd.DataFrame(lines, columns = ["infile","seq","mean_covg"] + ["covg_{}".format(i) for i in covg_cutoff_test])
    df_out.to_csv(outfile, index=False)
