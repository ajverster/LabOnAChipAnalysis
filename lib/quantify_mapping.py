
import pandas as pd
from collections import defaultdict
import csv
import sys

def count_sam(infile):
    counts = defaultdict(int)
    with open(infile) as f:
        for line in csv.reader(f, delimiter = "\t"):
            if line[2] == "*":
                continue
            counts[line[2]] += 1
    return counts


if __name__ == "__main__":
    infile = sys.argv[1]
    counts = count_sam(infile)
    df_out = pd.DataFrame(counts, index=["Counts"]).T
    df_out["percentage"] = df_out["Counts"] / df_out["Counts"].sum() * 100.
    print(df_out)
