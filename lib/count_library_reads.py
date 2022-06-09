
import sys
import subprocess
from pathlib import Path
import pandas as pd

if __name__ == "__main__":
    infiles = sys.argv[1:-1]
    outfile = sys.argv[-1]

    r=[]
    for inf in infiles:
        ps = subprocess.Popen(["zcat", inf], stdout=subprocess.PIPE)
        output = int(subprocess.check_output(["wc", "-l"], stdin=ps.stdout).split()[0]) / 4
        r.append([Path(inf).name, output])
    df_out = pd.DataFrame(r, columns = ["file","n_reads"])
    df_out.to_csv(outfile, sep = ",",index=False)
