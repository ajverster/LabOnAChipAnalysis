import ete3
import pandas as pd
from collections import defaultdict
import argparse
import os
os.environ['QT_QPA_PLATFORM']='offscreen'


def update_meta(infile, name_dict):
    df_meta = pd.read_csv(infile, sep = ",").fillna("").query('Level == "Complete"')

    dict_info_oi = defaultdict(list)
    for (i, dat) in df_meta.iterrows():
        if ("4b" in dat["#Organism Name"]) | ("4b" in dat["Strain"]):
            dict_info_oi["4b"].append(dat["Assembly"])
        if dat["Strain"] not in dat["#Organism Name"]:
            name_dict[dat["Assembly"]] = (dat["#Organism Name"] + " " + dat["Strain"]).replace("Listeria monocytogenes ", "")
        else:
            name_dict[dat["Assembly"]] = dat["#Organism Name"].replace("Listeria monocytogenes ", "")
    return name_dict, dict_info_oi

def load_tree(infile_tree):
    with open(infile_tree) as f:
        tree_str = f.read().strip("\n")

    t = ete3.Tree(tree_str)
    return t

def plot_tree(t, name_dict, dict_info_oi, outfile):
    # Try re-rooting the tree to make it look nicer

    nstyle = ete3.NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 8
    nstyle["fgcolor"] = "darkred"

    nstyle2 = ete3.NodeStyle()
    nstyle2["shape"] = "sphere"
    nstyle2["size"] = 8
    nstyle2["fgcolor"] = "steelblue"

    nstyle3 = ete3.NodeStyle()
    nstyle3["size"] = 1
    nstyle3["fgcolor"] = "blue"

    nstyle4 = ete3.NodeStyle()
    nstyle4["shape"] = "sphere"
    nstyle4["size"] = 8
    nstyle4["fgcolor"] = "darkgreen"

    nstyle5 = ete3.NodeStyle()
    nstyle5["shape"] = "sphere"
    nstyle5["size"] = 8
    nstyle5["fgcolor"] = "darkorange"


    max_use = 0
    n_root = None
    for n in t.traverse():
        if n.dist > max_use:
            n_root = n
            max_use = n.dist
        if "BMH" in n.name:
            print(n.name, name_dict[n.name])
            n.name = name_dict[n.name]
            n.set_style(nstyle2)
        elif "serotype" in name_dict[n.name]:
            n.name = name_dict[n.name]
            if "serotype 4" in n.name:
                n.set_style(nstyle)
            elif "serotype 1/2" in n.name:
                n.set_style(nstyle4)
            elif "serotype 3" in n.name:
                n.set_style(nstyle5)
        else:
            n.set_style(nstyle3)
            n.name = ""

    t.set_outgroup(n_root)
    ts = ete3.TreeStyle()
    t.render(outfile,tree_style=ts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile_tree")
    parser.add_argument("--outfile")
    parser.add_argument("--infile_genomes", default="/media/Data/LabOnAChip/PhylogeneticTrees/Genomes/ListeriaGenomes.csv")
    parser.add_argument("--samples",action="append",default=[])
    parser.add_argument("--labels",action="append",default=[])
    args = parser.parse_args()

    # Create name dict
    assert len(args.samples) == len(args.labels)
    name_dict = defaultdict(str)
    for i in range(len(args.samples)):
        name_dict[args.samples[i]] = args.labels[i]

    # Update the name_dict
    name_dict, dict_info_oi = update_meta(args.infile_genomes, name_dict)
    
    t = load_tree(args.infile_tree)
    plot_tree(t, name_dict, dict_info_oi, args.outfile)
