
import csv
import sys
import json
import pandas as pd
from Bio import SeqIO

csv.field_size_limit(sys.maxsize)

infile = "mpa_v30_CHOCOPhlAn_201901_marker_info.txt"

def load_ab_species(infile, infile_marker_info = "mpa_v30_CHOCOPhlAn_201901_marker_info.txt"):
	df_ab = pd.read_csv(infile, sep = "\t", comment="#", header=None) 
	df_ab.columns = ["marker","abundance"]

	data = [] 
	with open(infile_marker_info) as f: 
		for line in csv.reader(f, delimiter="\t"): 
			if line[0] in df_ab["marker"].values: 
				species = json.loads(line[1].replace("'",'"'))['clade'] 
				marker = line[0]
				data.append([marker, species])
	df_species = pd.DataFrame(data, columns = ["marker","species"])
	df_ab = df_ab.merge(df_species, how = "left")
	return df_ab
	

df_ab1 = load_ab_species("571834/571834_metphlan_geneabundance.txt")
df_ab1.rename(columns={"abundance":"571834"}, inplace=True)

df_ab2 = load_ab_species("571835/571835_metphlan_geneabundance.txt")
df_ab2.rename(columns={"abundance":"571835"}, inplace=True)

df_ab3 = load_ab_species("571836/571836_metphlan_geneabundance.txt")
df_ab3.rename(columns={"abundance":"571836"}, inplace=True)

df_full = df_ab1.merge(df_ab2, how = "outer").merge(df_ab3, how = "outer").fillna(0)

samples_oi = ["571834","571835","571836"]
samples_oi_nrom = ["571834_norm","571835_norm","571836_norm"]
for key in samples_oi:
	df_full[key + "_norm"] = df_full[key] / df_full[key].sum()

df_full.query('species.str.contains("s__Listeria_monocytogenes")', engine="python")


df_full.query('species.str.contains("s__Listeria_monocytogenes")', engine="python").sort_values("571834_norm")[samples_oi_nrom]


df_full.query('species.str.contains("s__Escherichia_coli")', engine="python").sort_values("571834_norm")[samples_oi_nrom]

