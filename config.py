#RUNS = ["571834","571835","571836"]
#RUNS = ["BMH-2022-000098","PB-AO2"]
#INDIR = "/media/Data/LabOnAChip/Mar23-2022"
RUNS = ["BMH-2022-000100","BMH-2022-000101"]
INDIR = "/media/Data/LabOnAChip/Apr26-2022"

# Mapping
mapping_db="Data/combined.fasta"
mapping_str="threegenomes"

#Quality trimming options
qtrim="rl" #Trim both sides
quality=12 #Trim bases below this quality score
quality_read=12 #Feature counts will ignore reads below this quality score
min_len=20 #Exclude reads below this length
adaptors='Data/adapters.fa' #Path to the adaptors that are trimmed by bbduk

species_trees=["s__Listeria_monocytogenes"]
species_trees_genomes={"s__Listeria_monocytogenes": "/media/Data/LabOnAChip/PhylogeneticTrees/Genomes/ListeriaGenomes/"}

mode="downsample"
