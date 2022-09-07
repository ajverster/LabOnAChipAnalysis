
RUNS = ["BMH-2022-000100","BMH-2022-000101"]
#RUNS = []
INDIR = "/media/Data/LabOnAChip/Apr26-2022"

# Mapping
#mapping_db="Data/combined.fasta"
#mapping_str="threegenomes"
mapping_db="Data/combined_reference.fasta"
mapping_str="threegenomesRef"

#Quality trimming options
qtrim="rl" #Trim both sides
quality=12 #Trim bases below this quality score
quality_read=12 #Feature counts will ignore reads below this quality score
min_len=20 #Exclude reads below this length
adaptors='Data/adapters.fa' #Path to the adaptors that are trimmed by bbduk

species_trees=["s__Listeria_monocytogenes"]
species_trees_genomes={"s__Listeria_monocytogenes": "/media/Data/LabOnAChip/PhylogeneticTrees/Genomes/ListeriaGenomes/"}
infile_markers="/media/Data/LabOnAChip/PhylogeneticTrees/Genomes/ListeriaCoreGenes99.fasta"

downsample=True
