
RUNS = ["BMH-2022-000226","BMH-2022-000227"]
#RUNS = []
INDIR = "/media/Data/LabOnAChip/230518_Aug3-2022"
OUTDIR = "/media/Data/LabOnAChip/230518_Aug3-2022/Results/"

# Mapping
mapping_db="Data/230518_straindoublecheck/combined_strain.fasta"
mapping_str="threegenomesstrain"
#mapping_db="Data/combined_reference.fasta"
#mapping_str="threegenomesRef"

#Quality trimming options
qtrim="rl" #Trim both sides
quality=12 #Trim bases below this quality score
quality_read=12 #Feature counts will ignore reads below this quality score
min_len=20 #Exclude reads below this length
adaptors='Data/adapters.fa' #Path to the adaptors that are trimmed by bbduk

species_trees=["s__Listeria_monocytogenes"]
species_trees_genomes={"s__Listeria_monocytogenes": "/media/Data/LabOnAChip/PhylogeneticTrees/Genomes/ListeriaGenomes/"}
#species_trees=["Lmonocytogenes_4b"]
#species_trees_genomes={"Lmonocytogenes_4b": "/media/Data/LabOnAChip/PhylogeneticTrees/Genomes/ListeriaGenomes/"}
infile_markers="/media/Data/LabOnAChip/PhylogeneticTrees/Genomes/ListeriaCoreGenes99.fasta"

downsample=True
