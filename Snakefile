# Adrian Verster

# conda activate 16S
# snakemake -j 3 -p

import config
from pathlib import Path
import pandas as pd

RUNS = config.RUNS
INDIR = config.INDIR

rule all:
    input:
        expand("%s/{run}/qhist.txt" %(INDIR), run=RUNS),
        expand("%s/{run}/{run}_metaphlan.txt" %(INDIR), run=RUNS),
        expand("%s/{run}/MetaSpades/" %(INDIR), run=RUNS),
        "%s/mapped_gc_content.csv" %(INDIR),
        "%s/assembly_stats.csv" %(INDIR),
        expand("%s/{run}/{run}_R1.GC.txt" %(INDIR),run=RUNS),
        expand("%s/{run}/{run}_%s_filter.vcf" %(INDIR, config.mapping_str), run=RUNS),
        expand("%s/{run}/{run}_%s_consensus.fasta" %(INDIR, config.mapping_str), run=RUNS),
        

rule gzip:
    input:
        R1="%s/{run}/{run}_R1.fastq" %(INDIR),
        R2="%s/{run}/{run}_R2.fastq" %(INDIR)
    output:
        R1="%s/{run}/{run}_R1.fastq.gz" %(INDIR),
        R2="%s/{run}/{run}_R2.fastq.gz" %(INDIR)
    run:
        shell("gzip {input.R1}"),
        shell("gzip {input.R2}")


rule qc_filter:
    input:
        R1=ancient("%s/{run}/{run}_R1.fastq.gz" %(INDIR)),
        R2=ancient("%s/{run}/{run}_R2.fastq.gz" %(INDIR)),
    output:
        R1="%s/{run}/{run}_R1.qc.fastq.gz" %(INDIR),
        R2="%s/{run}/{run}_R2.qc.fastq.gz" %(INDIR),
        qhist="%s/{run}/qhist.txt" %(INDIR),
        lhist="%s/{run}/lhist.txt" %(INDIR),
        aqhist="%s/{run}/aqhist.txt" %(INDIR),
    log:
        filter_stats="%s/{run}/QualityFiltering_Stats_bbduk.txt" %(INDIR)
    params:
        qtrim=config.qtrim,
        quality=config.quality,
        min_len=config.min_len,
        adaptors=config.adaptors,
    shell:
        '''
        bbduk.sh in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} maq={params.quality} ref={params.adaptors} qtrim={params.qtrim} trimq={params.quality} minlength={params.min_len} tpe tbo qhist={output.qhist} lhist={output.lhist} aqhist={output.aqhist} overwrite=t forcetrimright=150 2> {log.filter_stats}
        '''

rule count_reads:
    input:
        R1_all=expand("%s/{run}/{run}_R1.qc.fastq.gz" %(INDIR), run=RUNS),
        R2_all=expand("%s/{run}/{run}_R2.qc.fastq.gz" %(INDIR), run=RUNS),
    output:
        "%s/read_counts.csv" %(INDIR)
    run:
        shell("python lib/count_library_reads.py {input.R1_all} {input.R2_all} {output}")


rule downsample:
    input:
        R1="%s/{run}/{run}_R1.qc.fastq.gz" %(INDIR),
        R2="%s/{run}/{run}_R2.qc.fastq.gz" %(INDIR),
        read_counts="%s/read_counts.csv" %(INDIR)
    output:
        R1="%s/{run}/{run}_R1.qc.downsample.fastq.gz" %(INDIR),
        R2="%s/{run}/{run}_R2.qc.downsample.fastq.gz" %(INDIR),
    run:
        df_counts = pd.read_csv(input.read_counts)
        min_counts = int(df_counts["n_reads"].min())
        print(min_counts)
        outfile_R1_nonzip = output.R1.replace(".gz","")
        outfile_R2_nonzip = output.R2.replace(".gz","")
        shell("seqtk sample -s100 {input.R1} %i > %s" %(min_counts, outfile_R1_nonzip))
        shell("seqtk sample -s100 {input.R2} %i > %s" %(min_counts, outfile_R2_nonzip))
        shell("gzip %s" %(outfile_R1_nonzip))
        shell("gzip %s" %(outfile_R2_nonzip))

def input_reads(wildcards):
    if config.mode == "downsample":
        return {'R1':str(INDIR) + "/{wildcards.run}/{wildcards.run}_R1.qc.downsample.fastq.gz".format(wildcards=wildcards), 'R2': str(INDIR) + "/{wildcards.run}/{wildcards.run}_R2.qc.downsample.fastq.gz".format(wildcards=wildcards)}
    else:
        return {'R1':str(INDIR) + "/{wildcards.run}/{wildcards.run}_R1.qc.fastq.gz".format(wildcards=wildcards), 'R2': str(INDIR) + "/{wildcards.run}/{wildcards.run}_R2.qc.fastq.gz".format(wildcards=wildcards)}


rule metaphlan:
    input:
        unpack(input_reads)
    output:
        ab="%s/{run}/{run}_metaphlan.txt" %(INDIR),
        bowtie="%s/{run}/{run}.metaphlanbowtie.gz" %(INDIR),
    run:
        print(input.R1)
        print(input.R2)
        shell("/home/adrian/anaconda3/envs/16S/bin/metaphlan {input.R1},{input.R2} --input_type fastq -o {output.ab} --nproc 8 --bowtie2out {output.bowtie}")


# I didn't end up using this. I thought it might be helpful to look at the abundances of the individual marker genes, but it ended up leading me astray.
#rule metaphlan_marker_ab:
#    input:
#        unpack(input_reads)
#        output:
#                ab="%s/{run}/{run}_metphlan_geneabundance.txt" %(INDIR),
#                bowtie="%s/{run}/{run}.metaphlanbowtie_alt.gz" %(INDIR),
#    shell:
#        "/home/adrian/anaconda3/envs/16S/bin/metaphlan {input.R1},{input.R2} --input_type fastq -o {output.ab} --nproc 8 --bowtie2out {output.bowtie} -t marker_ab_table" 


# Basic assembly. Might want to try this with the --meta option as well.
rule metaspades:
    input:
        unpack(input_reads)
                #R1="%s/{run}/{run}_R1.qc.fastq.gz" %(INDIR),
                #R2="%s/{run}/{run}_R2.qc.fastq.gz" %(INDIR),
    output:
        "%s/{run}/MetaSpades/" %(INDIR)
    shell:
        "metaspades.py -1 {input.R1} -2 {input.R2} -o {output}"


rule metaspades_stats:
    input:
        expand("%s/{run}/MetaSpades/" %(INDIR), run=RUNS)
    output:
        "%s/assembly_stats.csv" %(INDIR)
    run:
        infiles = [Path(x) / "scaffolds.fasta" for x in input]
        infiles = [str(x) for x in infiles]
        shell("python lib/quantify_assemblies.py %s %s" %(" ".join(infiles), INDIR))


####
# Should count SNPs but I didn't use it
####
rule midas:
    input:
        unpack(input_reads)
                #R1="%s/{run}/{run}_R1.qc.fastq.gz" %(INDIR),
                #R2="%s/{run}/{run}_R2.qc.fastq.gz" %(INDIR),
    output:
        "%s/{run}/MIDAS" %(INDIR)
    shell:
        "/home/adrian/MIDAS/scripts/run_midas.py species {output} -1 {input.R1} -2 {input.R2} -t 16 -d /media/Data/LabOnAChip/MIDAS/midas_db_v1.2"

####
# Maps to the file with both listeria and ecoli in one fasta file
####
rule map:
    input:
        unpack(input_reads)
                #R1="%s/{run}/{run}_R1.qc.fastq.gz" %(INDIR),
                #R2="%s/{run}/{run}_R2.qc.fastq.gz" %(INDIR),
    output:
        "%s/{run}/{run}_%s_mapped.sam" %(INDIR, config.mapping_str)
    params:
        db=config.mapping_db
    shell:
        "bbmap.sh in={input.R1} in2={input.R2} out={output} ref={params.db} ambig=toss"

####
# This currently uses k=1, so it just counts GC content
####


def input_reads_and_bowtie(wildcards):
        if config.mode == "downsample":
                return {'R1':str(INDIR) + "/{wildcards.run}/{wildcards.run}_R1.qc.downsample.fastq.gz".format(wildcards=wildcards), 'R2': str(INDIR) + "/{wildcards.run}/{wildcards.run}_R1.qc.downsample.fastq.gz".format(wildcards=wildcards), 'metaphlan_bowtie': str(INDIR) + "/{wildcards.run}/{wildcards.run}.metaphlanbowtie.gz"}
        else:
                return {'R1':str(INDIR) + "/{wildcards.run}/{wildcards.run}_R1.qc.fastq.gz".format(wildcards=wildcards), 'R2': str(INDIR) + "/{wildcards.run}/{wildcards.run}_R1.qc.fastq.gz".format(wildcards=wildcards), 'metaphlan_bowtie': str(INDIR) + "/{wildcards.run}/{wildcards.run}.metaphlanbowtie.gz"}


rule count_kmers:
    input:
        unpack(input_reads_and_bowtie)
    output:
        GC="%s/{run}/{run}_R1.GC.txt" %(INDIR),
        GC_mapped="%s/{run}/{run}_R1.GC.mapped.txt" %(INDIR)
    run:
        shell("python lib/count_kmers.py {input.R1} {output.GC}"),
        shell("python lib/count_kmers.py {input.R1} {output.GC_mapped} {input.metaphlan_bowtie}")


rule pileup:
    input:
        "%s/{run}/{run}_%s_mapped.sam" %(INDIR, config.mapping_str)
    output:
        bam="%s/{run}/{run}_%s_mapped.bam" %(INDIR, config.mapping_str), 
        bam_sorted="%s/{run}/{run}_%s_mapped_sorted.bam" %(INDIR, config.mapping_str),
        pileup="%s/{run}/{run}_%s_mapped.pileup" %(INDIR, config.mapping_str)
    params:
        db=config.mapping_db
    run:
        shell("samtools view -S -b {input} > {output.bam}")
        shell("samtools sort {output.bam} -o {output.bam_sorted}")
        shell("samtools mpileup {output.bam_sorted} -f {params.db}  > {output.pileup}")

rule pileup_all:
    input:
        expand("%s/{run}/{run}_%s_mapped.pileup" %(INDIR, config.mapping_str), run=RUNS)

rule free_bayes:
    input:
        bam_sorted="%s/{run}/{run}_%s_mapped_sorted.bam" %(INDIR, config.mapping_str),
    output:
        "%s/{run}/{run}_%s.vcf" %(INDIR, config.mapping_str),
    params:
        db=config.mapping_db
    shell:
        "freebayes -f {params.db} --ploidy 1 {input} > {output}"

rule filter_free_bayes:
    input:
        "%s/{run}/{run}_%s.vcf" %(INDIR, config.mapping_str),
    output:
        gz="%s/{run}/{run}_%s_filter.vcf.gz" %(INDIR, config.mapping_str),
    run:
        shell("cat {input} | vcffilter -f 'QUAL > 20' -f 'DP > 10' > %s" %(output.gz.replace(".gz","")))
        shell("bgzip %s" %(output.gz.replace(".gz","")))
        shell("bcftools index {output}")

rule consensus_sequence:
    input:
        "%s/{run}/{run}_%s_filter.vcf.gz" %(INDIR, config.mapping_str),
    output:
        "%s/{run}/{run}_%s_consensus.fasta" %(INDIR, config.mapping_str),
    params:
        db=config.mapping_db
    shell:
        "bcftools consensus {input} --fasta-ref {params.db} > {output}"

rule phylophlan_db:
    output:
        "/media/Data/LabOnAChip/PhylogeneticTrees/phylophlan_databases/{species}"
    shell:
        "phylophlan_setup_database -g {wildcards.species} -o {output} --verbose"

rule phylophlan:
    input:
        consensus=expand("%s/{run}/{run}_%s_consensus.fasta" %(INDIR, config.mapping_str), run=RUNS),
        db="/media/Data/LabOnAChip/PhylogeneticTrees/phylophlan_databases/{species}"
    output:
        "%s/phylophlan_{species}/" %(INDIR)
    run:
        genome_folder = config.species_trees_genomes[wildcards.species]
        db_folder = Path(input.db).parents[0]
        temp_dir = "/tmp/%s_%s" %(Path(INDIR).name, wildcards.species)
        # copy the new genome into db_folder
        cmd_run = "python lib/prep_genome_directory.py --species {wildcards.species} --indir_genomes %s --outdir_use %s" %(genome_folder, temp_dir)
        for infile_consensus in input.consensus:
            cmd_run += " --infile_consensus " + infile_consensus
        shell(cmd_run)
        shell("phylophlan -i %s -o {output} -d {wildcards.species} -t a -f ../PhylogeneticTrees/isolates_config.cfg --nproc 16 --subsample twentyfivepercent --diversity low --fast --databases_folder %s" %(temp_dir, db_folder))

rule phylophlan_all:
    input:
        expand("%s/phylophlan_{species}/" %(INDIR), species=config.species_trees)

rule consensus_prokka:
    input:
        consensus="%s/{run}/{run}_%s_consensus.fasta" %(INDIR, config.mapping_str)
    output:
        prokka=directory("%s/{run}_%s_consensus_{species}_prokka/" %(INDIR, config.mapping_str)),
        tr="%s/{run}_%s_consensus_{species}_prokka/{run}_{species}.ffn" %(INDIR, config.mapping_str),
        db="%s/{run}_%s_consensus_{species}_prokka/{run}_{species}_nuc_BLASTDB.nhr" %(INDIR, config.mapping_str),
        seq="%s/{run}/{run}_%s_consensus_{species}.fasta" %(INDIR, config.mapping_str),
    resources:
        cpus=8
    run:
        prefix = wildcards.run + "_" + wildcards.species
        shell("python ../PhylogeneticTrees/bioinformatics_parsing/src/bioinformatics_parsing/parse_sequences.py --infile_seqs {input} --genes_oi {wildcards.species} --outfile {output.seq}")
        shell("prokka --outdir {output.prokka} --cpus {resources.cpus} {output.seq} --centre X --compliant --quiet --prefix %s --force" %(prefix))
        shell("makeblastdb -in {output.tr} -dbtype nucl -out %s" %(output.db.replace(".nhr","")))


rule phylogenetic_trees:
    input:
        #consensus=expand("%s/{run}/{run}_%s_consensus.fasta" %(INDIR, config.mapping_str), run=RUNS),
        prokka=expand("%s/{run}_%s_consensus_{{species}}_prokka/" %(INDIR, config.mapping_str), run=RUNS)
    output:
        outdir="%s/phylogenetic_tree_{species}/" %(INDIR),
        aln="%s/phylogenetic_tree_alignment_{species}.fasta" %(INDIR),
    params:
        infile_marker_genes=config.infile_markers
    run:
        genome_folder = config.species_trees_genomes[wildcards.species]
        #temp_dir = "/tmp/%s_%s" %(Path(INDIR).name, wildcards.species)
        # copy the new genome into db_folder
        #cmd_run = "python lib/prep_genome_directory.py --species {wildcards.species} --indir_genomes %s --outdir_use %s" %(genome_folder, temp_dir)
        #for infile_consensus in input.consensus:
        #    cmd_run += " --infile_consensus " + infile_consensus
        #shell(cmd_run)
        #print(temp_dir)
        
        shell("python ../PhylogeneticTrees/bioinformatics_parsing/src/bioinformatics_parsing/parse_phylogeny.py --outfile_aln {output.aln} --outdir_tree {output.outdir} --infile_markers {params.infile_marker_genes} --indir %s --indir_secondary %s --method fasttree" %(genome_folder, INDIR))

rule phylogenetic_trees_all:
    input:
        expand("%s/phylogenetic_tree_{species}/" %(INDIR), species=config.species_trees)


rule mapped_gc_bins:
    input:
        expand("%s/{run}/{run}_%s_mapped.pileup" %(INDIR, config.mapping_str), run=RUNS)
    output:
        "%s/mapped_gc_content.csv" %(INDIR)
    shell:
        "python lib/readpileup.py {input} {output}"

rule coverage_quantification:
    input:
        expand("%s/{run}/{run}_%s_mapped.pileup" %(INDIR, config.mapping_str), run=RUNS)
    output:
        "%s/mapped_coverage.csv" %(INDIR)
    params:
        db=config.mapping_db
    shell:
        "python lib/quantify_pileup_coverage.py {params.db} {input} {output}"
