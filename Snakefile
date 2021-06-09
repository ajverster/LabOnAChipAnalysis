# Adrian Verster

import config

RUNS = ["571834","571835","571836"]

rule all:
	input:
		expand("{run}/qhist.txt", run=RUNS),
		expand("{run}/{run}_metphlan.txt", run=RUNS),
                expand("{run}/MetaSpades/", run=RUNS),
                expand("{run}/{run}_metphlan_geneabundance.txt",run=RUNS)

rule gzip:
	input:
		R1="{run}/{run}_R1.fastq",
		R2="{run}/{run}_R2.fastq"
	output:
		R1="{run}/{run}_R1.fastq.gz",
		R2="{run}/{run}_R2.fastq.gz"
	run:
		shell("gzip {input.R1}"),
		shell("gzip {input.R2}")


rule qc_filter:
	input:
		R1=ancient("{run}/{run}_R1.fastq.gz"),
		R2=ancient("{run}/{run}_R2.fastq.gz"),
	output:
		R1="{run}/{run}_R1.qc.fastq.gz",
		R2="{run}/{run}_R2.qc.fastq.gz",
		qhist="{run}/qhist.txt",
		lhist="{run}/lhist.txt",
		aqhist="{run}/aqhist.txt",
	log:
		filter_stats="{run}/QualityFiltering_Stats_bbduk.txt"
	params:
		qtrim=config.qtrim,
		quality=config.quality,
		min_len=config.min_len,
		adaptors=config.adaptors,
	shell:
		'''
		bbduk.sh in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} maq={params.quality} ref={params.adaptors} qtrim={params.qtrim} trimq={params.quality} minlength={params.min_len} tpe tbo qhist={output.qhist} lhist={output.lhist} aqhist={output.aqhist} overwrite=t 2> {log.filter_stats}
		'''

rule metaphlan:
	input:
                R1="{run}/{run}_R1.qc.fastq.gz",
                R2="{run}/{run}_R2.qc.fastq.gz",
	output:
		ab="{run}/{run}_metphlan.txt",
		bowtie="{run}/{run}.metaphlanbowtie.gz",
	shell:
		"/home/adrian/anaconda3/envs/16S/bin/metaphlan {input.R1},{input.R2} --input_type fastq -o {output.ab} --nproc 8 --bowtie2out {output.bowtie}" 

rule metaphlan_marker_ab:
        input:
                R1="{run}/{run}_R1.qc.fastq.gz",
                R2="{run}/{run}_R2.qc.fastq.gz",
        output:
                ab="{run}/{run}_metphlan_geneabundance.txt",
                bowtie="{run}/{run}.metaphlanbowtie.gz",
	shell:
		"/home/adrian/anaconda3/envs/16S/bin/metaphlan {input.R1},{input.R2} --input_type fastq -o {output.ab} --nproc 8 --bowtie2out {output.bowtie} -t marker_ab_table" 

rule metaspades:
	input:
                R1="{run}/{run}_R1.qc.fastq.gz",
                R2="{run}/{run}_R2.qc.fastq.gz",
	output:
		"{run}/MetaSpades/"
	shell:
		"metaspades.py -1 {input.R1} -2 {input.R2} -o {output}"


rule midas:
	input:
                R1="{run}/{run}_R1.qc.fastq.gz",
                R2="{run}/{run}_R2.qc.fastq.gz",
	output:
		"{run}/MIDAS"
	shell:
		"/home/adrian/MIDAS/scripts/run_midas.py species {output} -1 {input.R1} -2 {input.R2} -t 16 -d /media/Data/LabOnAChip/MIDAS/midas_db_v1.2"


rule map:
	input:
                R1="{run}/{run}_R1.qc.fastq.gz",
                R2="{run}/{run}_R2.qc.fastq.gz",
	output:
		"{run}/{run}_R1_ecoli_listeria_mapped.sam"
	shell:
		"bbmap.sh in={input.R1} out={output} ref=Data/combined_ecoli_listeria.fasta ambig=toss"


rule count_kmers:
	input:
                R1="{run}/{run}_R1.qc.fastq.gz",
                R2="{run}/{run}_R2.qc.fastq.gz",
		metaphlan_bowtie="{run}/{run}.metaphlanbowtie.gz",
	output:
		"{run}/{run}_R1.GC.mapped.txt"
	shell:
		"python count_kmers.py {input.R1} {output} {input.metaphlan_bowtie}"


rule pileup:
	input:
		"{run}/{run}_R1_ecoli_listeria_mapped.sam"
	output:
                bam="{run}/{run}_R1_ecoli_listeria_mapped.bam",
                bam_sorted="{run}/{run}_R1_ecoli_listeria_mapped_sorted.bam",
                pileup="{run}/{run}_R1_ecoli_listeria_mapped.pileup"
	run:
		shell("samtools view -S -b {input} > {output.bam}")
		shell("samtools sort {output.bam} -o {output.bam_sorted}")
		shell("samtools mpileup {output.bam_sorted} -f Data/combined_ecoli_listeria.fasta  > {output.pileup}")

rule mapped_gc_bins:
	input:
		expand("{run}/{run}_R1_ecoli_listeria_mapped.pileup", run=RUNS)
	output:
		"mapped_gc_content.csv"
	shell:
		"python lib/readpileup.py {input} {output}"
