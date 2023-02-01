python lib/combine_for_db.py -i GCA_009990925.1_PDT000129940.2_genomic.fna.gz -n Lmonocytogenes_4b -i GCF_000732965.1_ASM73296v1_genomic.fna.gz -n EcoliO157EDL933 -i GCF_016864495.1_ASM1686449v1_genomic.fna.gz -n Senterica_ATCC14028 -o combined.fasta

python lib/combine_for_db.py -i GCA_009990925.1_PDT000129940.2_genomic.fna.gz -n s__Listeria_monocytogenes -i GCF_000732965.1_ASM73296v1_genomic.fna.gz -n s__Escherichia_coli -i GCF_016864495.1_ASM1686449v1_genomic.fna.gz -n s__Salmonella_enterica -o combined.fasta



python lib/combine_for_db.py -i GCF_000196035.1_ASM19603v1_genomic.fna.gz -n s__Listeria_monocytogenes -i GCF_000005845.2_ASM584v2_genomic.fna.gz -n s__Escherichia_coli -i GCF_000006945.2_ASM694v2_genomic.fna.gz -n s__Salmonella_enterica -o combined_reference.fasta --mode first
bowtie2-build combined_reference.fasta combined_reference.fasta

python lib/combine_for_db.py -i GCA_003031895.1_ASM303189v1_genomic.fna -n Lmonocytogenes_4b -i GCF_000732965.1_ASM73296v1_genomic.fna -n EcoliO157EDL933 -i GCF_016864495.1_ASM1686449v1_genomic.fna -n Senterica_ATCC14028 -o combined_strain.fasta --mode first
bowtie2-build combined_strain.fasta combined_strain.fasta

