
library('tidyverse')
library('scales')
library('yaml')
library('eulerr')
library('gridExtra')

####
# AQHist Plot
####

params <- read_yaml("params_graphing.yaml")
indir <- params$indir
samples <- params$samples #c("BMH-2022-000100","BMH-2022-000101")
samples.names <- params$samples_names #c("LabOnAChip", "Standard")
colors.use <- params$colors_use #c("#4682b4", "#D74B4B")
w.use <- params$w_use #800
h.use <- params$h_use #800
outfile.quality <- file.path(indir, params$outfile_quality) #"QualityBins.png"
outfile.qhist <- file.path(indir, params$outfile_qhist) #"QualityHist.png"
min_ab_metaphlan <- params$min_ab_metaphlan


prefix <- paste0(basename(indir),"-")
outfile.metaphlan.genus <- file.path(indir, paste0(prefix,params$outfile_metaphlan_genus)) #"MetaphlanGenus.png"
outfile.metaphlan.species <- file.path(indir, paste0(prefix,params$outfile_metaphlan_species)) #"MetaphlanSpecies.png"
outfile.gc <- file.path(indir, paste0(prefix,params$outfile_gc)) #"GC_Coarse.png"
outfile.assembly <- file.path(indir, paste0(prefix,params$outfile_assembly))
outfile.snps <- file.path(indir, paste0(prefix,params$outfile_snps))

outfile.aqhist.df <- file.path(indir, paste0(prefix,params$outfile_aqhist_df))
outfile.qhist.df <- file.path(indir, paste0(prefix,params$outfile_qhist_df))

# Convert ab_use from grams to molar
ab.use <- params$ab_use / params$genome_size
ab.use <- ab.use / sum(ab.use) * 100

indir_root <- params$indir_root
outfile_coverage_correlation <- params$outfile_coverage_correlation

load_bbduk_files <- function(infile, label) {
	df <- read_tsv(infile)
	df$type <- label
	colnames(df)[1] <- str_replace(colnames(df)[1], "#", "")
	return(df)
}


quality_plots <- function() {
  df_plot <- data.frame()
  for (i in 1:length(samples)) {
    df1 <- load_bbduk_files(file.path(indir, samples[i], "aqhist.txt"), samples.names[i])
    df_plot <- bind_rows(df_plot, df1)
  }


  df_plot_group <- df_plot %>% filter(Quality > 0) %>% group_by(type,bin=cut(Quality, breaks = seq(from=0,to=40,by=5), labels = paste(seq(from=0,to=35,by=5), "-", seq(from=5,to=40,by=5)))) %>% summarize(sum=sum(fraction1))

  # Write to file
  df_plot_group %>% add_column(run=basename(indir)) %>% write_csv(outfile.aqhist.df)

  png(outfile.quality, width=w.use, height=h.use)
  ggplot(df_plot_group, aes(x = bin, y = sum)) + geom_bar(stat="identity", aes(fill = type), position="dodge", color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + xlab("Quality score") + ylab("Percentage of reads") + theme(axis.text = element_text(color="black"), axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=colors.use)
  dev.off()
}

####
# QHist
####

qhist_plots <- function() {
  df_plot <- data.frame()
  for (i in 1:length(samples)) {
    df1 <- load_bbduk_files(file.path(indir,samples[i],"qhist.txt"), samples.names[i])
    df_plot <- bind_rows(df_plot, df1)
  }


  df_plot_sub <- df_plot %>% select(BaseNum, type, Read1_linear, Read2_linear) %>% pivot_longer(c(Read1_linear, Read2_linear))
  df_plot_sub$name <- factor(df_plot_sub$name)
  levels(df_plot_sub$name) <- c("Paired End #1", "Paired End #2")

  # Write to file
  df_plot_sub %>% add_column(run=basename(indir)) %>% write_csv(outfile.qhist.df)

  png(outfile.qhist, width=w.use*1.6, height=h.use)
  ggplot(df_plot_sub,aes(x = BaseNum, y = value))  + geom_line(aes(group=type, color=type)) + facet_grid(. ~ name) + theme_bw(16) + scale_color_manual(values=colors.use) +
    ylim(c(0,40)) + ylab("Quality score") + xlab("Base position") + theme(axis.text = element_text(color="black"), legend.title=element_blank())
  dev.off()
}
####
# Metaphlan Plots
####

load_metaphlan <- function(infile, metadata, level="genus") {
  df <- read_tsv(infile, comment = "#", col_names=FALSE)
  
  if (level == "genus") {
    df <- df[str_detect(df$X1,"g__[A-Za-z]+$"),]
    df$X1 <- sapply(df$X1, function(x) str_split(x, fixed("|"))[[1]][6])
  } else {
    df <- df[str_detect(df$X1,"s__[A-Za-z]+"),]
    df$X1 <- sapply(df$X1, function(x) str_split(x, fixed("|"))[[1]][7])
    
  }
  df$type <- metadata
  return(df)
}

metaphlan_plots <- function() {
  df_plot <- data.frame()
  for (i in 1:length(samples)) {
    df1 <- load_metaphlan(file.path(indir, samples[i], sprintf("%s_metaphlan.txt",samples[i])), samples.names[i])
    df_plot <- bind_rows(df_plot, df1)
  }


  df_plot <- df_plot %>% select(X1, X3, type) %>% rename(genus=X1, abundance=X3)

  i <- 1
  df1 <- load_metaphlan(file.path(indir, samples[i], sprintf("%s_metaphlan.txt",samples[i])), samples.names[i])

  if (length(ab.use) == 3) {
    tax <- c("g__Listeria", "g__Escherichia", "g__Salmonella")
  } else
    tax <- c("g__Listeria", "g__Escherichia", "g__Salmonella", "g__Staphylococcus")

  df_plot <- bind_rows(df_plot, data.frame(genus = tax, abundance = ab.use, type = "Theoretical"))
  df_plot$abundance <- df_plot$abundance / 100
  df_plot <- dplyr::filter(df_plot, abundance > min_ab_metaphlan)

  png(outfile.metaphlan.genus, width=w.use, height=h.use)
  ggplot(df_plot, aes(x = type, y = abundance)) + geom_bar(stat = "identity", aes(fill = genus), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + 
        theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.title=element_blank()) + xlab("") + ylab("MetaPhlAn abundance")
  dev.off()

  # Species level

  df_plot <- data.frame()
  for (i in 1:length(samples)) {
    df1 <- load_metaphlan(file.path(indir, samples[i], sprintf("%s_metaphlan.txt",samples[i])), samples.names[i], "species")
    df_plot <- bind_rows(df_plot, df1)
  }

  df_plot <- df_plot %>% select(X1, X3, type) %>% rename(species=X1, abundance=X3)

  if (length(ab.use) == 3) {
    tax <- c("s__Listeria_monocytogenes", "s__Escherichia_coli", "s__Salmonella_enterica")
  } else
    tax <- c("s__Listeria_monocytogenes", "s__Escherichia_coli", "s__Salmonella_enterica", "s__Staphylococcus_epidermidis")


  df_plot <- bind_rows(df_plot, data.frame(species = tax, abundance = ab.use, type = "Theoretical"))
  df_plot$abundance <- df_plot$abundance / 100
  df_plot <- dplyr::filter(df_plot, abundance > min_ab_metaphlan)

  png(outfile.metaphlan.species, width=w.use, height=h.use)
  ggplot(df_plot, aes(x = type, y = abundance)) + geom_bar(stat = "identity", aes(fill = species), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + 
    theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.title=element_blank()) + xlab("") + ylab("MetaPhlAn abundance")
  dev.off()
}

####
# Counts the GC content
####

count_gc <- function(infile) {
  df <- read.csv(infile, sep = ",", header=TRUE)
  return(filter(df, X%in%c("C","G")) %>% pull(percentage) %>% sum)
}

gc_content_plots <- function() {
  results <- data.frame()
  for (ind in samples) {
    infile <- file.path(indir, ind, sprintf("%s_R1.GC.txt", ind))
    infile.mapped <- file.path(indir, ind, sprintf("%s_R1.GC.mapped.txt",ind))
    results <- bind_rows(results, data.frame(ind=ind, gc_full=count_gc(infile), gc_mapped=count_gc(infile.mapped)))
  }
  results$gc_full <- results$gc_full / 100
  results$gc_mapped <- results$gc_mapped / 100
  results$group <- samples.names #c("Standard","LabOnAChip1","LabOnAChip2")

  df_plot <- pivot_longer(select(results,-ind), cols = gc_full:gc_mapped)
  df_plot$name <- factor(df_plot$name)
  levels(df_plot$name) <- c("All reads", "MetaPhlAn mapped reads")


  png(outfile.gc, width=w.use, height=h.use)
  ggplot(df_plot, aes(x = group, y = value)) + geom_bar(stat="identity", aes(fill=group), color="black") + facet_grid(. ~ name) + 
    theme_bw(16) + ylab("GC%") + scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = colors.use) + xlab("") + theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1)) +
    theme(legend.position="none") + geom_text(fontface="bold",aes(y = value+.02, label=100 * round(value, digits=2)))
  dev.off()
}

####
# GC content of the genomes
####

gc_content_genome_plots <- function() {
  df_plot <- data.frame(species = c("Listeria monocytogenes","Escherichia coli"), gc = c(38.19227139522177, 50.36487291548511))
  df_plot$gc <- df_plot$gc / 100

  ggplot(df_plot, aes(x=species, y = gc)) + geom_bar(stat="identity", aes(fill = species), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + ylab("GC %") + theme(axis.text = element_text(color="black"), legend.position = "none") + scale_fill_brewer(palette="Set2") + xlab("")
}

####
# Assembly statistics
####

assembly_plots <- function() {
  df_assembly <- read.table(file.path(indir,"assembly_stats.csv"), sep = ",", header=TRUE)
  df_meta <- data.frame(infile=samples, group=samples.names, colors=colors.use)
  df_assembly <- left_join(df_assembly, df_meta)

  png(outfile.assembly, width=w.use*1.2, height=h.use)
  p1 <- ggplot(filter(df_assembly, metric=="l50"), aes(x = group, y = value)) + geom_bar(stat="identity",color="black", aes(fill = group)) + theme_bw(16) + ylab("L50") + theme(legend.position = "none", axis.text = element_text(color="black"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + scale_fill_manual(values=unique(df_assembly$colors))
  p2 <- ggplot(filter(df_assembly, metric=="n50"), aes(x = group, y = value, fill=color)) + geom_bar(stat="identity",color="black", aes(fill = group)) + theme_bw(16) + ylab("N50") + theme(legend.position = "none", axis.text = element_text(color="black"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + scale_fill_manual(values=unique(df_assembly$colors))
  grid.arrange(p1,p2, ncol=2)
  dev.off()
}

####
# SNPs
####

snp_plots <- function() {
  infile_1 <- file.path(indir, samples[1], sprintf("%s_threegenomesstrain_filter.vcf.gz", samples[1]))
  infile_2 <- file.path(indir, samples[2], sprintf("%s_threegenomesstrain_filter.vcf.gz", samples[2]))
  df_snps_1 <- read_tsv(infile_1, comment="##") %>% dplyr::rename(CHROM=`#CHROM`) %>% filter(QUAL >= 20)
  df_snps_2 <- read_tsv(infile_2, comment="##") %>% dplyr::rename(CHROM=`#CHROM`) %>% filter(QUAL >= 20)

  df_snps_2_sub <- dplyr::select(df_snps_2, CHROM, POS, REF, ALT) %>% dplyr::rename(REF2=REF,ALT2=ALT)
  df_snps_merge <- dplyr::select(df_snps_1, CHROM, POS, REF, ALT) %>% full_join(df_snps_2_sub, on=c("CHROM","POS")) %>% 
    filter(((nchar(REF) == 1) | (is.na(REF))) & ((nchar(REF2) == 1) | (is.na(REF2)))) %>%
    filter(((nchar(ALT) == 1) | (is.na(ALT))) & ((nchar(ALT2) == 1) | (is.na(ALT2))))

  #Filter for species
  for (sp in unique(pull(df_snps_merge, CHROM))) {
      df_snps_merge_sp <- filter(df_snps_merge, CHROM==sp) 
      n_same <- df_snps_merge_sp %>% filter(!((ALT != ALT2) | (is.na(ALT) & !is.na(ALT2)) | (!is.na(ALT) & is.na(ALT2)))) %>% nrow
      n_one <- df_snps_merge_sp %>% filter(!is.na(ALT) & is.na(ALT2)) %>% nrow
      n_two <- df_snps_merge_sp %>% filter(is.na(ALT) & !is.na(ALT2)) %>% nrow
      n_diff <- df_snps_merge_sp %>% filter(ALT != ALT2) %>% nrow

      n_A <- n_one + n_diff
      n_B <- n_two + n_diff
      n_AB <- n_same


      VennDiag <- euler(c("A"=n_A, "B"=n_B, "A&B"=n_AB))
      png(str_replace(outfile.snps,"_sp",paste0("_",sp)), width=w.use*1.15, height=h.use)
      par(mar = c(20, 20, 20, 20))
      print(plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5,
          fill=c("#4682b4", "#D74B4B"), quantities=TRUE, labels = c(samples.names[1], samples.names[2])))
      dev.off()
  }

}

####
# Genome coveage - GC relationship
####

covg_gc_plots <- function() {
  expt_dirs <- c("Aug3-2022/","Jul26-2022/","Apr26-2022/")
    
  df_full <- data.frame()
  for (expt in expt_dirs) {
    df <- read_csv(file.path(indir_root, paste0(expt,"Results/mapped_gc_content_threegenomesstrain.csv")))
    df_full <- bind_rows(df_full, df)
    
  }
  df_full$Name <- sapply(df_full$infile,function(x) str_split(x,"_")[[1]][1])

  df_meta <- read_csv(file.path(indir_root, "Pipeline/metadata.txt"), col_names=FALSE)
  colnames(df_meta) <- c("Dataset","Name","Type")
    

  for (sp in unique(pull(df_full, species))) {
      df_plot <- df_full %>% filter((window_size == 1000) & (species == sp))
      df_plot <- left_join(df_plot, df_meta)
      df_plot$Dataset <- factor(df_plot$Dataset, ordered=TRUE, levels=c("Apr26","Jul26","Aug3"))

      df_mapped_nonoverlapping <- df_plot[seq(from=1,to=nrow(df_plot), by = 20),]

      df_text <- filter(df_mapped_nonoverlapping, (covg > 10)) %>% group_by(Name) %>% dplyr::summarize(cor=cor(gc, covg, method="spearman")) %>% left_join(df_meta)
df_text$gc <- (df_mapped_nonoverlapping %>% pull(gc) %>% max) - 0.3
      df_text$covg <- df_mapped_nonoverlapping %>% pull(covg) %>% max
      df_text$label <- sprintf("r = %.2f", df_text$cor)

      df_text$Dataset <- factor(df_text$Dataset, levels=c("Apr26","Jul26","Aug3"))
      df_mapped_nonoverlapping$Dataset <- factor(df_mapped_nonoverlapping$Dataset, levels=c("Apr26","Jul26","Aug3"))

      p <- ggplot(filter(df_mapped_nonoverlapping, (covg > 10)),aes(x = gc, y = covg)) + geom_point(aes(fill=Type), alpha=0.5, color="black",pch=21) + geom_text(data=df_text, size=8, aes(label=label)) +
        theme_bw(24) + facet_grid(Type ~ Dataset) +
        scale_x_continuous(labels = scales::percent) + ylab("Bin Coverage") + xlab("Bin GC%") + scale_fill_manual(values = c("#4682b4","#D74B4B")) + theme(legend.position="none")

      png(str_replace(outfile_coverage_correlation,"_sp",paste0("_",sp)), width=w.use*1.75, height=h.use)
      print(p)
      dev.off()
      

  }
}

covg_gc_plots()