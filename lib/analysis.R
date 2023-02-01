
library('tidyverse')
library('scales')
library('yaml')
library('eulerr')
library('gridExtra')

####
# Parameters
####

params <- read_yaml("params_graphing.yaml")
#indir <- params$indir
#outdir <- params$outdir
#samples <- params$samples #c("BMH-2022-000100","BMH-2022-000101")
#samples.names <- params$samples_names #c("LabOnAChip", "Standard")
colors.use <- params$colors_use #c("#4682b4", "#D74B4B")
w.use <- params$w_use #800
h.use <- params$h_use #800
min_ab_metaphlan <- params$min_ab_metaphlan

infile_metadata <- params$infile_metadata

# Convert ab_use from grams to molar
ab.use <- params$ab_use / params$genome_size
ab.use <- ab.use / sum(ab.use) * 100

indir_root <- params$indir_root

load_bbduk_files <- function(infile, label) {
	df <- read_tsv(infile)
	df$type <- label
	colnames(df)[1] <- str_replace(colnames(df)[1], "#", "")
	return(df)
}

####
# AQHist
####


quality_plots <- function(indir, samples, samples.names, outfile.plot, outfile.df) {
  df_plot <- data.frame()
  for (i in 1:length(samples)) {
    df1 <- load_bbduk_files(file.path(indir, samples[i], "aqhist.txt"), samples.names[i])
    df_plot <- bind_rows(df_plot, df1)
  }

  df_plot_group <- df_plot %>% filter(Quality > 0) %>% group_by(type,bin=cut(Quality, breaks = seq(from=0,to=40,by=5), labels = paste(seq(from=0,to=35,by=5), "-", seq(from=5,to=40,by=5)))) %>% summarize(sum=sum(fraction1))

  # Write to file
  df_plot_group %>% add_column(run=basename(indir)) %>% write_csv(outfile.df)

  p <- ggplot(df_plot_group, aes(x = bin, y = sum)) + geom_bar(stat="identity", aes(fill = type), position="dodge", color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + xlab("Quality score") + ylab("Percentage of reads") + theme(axis.text = element_text(color="black"), axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=colors.use)
  png(outfile.plot, width=w.use, height=h.use)
  print(p)
  dev.off()
}

####
# QHist
####

qhist_plots <- function(indir, samples, samples.names, outfile.plot, outfile.df) {
  df_plot <- data.frame()
  for (i in 1:length(samples)) {
    df1 <- load_bbduk_files(file.path(indir,samples[i],"qhist.txt"), samples.names[i])
    df_plot <- bind_rows(df_plot, df1)
  }


  df_plot_sub <- df_plot %>% select(BaseNum, type, Read1_linear, Read2_linear) %>% pivot_longer(c(Read1_linear, Read2_linear))
  df_plot_sub$name <- factor(df_plot_sub$name)
  levels(df_plot_sub$name) <- c("Paired End #1", "Paired End #2")

  # Write to file
  df_plot_sub %>% add_column(run=basename(indir)) %>% write_csv(outfile.df)

  p <- ggplot(df_plot_sub,aes(x = BaseNum, y = value))  + geom_line(aes(group=type, color=type)) + facet_grid(. ~ name) + theme_bw(16) + scale_color_manual(values=colors.use) +
    ylim(c(0,40)) + ylab("Quality score") + xlab("Base position") + theme(axis.text = element_text(color="black"), legend.title=element_blank())
  png(outfile.plot, width=w.use*1.6, height=h.use)
  print(p)
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

metaphlan_plots <- function(indir, samples, samples.names, outfile.plot.species, outfile.plot.genus) {
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

  p <- ggplot(df_plot, aes(x = type, y = abundance)) + geom_bar(stat = "identity", aes(fill = genus), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + 
        theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.title=element_blank()) + xlab("") + ylab("MetaPhlAn abundance")
  png(outfile.plot.genus, width=w.use, height=h.use)
  print(p)
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

  p <- ggplot(df_plot, aes(x = type, y = abundance)) + geom_bar(stat = "identity", aes(fill = species), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + 
    theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.title=element_blank()) + xlab("") + ylab("MetaPhlAn abundance")
  png(outfile.plot.species, width=w.use, height=h.use)
  print(p)
  dev.off()
}

####
# Counts the GC content
####

count_gc <- function(infile) {
  df <- read.csv(infile, sep = ",", header=TRUE)
  return(filter(df, X%in%c("C","G")) %>% pull(percentage) %>% sum)
}

gc_content_plots <- function(indir, samples, samples.names, outfile.plot) {
  results <- data.frame()
  for (ind in samples) {
    infile <- file.path(indir, ind, sprintf("%s_R1.GC.txt", ind))
    infile.mapped <- file.path(indir, ind, sprintf("%s_R1.GC.mapped.txt",ind))
    results <- bind_rows(results, data.frame(ind=ind, gc_full=count_gc(infile), gc_mapped=count_gc(infile.mapped)))
  }
  results$gc_full <- results$gc_full / 100
  results$gc_mapped <- results$gc_mapped / 100
  results$group <- samples.names 

  df_plot <- pivot_longer(select(results,-ind), cols = gc_full:gc_mapped)
  df_plot$name <- factor(df_plot$name)
  levels(df_plot$name) <- c("All reads", "MetaPhlAn mapped reads")


  p <- ggplot(df_plot, aes(x = group, y = value)) + geom_bar(stat="identity", aes(fill=group), color="black") + facet_grid(. ~ name) + 
    theme_bw(16) + ylab("GC%") + scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = colors.use) + xlab("") + theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1)) +
    theme(legend.position="none") + geom_text(fontface="bold",aes(y = value+.02, label=100 * round(value, digits=2)))
  png(outfile.plot, width=w.use, height=h.use)
  print(p)
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

assembly_plots <- function(indir, samples, samples.names, outfile.plot) {
  df_assembly <- read.table(file.path(outdir,"assembly_stats.csv"), sep = ",", header=TRUE)
  df_meta <- data.frame(infile=samples, group=samples.names, colors=colors.use)
  df_assembly <- left_join(df_assembly, df_meta)

  p1 <- ggplot(filter(df_assembly, metric=="l50"), aes(x = group, y = value)) + geom_bar(stat="identity",color="black", aes(fill = group)) + theme_bw(16) + ylab("L50") + theme(legend.position = "none", axis.text = element_text(color="black"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + scale_fill_manual(values=unique(df_assembly$colors))
  p2 <- ggplot(filter(df_assembly, metric=="n50"), aes(x = group, y = value, fill=color)) + geom_bar(stat="identity",color="black", aes(fill = group)) + theme_bw(16) + ylab("N50") + theme(legend.position = "none", axis.text = element_text(color="black"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + scale_fill_manual(values=unique(df_assembly$colors))
  png(outfile.plot, width=w.use*1.2, height=h.use)
  grid.arrange(p1,p2, ncol=2)
  dev.off()
}

####
# SNPs
####

snp_plots <- function(indir, samples, samples.names, outfile.plot) {
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

      #save(list = ls(all.names = TRUE), file = "debug.RData", envir = environment())

      VennDiag <- euler(c("A"=n_A, "B"=n_B, "A&B"=n_AB))
      png(str_replace(outfile.plot,"_sp",paste0("_",sp)), width=w.use*1.15, height=h.use)
      par(mar = c(20, 20, 20, 20))
      print(plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5,
          fill=c("#4682b4", "#D74B4B"), quantities=TRUE, labels = c(samples.names[1], samples.names[2])))
      dev.off()
  }

}

####
# Genome coveage - GC relationship
####

covg_gc_plots <- function(df_meta, outfile.plot, min_covg=10, window_size=1000) {
    
  df_full <- data.frame()
  for (expt in unique(df_meta$Directory)) {
    df <- read_csv(file.path(indir_root, expt,"Results/mapped_gc_content_threegenomesstrain.csv"))
    df_full <- bind_rows(df_full, df)

  }

  df_full$Samples <- sapply(df_full$infile,function(x) str_split(x,"_")[[1]][1])
  levels_use <- sort(unique(df_meta$Name))

  for (sp in unique(pull(df_full, species))) {
      df_plot <- df_full %>% filter((window_size == window_size) & (species == sp))
      df_plot <- left_join(df_plot, df_meta)
      df_plot$Name <- factor(df_plot$Name, ordered=TRUE, levels=levels_use)

      df_mapped_nonoverlapping <- df_plot[seq(from=1,to=nrow(df_plot), by = window_size/50),]

      df_text <- filter(df_mapped_nonoverlapping, (covg > min_covg)) %>% group_by(Samples) %>% dplyr::summarize(cor=cor(gc, covg, method="spearman")) %>% left_join(df_meta)
      df_text$gc <- (df_mapped_nonoverlapping %>% pull(gc) %>% max) - 0.3
      df_text$covg <- df_mapped_nonoverlapping %>% pull(covg) %>% max
      df_text$label <- sprintf("r = %.2f", df_text$cor)

      df_text$Name <- factor(df_text$Name, levels=levels_use)
      df_mapped_nonoverlapping$Name <- factor(df_mapped_nonoverlapping$Name, levels=levels_use)

      p <- ggplot(filter(df_mapped_nonoverlapping, (covg > min_covg)),aes(x = gc, y = covg)) + geom_point(aes(fill=Type), alpha=0.5, color="black",pch=21) + geom_text(data=df_text, size=8, aes(label=label)) +
        theme_bw(24) + facet_grid(Type ~ Name) +
        scale_x_continuous(labels = scales::percent) + ylab("Bin Coverage") + xlab("Bin GC%") + scale_fill_manual(values = c("#4682b4","#D74B4B")) + theme(legend.position="none")

      outfile.use <- str_replace(outfile.plot,"_sp",paste0("_",sp))
      print(sprintf("Writing to %s", outfile.use))

      png(outfile.use, width=w.use*1.75, height=h.use)
      print(p)
      dev.off()

  }
}

####
# Make the plots
####

df_metadata <- read_csv(infile_metadata)
df_metadata <- df_metadata %>% arrange(Type)


dir <- unique(df_metadata$Directory)[1]
indir <- file.path(indir_root, dir)
outdir <- file.path(indir,"Results/")
prefix <- paste0(basename(indir),"-")

outfile_coverage_correlation <- file.path(outdir, paste0(prefix,params$outfile_coverage_correlation))
covg_gc_plots(df_metadata, outfile_coverage_correlation)


# Iterate through each Directory and create plots 
for (dir in unique(df_metadata$Directory)) {
    group_df <- filter(df_metadata, Directory == dir)
    indir <- file.path(indir_root, dir)
    outdir <- file.path(indir,"Results/")
    prefix <- paste0(basename(indir),"-")

    outfile.quality <- file.path(outdir, paste0(prefix,params$outfile_quality)) #"QualityBins.png"
    outfile.qhist <- file.path(outdir, paste0(prefix,params$outfile_qhist)) #"QualityHist.png"
    outfile.aqhist.df <- file.path(outdir, paste0(prefix,params$outfile_aqhist_df))
    outfile.qhist.df <- file.path(outdir, paste0(prefix,params$outfile_qhist_df))

    outfile.metaphlan.genus <- file.path(outdir, paste0(prefix,params$outfile_metaphlan_genus)) #"MetaphlanGenus.png"
    outfile.metaphlan.species <- file.path(outdir, paste0(prefix,params$outfile_metaphlan_species)) #"MetaphlanSpecies.png"
    outfile.gc <- file.path(outdir, paste0(prefix,params$outfile_gc)) #"GC_Coarse.png"
    outfile.assembly <- file.path(outdir, paste0(prefix,params$outfile_assembly))
    outfile.snps <- file.path(outdir, paste0(prefix,params$outfile_snps))
    
    #quality_plots(indir, group_df$Samples, group_df$Type, outfile.quality, outfile.aqhist.df)
    #qhist_plots(indir, group_df$Samples, group_df$Type, outfile.qhist, outfile.qhist.df)
    #metaphlan_plots(indir, group_df$Samples, group_df$Type, outfile.metaphlan.species, outfile.metaphlan.genus)
    #gc_content_plots(indir, group_df$Samples, group_df$Type, outfile.gc)
    #assembly_plots(indir, group_df$Samples, group_df$Type, outfile.assembly)
    snp_plots(indir, group_df$Samples, group_df$Type, outfile.snps)


}
