
library('tidyverse')
library('scales')
library('yaml')
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

outfile.metaphlan.genus <- file.path(indir, params$outfile_metaphlan_genus) #"MetaphlanGenus.png"
outfile.metaphlan.species <- file.path(indir, params$outfile_metaphlan_species) #"MetaphlanSpecies.png"
outfile.gc <- file.path(indir, params$outfile_gc) #"GC_Coarse.png"
outfile.assembly <- file.path(indir, params$outfile_assembly)

load_bbduk_files <- function(infile, label) {
	df <- read_tsv(infile)
	df$type <- label
	colnames(df)[1] <- str_replace(colnames(df)[1], "#", "")
	return(df)
}

df_plot <- data.frame()
for (i in 1:length(samples)) {
  df1 <- load_bbduk_files(file.path(indir, samples[i], "aqhist.txt"), samples.names[i])
  df_plot <- bind_rows(df_plot, df1)
}


df_plot_group <- df_plot %>% filter(Quality > 0) %>% group_by(type,bin=cut(Quality, breaks = seq(from=0,to=40,by=5), labels = paste(seq(from=0,to=35,by=5), "-", seq(from=5,to=40,by=5)))) %>% summarize(sum=sum(fraction1))

png(outfile.quality, width=w.use, height=h.use)
ggplot(df_plot_group, aes(x = bin, y = sum)) + geom_bar(stat="identity", aes(fill = type), position="dodge", color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + xlab("Quality score") + ylab("Percentage of reads") + theme(axis.text = element_text(color="black"), axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=colors.use)
dev.off()


####
# QHist
####

df_plot <- data.frame()
for (i in 1:length(samples)) {
  df1 <- load_bbduk_files(file.path(indir,samples[i],"qhist.txt"), samples.names[i])
  df_plot <- bind_rows(df_plot, df1)
}


df_plot_sub <- df_plot %>% select(BaseNum, type, Read1_linear, Read2_linear) %>% pivot_longer(c(Read1_linear, Read2_linear))
df_plot_sub$name <- factor(df_plot_sub$name)
levels(df_plot_sub$name) <- c("Paired End #1", "Paired End #2")


png(outfile.qhist, width=w.use*1.6, height=h.use)
ggplot(df_plot_sub,aes(x = BaseNum, y = value))  + geom_line(aes(group=type, color=type)) + facet_grid(. ~ name) + theme_bw(16) + scale_color_manual(values=colors.use) +
  ylim(c(0,40)) + ylab("Quality score") + xlab("Base position") + theme(axis.text = element_text(color="black"), legend.title=element_blank())
dev.off()

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


df_plot <- data.frame()
for (i in 1:length(samples)) {
  df1 <- load_metaphlan(file.path(indir, samples[i], sprintf("%s_metaphlan.txt",samples[i])), samples.names[i])
  df_plot <- bind_rows(df_plot, df1)
}


df_plot <- df_plot %>% select(X1, X3, type) %>% rename(genus=X1, abundance=X3)

i <- 1
df1 <- load_metaphlan(file.path(indir, samples[i], sprintf("%s_metaphlan.txt",samples[i])), samples.names[i])

ab.use <- params$ab_use
if (length(ab.use) == 3) {
	tax <- c("g__Listeria", "g__Escherichia", "g__Salmonella")
} else
	tax <- c("g__Listeria", "g__Escherichia", "g__Salmonella", "g__Staphylococcus")

df_plot <- bind_rows(df_plot, data.frame(genus = tax, abundance = ab.use, type = "Theoretical"))
df_plot$abundance <- df_plot$abundance / 100


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

png(outfile.metaphlan.species, width=w.use, height=h.use)
ggplot(df_plot, aes(x = type, y = abundance)) + geom_bar(stat = "identity", aes(fill = species), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + 
  theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.title=element_blank()) + xlab("") + ylab("MetaPhlAn abundance")
dev.off()

####
# Counts the GC content
####

count_gc <- function(infile) {
  df <- read.csv(infile, sep = ",", header=TRUE)
  return(filter(df, X%in%c("C","G")) %>% pull(percentage) %>% sum)
}

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

####
# GC content of the genomes
####

df_plot <- data.frame(species = c("Listeria monocytogenes","Escherichia coli"), gc = c(38.19227139522177, 50.36487291548511))
df_plot$gc <- df_plot$gc / 100

ggplot(df_plot, aes(x=species, y = gc)) + geom_bar(stat="identity", aes(fill = species), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + ylab("GC %") + theme(axis.text = element_text(color="black"), legend.position = "none") + scale_fill_brewer(palette="Set2") + xlab("")


####
# Assembly statistics
####

df_assembly <- read.table(file.path(indir,"assembly_stats.csv"), sep = ",", header=TRUE)
df_meta <- data.frame(infile=samples, group=samples.names, colors=colors.use)
df_assembly <- left_join(df_assembly, df_meta)
print(df_assembly)

png(outfile.assembly, width=w.use*1.2, height=h.use)
p1 <- ggplot(filter(df_assembly, metric=="l50"), aes(x = group, y = value)) + geom_bar(stat="identity",color="black", aes(fill = group)) + theme_bw(16) + ylab("L50") + theme(legend.position = "none", axis.text = element_text(color="black"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + scale_fill_manual(values=unique(df_assembly$colors))
p2 <- ggplot(filter(df_assembly, metric=="n50"), aes(x = group, y = value, fill=color)) + geom_bar(stat="identity",color="black", aes(fill = group)) + theme_bw(16) + ylab("N50") + theme(legend.position = "none", axis.text = element_text(color="black"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + scale_fill_manual(values=unique(df_assembly$colors))
grid.arrange(p1,p2, ncol=2)
dev.off()

quit(status=0)

#####
# GC content vs coverage in bins. Slides 7-10 in the slideshow
#####

df_mapped <- read_csv(file.path("mapped_gc_content.csv"))
df_mapped$infile <- factor(df_mapped$infile)
levels(df_mapped$infile)  <- c("Standard","LabOnAChip1","LabOnAChip2")

#ggplot(filter(df_mapped, species == "listeria"), aes(x = pos, y = covg)) + geom_line(alpha=0.5, aes(color=infile, group=infile)) + scale_color_manual(values = c("#D74B4B","#4682b4", "#5e99c5")) + 
#                                             theme_bw(16) + theme(legend.title = element_blank())

ggplot(filter(df_mapped, (species == "listeria") & (covg > 10)),aes(x = gc, y = covg)) + geom_point(aes(fill=infile), alpha=0.5, color="black",pch=21) + facet_wrap(infile ~ ., ncol=2) + theme_bw(16) + 
  scale_x_continuous(labels = scales::percent) + ylab("Bin Coverage") + xlab("Bin GC%") + ggtitle("1000 bp bins, Listeria monocytogenes") + scale_fill_manual(values = c("#D74B4B","#4682b4", "#5e99c5"))

# Correlation between GC and coverage
for (infile_sub in unique(df_mapped$infile)) {
  df_mapped_sub <- filter(df_mapped, (species == "listeria") & (covg > 10) & (infile == infile_sub))
  print(infile_sub)
  print(cor(df_mapped_sub$gc, df_mapped_sub$covg, method = "spearman"))
}

df_mapped_nonoverlapping <- df_mapped[seq(from=1,to=nrow(df_mapped), by = 20),]

ggplot(filter(df_mapped_nonoverlapping, (species == "listeria") & (covg > 10)),aes(x = gc, y = covg)) + geom_point(aes(fill=infile), alpha=0.5, color="black",pch=21) + facet_wrap(infile ~ ., ncol=3) + theme_bw(16) +
scale_x_continuous(labels = scales::percent) + ylab("Bin Coverage") + xlab("Bin GC%") + ggtitle("1000 bp bins, non-overlapping, Listeria monocytogenes") + scale_fill_manual(values = c("#D74B4B","#4682b4", "#5e99c5")) + theme(legend.title=element_blank(), axis.text=element_text(color="black"))


for (infile_sub in unique(df_mapped$infile)) {
  df_mapped_sub <- filter(df_mapped, (species == "Ecoli") & (covg > 10) & (infile == infile_sub))
  print(infile_sub)
  print(cor(df_mapped_sub$gc, df_mapped_sub$covg, method = "spearman"))
}

ggplot(filter(df_mapped, (species == "Ecoli") & (covg > 10)),aes(x = gc, y = covg)) + geom_point(aes(fill=infile), alpha=0.5, color="black",pch=21) + facet_wrap(infile ~ ., ncol=3) + theme_bw(16) + 
  scale_x_continuous(labels = scales::percent) + ylab("Bin Coverage") + xlab("Bin GC%") + ggtitle("1000 bp bins, Escherichia coli
") + scale_fill_manual(values = c("#D74B4B","#4682b4", "#5e99c5")) + theme(legend.title=element_blank())


ggplot(filter(df_mapped_nonoverlapping, (species == "Ecoli") & (covg > 10)),aes(x = gc, y = covg)) + geom_point(aes(fill=infile), alpha=0.5, color="black",pch=21) + facet_wrap(infile ~ ., ncol=3) + theme_bw(16) + 
  scale_x_continuous(labels = scales::percent) + ylab("Bin Coverage") + xlab("Bin GC%") + ggtitle("1000 bp bins, non-overlapping, Escherichia coli
")  + scale_fill_manual(values = c("#D74B4B","#4682b4", "#5e99c5"))  + theme(legend.title=element_blank(), axis.text=element_text(color="black"))


ggplot(filter(df_mapped, (species == "Ecoli") & (covg > 10)),aes(x = gc, y = covg)) + stat_bin_2d() + facet_wrap(infile ~ ., ncol=2) + theme_bw(16) + 
  scale_x_continuous(labels = scales::percent) + ylab("Bin Coverage") + xlab("Bin GC%") + ggtitle("1000 bp bins, Ecoli")


#### Show a full coverage plot, over the course of the genome
df_plot <- filter(df_mapped, (species == "listeria") & (infile != "GC"))
ggplot(df_plot, aes(x = pos)) + geom_line(alpha=0.5, aes(y=covg, color=infile, group=infile)) + scale_color_manual(values = c("#D74B4B","#4682b4", "#5e99c5"))  + theme_bw(16) + theme(legend.title=element_blank()) + scale_x_continuous(labels=comma) + theme(axis.text = element_text(color="black")) + xlab("") + ylab("Coverage") + ggtitle("Listeria monocytogenes")

df_plot <- filter(df_mapped, (species == "Ecoli") & (infile != "GC"))
ggplot(df_plot, aes(x = pos)) + geom_line(alpha=0.5, aes(y=covg, color=infile, group=infile)) + scale_color_manual(values = c("#D74B4B","#4682b4", "#5e99c5"))  + theme_bw(16) + 
  theme(legend.title=element_blank()) + scale_x_continuous(labels=comma) + theme(axis.text = element_text(color="black")) + xlab("") + ylab("Coverage") + ggtitle("Escherichia coli") + ylim(c(0,2000))



start <- 120000
width <- 10000
df_plot <- filter(df_mapped, species == "listeria", pos > start, pos < (start + width))
var(df_plot$gc)

trans <- 5000
ggplot(df_plot, aes(x = pos)) + geom_line(data=filter(df_plot, infile=="Standard"), lwd=1.5,aes(y=gc * trans, color="GC")) + geom_line(alpha=0.5, aes(y=covg, color=infile, group=infile)) + scale_y_continuous(sec.axis = sec_axis(~./trans, name="GC", labels = scales::percent)) + scale_color_manual(values = c("black","#4682b4", "#5e99c5","#D74B4B"))  + theme_bw(16) +
  theme(legend.title=element_blank()) + scale_x_continuous(labels=comma) + theme(axis.text = element_text(color="black")) + xlab("") + ylab("Coverage") + ggtitle("Listeria monocytogenes, arbitrary region")


# This is in one of the genomic islands
df_plot <- filter(df_mapped, species == "Ecoli", pos > 460000, pos < 470000)

trans <- 1500
ggplot(df_plot, aes(x = pos)) + geom_line(data=filter(df_plot, infile=="Standard"), lwd=1.5,aes(y=gc * trans, color="GC")) + geom_line(alpha=0.5, aes(y=covg, color=infile, group=infile)) + scale_y_continuous(name="coverage", sec.axis = sec_axis(~./trans, name="GC", labels = scales::percent)) + scale_color_manual(values = c("black","#4682b4", "#5e99c5","#D74B4B"))  + theme_bw(16) + 
  theme(legend.title=element_blank()) + scale_x_continuous(labels=comma) + theme(axis.text = element_text(color="black")) + xlab("") + ylab("Coverage") + ggtitle("Escherichia coli, genomic island")



df_plot <- filter(df_mapped, species == "Ecoli", pos > 390000, pos < 400000)
var(df_plot$gc)

trans <- 1500
ggplot(df_plot, aes(x = pos)) + geom_line(data=filter(df_plot, infile=="Standard"), lwd=1.5,aes(y=gc * trans, color="GC")) + geom_line(alpha=0.5, aes(y=covg, color=infile, group=infile)) + scale_y_continuous(name="coverage", sec.axis = sec_axis(~./trans, name="GC", labels = scales::percent)) + scale_color_manual(values = c("black","#4682b4", "#5e99c5","#D74B4B"))  + theme_bw(16)

