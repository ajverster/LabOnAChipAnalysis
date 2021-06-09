
library('tidyverse')
library('scales')

####
# AQHist Plot
####

infile <- "571834_aqhist.txt"
df1 <- read_tsv(infile)
df1$type <- "standard"
colnames(df1)[1] <- str_replace(colnames(df1)[1], "#", "")

infile <- "571835_aqhist.txt"
df2 <- read_tsv(infile)
df2$type <- "LabOnAChip1"
colnames(df2)[1] <- str_replace(colnames(df2)[1], "#", "")

infile <- "571836_aqhist.txt"
df3 <- read_tsv(infile)
df3$type <- "LabOnAChip2"
colnames(df3)[1] <- str_replace(colnames(df3)[1], "#", "")

df_plot <- bind_rows(df1, df2, df3)

df_plot_group <- df_plot %>% filter(Quality > 0) %>% group_by(type,bin=cut(Quality, breaks = seq(from=0,to=40,by=5), labels = paste(seq(from=0,to=35,by=5), "-", seq(from=5,to=40,by=5)))) %>% summarize(sum=sum(fraction1))

ggplot(df_plot_group, aes(x = bin, y = sum)) + geom_bar(stat="identity", aes(fill = type), position="dodge", color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + xlab("Quality score") + ylab("Percentage of reads") + theme(axis.text = element_text(color="black"), axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c("#4682b4", "#5e99c5", "#D74B4B"))

####
# QHist
####

infile <- "571834_qhist.txt"
df1 <- read_tsv(infile)
df1$type <- "standard"
colnames(df1)[1] <- str_replace(colnames(df1)[1], "#", "")

infile <- "571835_qhist.txt"
df2 <- read_tsv(infile)
df2$type <- "LabOnAChip1"
colnames(df2)[1] <- str_replace(colnames(df2)[1], "#", "")


infile <- "571836_qhist.txt"
df3 <- read_tsv(infile)
df3$type <- "LabOnAChip2"
colnames(df3)[1] <- str_replace(colnames(df3)[1], "#", "")


df_plot <- bind_rows(df1, df2, df3)


df_plot_sub <- df_plot %>% select(BaseNum, type, Read1_linear, Read2_linear) %>% pivot_longer(c(Read1_linear, Read2_linear))
df_plot_sub$name <- factor(df_plot_sub$name)
levels(df_plot_sub$name) <- c("Paired End #1", "Paired End #2")


ggplot(df_plot_sub,aes(x = BaseNum, y = value))  + geom_line(aes(group=type, color=type)) + facet_grid(. ~ name) + theme_bw(16) + scale_color_manual(values=c("#4682b4", "#5e99c5", "#D74B4B")) +
  ylim(c(0,35)) + ylab("Quality score") + xlab("Base position") + theme(axis.text = element_text(color="black"), legend.title=element_blank())

####
# Metaphlan Plots
####

load_metaphlan <- function(infile, type="genus") {
  df <- read_tsv(infile, comment = "#", col_names=FALSE)
  
  if (type == "genus") {
    df <- df[str_detect(df$X1,"g__[A-Za-z]+$"),]
    df$X1 <- sapply(df$X1, function(x) str_split(x, fixed("|"))[[1]][6])
  } else {
    df <- df[str_detect(df$X1,"s__[A-Za-z]+"),]
    df$X1 <- sapply(df$X1, function(x) str_split(x, fixed("|"))[[1]][7])
    
  }
  return(df)
}


# Genus level
infile <- "571834_metphlan.txt"
df1 <- load_metaphlan(infile)
df1$type <- "Standard"

infile <- "571835_metphlan.txt"
df2 <- load_metaphlan(infile)
df2$type <- "LabOnAChip1"


infile <- "571836_metphlan.txt"
df3 <- load_metaphlan(infile)
df3$type <- "LabOnAChip2"

df_plot <- bind_rows(df1, df2, df3)
df_plot <- df_plot %>% select(X1, X3, type) %>% rename(genus=X1, abundance=X3)
df_plot <- bind_rows(df_plot, data.frame(genus = df1$X1, abundance = c(66.6, 33.3, 0.09, 0.01), type = "Theoretical"))

df_plot$abundance <- df_plot$abundance / 100

ggplot(df_plot, aes(x = type, y = abundance)) + geom_bar(stat = "identity", aes(fill = genus), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + 
      theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.title=element_blank()) + xlab("") + ylab("MetaPhlAn abundance")


# Species level
infile <- "571834_metphlan.txt"
df1 <- load_metaphlan(infile, "species")
df1$type <- "Standard"

infile <- "571835_metphlan.txt"
df2 <- load_metaphlan(infile, "species")
df2$type <- "LabOnAChip1"


infile <- "571836_metphlan.txt"
df3 <- load_metaphlan(infile, "species")
df3$type <- "LabOnAChip2"

df_plot <- bind_rows(df1, df2, df3)
df_plot <- df_plot %>% select(X1, X3, type) %>% rename(species=X1, abundance=X3)
df_plot <- bind_rows(df_plot, data.frame(species = df1$X1, abundance = c(66.6, 33.3, 0.09, 0.01), type = "Theoretical"))

df_plot$abundance <- df_plot$abundance / 100

ggplot(df_plot, aes(x = type, y = abundance)) + geom_bar(stat = "identity", aes(fill = species), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + 
  theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.title=element_blank()) + xlab("") + ylab("MetaPhlAn abundance")



####
# Counts the GC content
####

count_gc <- function(infile) {
  df <- read.csv(infile, sep = ",", header=TRUE)
  return(filter(df, X%in%c("C","G")) %>% pull(percentage) %>% sum)
}

indirs <- c("571834","571835","571836")
results <- data.frame()
for (ind in indirs) {
  infile <- sprintf("%s/%s_R1.GC.txt", ind, ind)
  infile.mapped <- sprintf("%s/%s_R1.GC.mapped.txt", ind, ind)
  results <- bind_rows(results, data.frame(ind=ind, gc_full=count_gc(infile), gc_mapped=count_gc(infile.mapped)))
}
results$type <- c("Standard","LabOnAChip1","LabOnAChip2")


results <- read.table("gc_data.csv", sep = ",", header = TRUE)
results$gc_full <- results$gc_full / 100
results$gc_mapped <- results$gc_mapped / 100
results$group <- c("Standard","LabOnAChip","LabOnAChip")

df_plot <- pivot_longer(select(results,-ind), cols = gc_full:gc_mapped)
df_plot$name <- factor(df_plot$name)
levels(df_plot$name) <- c("All reads", "MetaPhlAn mapped reads")

ggplot(df_plot, aes(x = type, y = value)) + geom_bar(stat="identity", aes(fill=type), color="black") + facet_grid(. ~ name) + 
  theme_bw(16) + ylab("GC%") + scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = c("#4682b4", "#5e99c5", "#D74B4B")) + xlab("") + theme(axis.text = element_text(color = "black"), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1)) +
  theme(legend.position="none") + geom_text(fontface="bold",aes(y = value+.02, label=100 * round(value, digits=2)))


####
# GC content of the genomes
####

df_plot <- data.frame(species = c("Listeria monocytogenes","Escherichia coli"), gc = c(38.19227139522177, 50.36487291548511))
df_plot$gc <- df_plot$gc / 100

ggplot(df_plot, aes(x=species, y = gc)) + geom_bar(stat="identity", aes(fill = species), color = "black") + theme_bw(16) + scale_y_continuous(labels = scales::percent) + ylab("GC %") + theme(axis.text = element_text(color="black"), legend.position = "none") + scale_fill_brewer(palette="Set2") + xlab("")



####
# Assembly statistics
####

df_assembly <- read.table("assembly_stats.csv", sep = ",", header=TRUE)
df_assembly$group <- c("Standard","Standard","LabOnAChip1","LabOnAChip1","LabOnAChip2","LabOnAChip2")

ggplot(filter(df_mapped, species == "listeria"), aes(x = pos, y = covg)) + geom_line(alpha=0.5, aes(color=infile, group=infile)) + scale_color_manual(values = c("#D74B4B","#4682b4", "#5e99c5")) + 
                                             theme_bw(16) + theme(legend.title = element_blank())

ggplot(filter(df_assembly, metric=="l50"), aes(x = group, y = value)) + geom_bar(stat="identity",color="black", aes(fill = group)) + theme_bw(16) + scale_fill_manual(values=c("#4682b4", "#5e99c5", "#D74B4B")) + ylab("L50") + theme(legend.position = "none", axis.text = element_text(color="black"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")


#####
# GC content vs coverage in bins
#####

df_mapped <- read_csv("mapped_gc_content.csv")
df_mapped$infile <- factor(df_mapped$infile)
levels(df_mapped$infile)  <- c("Standard","LabOnAChip1","LabOnAChip2")

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
