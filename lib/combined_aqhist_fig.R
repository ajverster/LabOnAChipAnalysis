library('tidyverse')
library('scales')
library('yaml')
library('eulerr')
library('gridExtra')

params <- read_yaml("params_graphing.yaml")
#indir <- params$indir
#outdir <- params$outdir
#samples <- params$samples #c("BMH-2022-000100","BMH-2022-000101")
#samples.names <- params$samples_names #c("LabOnAChip", "Standard")
colors.use <- params$colors_use #c("#4682b4", "#D74B4B")
w.use <- params$w_use #800
h.use <- params$h_use #800

infile_metadata <- params$infile_metadata
indir_root <- params$indir_root

load_bbduk_files <- function(infile, label) {
    df <- read_tsv(infile)
    df$type <- label
    colnames(df)[1] <- str_replace(colnames(df)[1], "#", "")
    return(df)
}

quality_plots <- function(indir, samples, samples.names, outfile.plot, outfile.df) {
  df_plot <- data.frame()
  for (i in 1:length(samples)) {
    df1 <- load_bbduk_files(file.path(indir, samples[i], "aqhist.txt"), samples.names[i])
    df1_cleaned <- df1 %>% filter(Quality > 0) %>% group_by(type,bin=cut(Quality, breaks = seq(from=0,to=40,by=5), labels = paste(seq(from=0,to=35,by=5), "-", seq(from=5,to=40,by=5)))) %>% summarize(sum=sum(fraction1))
    df_plot <- bind_rows(df_plot, df1_cleaned)
  }

  # Write to file
  df_plot_group <- df_plot %>% group_by(type,bin) %>% summarize(mean = mean(sum), sd = sd(sum), .groups = "drop") %>% write_csv(outfile.df)

  p <- ggplot(df_plot_group, aes(x = bin, y = mean, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.8, position = position_dodge(0.6)) +
    theme_classic(16) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .03))) +
    xlab("Quality score") +
    ylab("Percentage of reads") +
    theme(axis.text = element_text(color = "black", size=12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position=c(0.2,0.8), legend.title= element_blank(), axis.ticks.x = element_blank(), axis.title=element_text(size=14,face="bold"), panel.border = element_rect(color = "#000000", fill = NA, size = 1)) +
    scale_fill_manual(values = colors.use)
  
  png(outfile.plot, width = w.use, height = h.use)
  print(p)
  dev.off()
}

df_metadata <- read_csv(infile_metadata)
df_metadata <- df_metadata %>% arrange(Type)


dir <- unique(df_metadata$Directory)[1]
indir <- file.path(indir_root, dir)
outdir <- file.path(indir,"Results/")
prefix <- paste0(basename(indir),"-")

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

    quality_plots(indir, group_df$Samples, group_df$Type, outfile.quality, outfile.aqhist.df)
}
