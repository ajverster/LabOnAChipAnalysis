library('htmlwidgets')
library('plotly')
library('tidyverse')
library('scales')
library('Biostrings')

args <- commandArgs(trailingOnly=TRUE)

infile_meta <- args[1]
df_meta <- read_csv(infile_meta)
type_tag <- args[2]

indir <- "/media/Data/LabOnAChip/"
window_size_use <- 1000
min_qual <- 20

translation_dict <- list()
translation_dict[["EcoliO157EDL933"]] <- "s__Escherichia_coli"
translation_dict[["Lmonocytogenes_4b"]] <- "s__Listeria_monocytogenes"
translation_dict[["Senterica_ATCC14028"]] <- "s__Salmonella_enterica"

ylimits <- c("Lmonocytogenes_4b"=490, "EcoliO157EDL933"=150, "Senterica_ATCC14028"=80)
do_line_plot <- function(df_mapped_nonoverlapping, name) {
  levels(df_mapped_nonoverlapping$infile)  <- c("Powerblade","Standard")

  for (sp in unique(df_mapped_nonoverlapping$species)) {
    if (sp%in%names(translation_dict)) {
      sp_write <- translation_dict[[sp]]
    } else {
      sp_write <- sp
    }
    sp_use <- str_replace_all(sp,"s__","") %>% str_replace_all("_"," ")
    p <- ggplot(filter(df_mapped_nonoverlapping, (species == sp)), aes(x = pos)) + geom_line(alpha=0.5, aes(y=covg, color=infile, group=infile)) + scale_color_manual(values = c("#4682b4", "#D74B4B"))  + theme_bw(16) + scale_x_continuous(labels=comma) + theme(panel.grid=element_blank(), axis.text = element_text(color="black"), legend.position="none") + xlab("") + ylab("Coverage") + ggtitle(sp_use) + ylim(0, ylimits[sp])
  
    # This was showing line breaks according to the location of the contig breaks
    # p <- add_contig_line_breaks(p, contig_lengths)

    pp <- ggplotly(p) %>% layout(legend=list(x=0.015,y=0.97,title=list(text=""),bgcolor="#FFFFFFFF",borderwidth=1, bordercolor="grey"))
    saveWidget(pp, file.path(outdir,sprintf("%s_%s_%s.html",sp_write,name,type_tag)), selfcontained = F, libdir = "lib")
    png(file.path(outdir,sprintf("%s_%s_%s.png",sp_write,name,type_tag)), width=800*1.75, height=800)
    print(p)
    dev.off()
  }
}

add_contig_line_breaks <- function(p, contig_lengths) {
      if (length(contig_lengths) > 0) { 
      s <- 0
      for (i in 1:length(contig_lengths)) {
        s <- s + contig_lengths[[i]]
        if (contig_lengths[[i]] < 10000)
          break
        p <- p + geom_vline(xintercept=s, linetype="dashed", size=0.1, alpha=0.5)
      }
    }
    return(p)
}


#####
# Do the SNP Plot
####
do_snp_plot <- function(df_mapped_nonoverlapping, df_meta, name) {

  indir <- dirname(dirname(infile_mapped))
  for (sp in unique(df_mapped_nonoverlapping$species)) {
    #if (sp != "EcoliO157EDL933")
    #  next

   df_plot_line <- filter(df_mapped_nonoverlapping, species==sp)

    if (sp%in%names(translation_dict)) {
      sp_write <- translation_dict[[sp]]
    } else {
      sp_write <- sp
    }
    sp_use <- str_replace_all(sp,"s__","") %>% str_replace_all("_"," ")

    ## Load up the SNP vcf files
    df_list <- list()
    for (subdir in list.dirs(indir, recursive=FALSE)) {
      if (str_detect(subdir, "Result"))
        next
      
      sample <- basename(subdir)
      infile <- file.path(subdir, sprintf("%s_threegenomes%s_filter.vcf.gz", sample, type_tag))
      df <- read_tsv(infile, comment="##") %>% dplyr::rename(CHROM=`#CHROM`) %>% filter(QUAL >= min_qual) %>% filter(CHROM==sp)
      print(paste0("Number of SNP lines ",nrow(df)))
      df_list[[subdir]] <- df
    }
    
    # Find unique SNPs
    df_snps_2_sub <- dplyr::select(df_list[[2]], CHROM, POS, REF, ALT) %>% dplyr::rename(REF2=REF,ALT2=ALT)
    df_snps_merge <- dplyr::select(df_list[[1]], CHROM, POS, REF, ALT) %>% full_join(df_snps_2_sub, on=c("CHROM","POS")) %>%
      filter(((nchar(REF) == 1) | (is.na(REF))) & ((nchar(REF2) == 1) | (is.na(REF2)))) %>%
      filter(((nchar(ALT) == 1) | (is.na(ALT))) & ((nchar(ALT2) == 1) | (is.na(ALT2))))
    
    df_snps_merge_sp <- filter(df_snps_merge, CHROM==sp)
    
    pos_both <- df_snps_merge_sp %>% filter(!((ALT != ALT2) | (is.na(ALT) & !is.na(ALT2)) | (!is.na(ALT) & is.na(ALT2)))) %>% pull(POS)
    pos_one <- df_snps_merge_sp %>% filter(!is.na(ALT) & is.na(ALT2)) %>% pull(POS)
    pos_two <- df_snps_merge_sp %>% filter(is.na(ALT) & !is.na(ALT2)) %>% pull(POS)
    
    #save(list = ls(all.names = TRUE), file = "image.RData", envir = environment())
    #qwe

    stopifnot(length(unique(c(pos_one, pos_two))) == length(pos_one) + length(pos_two))
    
    df_snps <- data.frame(pos=c(pos_one,pos_two), Samples=c(rep(basename(names(df_list)[1]), times=length(pos_one)), rep(basename(names(df_list)[2]), times=length(pos_two))))
    
    # Create the two dataframes for plotting
    df_plot_line <- filter(df_mapped_nonoverlapping, species==sp)
    df_plot_line$Samples <- str_replace(df_plot_line$infile, sprintf("_threegenomes%s_mapped.pileup",type_tag),"")
    
    #save(list = ls(all.names = TRUE), file = "image.RData", envir = environment())
    #qwe

    p <- NULL
    pl <- NULL

    print(sprintf("Writing to %s",file.path(outdir,sprintf("%s_%s_%s_withSNPs.html",sp_write,name,type_tag))))
    p <- ggplot() + geom_line(data=df_plot_line, aes(x = pos, y = covg, color=Samples), alpha=0.5)

    if (nrow(df_snps) > 0) {
      df_snps_full <- df_snps %>% left_join(df_meta)
      df_summary_all <- data.frame()
        # Count the unique SNPs in each bin
        for (i in seq(from=0,to=max(df_snps_full$pos),by=window_size_use)) {
          df_snps_sub <- filter(df_snps_full, (Name == name) & (pos > i) & (pos <= i + window_size_use))
          if (nrow(df_snps_sub) == 0)
            next
          df_summary <- df_snps_sub %>% group_by(Samples) %>% dplyr::summarize(n_snps=length(pos))
          df_summary$pos <- (i + window_size_use/2)
          df_summary_all <- bind_rows(df_summary_all, df_summary)
        }
    df_plot_point <- df_summary_all %>% left_join(df_plot_line)
    p <- p + geom_point(data=df_plot_point, aes(x=pos, y=covg, color=Samples, size=n_snps), pch=19, show.legend=FALSE)

    } else {
        print(sprintf("No unique SNPs found, only making a line plot for species %s name %s", sp, name))
    }
  
  p <- p + theme_bw(16) + scale_color_manual(values = c("#4682b4", "#D74B4B")) + scale_fill_manual(values = c("#4682b4", "#D74B4B")) + scale_x_continuous(labels=comma) + theme(panel.grid=element_blank(), axis.text = element_text(color="black"), legend.position=c(.9,.75), legend.title=element_blank()) + xlab("") + ylab("Coverage") + ggtitle(sp_use)
  
    pl <- ggplotly(p, tooltip=c("n_snps","pos")) %>% style(hoverinfo = "skip", traces = c(1,2))
  pl <- style(pl, name = "Powerblade", traces = 1) %>% style(pl, name = "Standard", traces = 2)
  pl <- pl %>% layout(legend=list(x=0.015,y=0.97,title=list(text=""),bgcolor="#FFFFFFFF",borderwidth=1, bordercolor="grey"))
    saveWidget(pl, file.path(outdir,sprintf("%s_%s_%s_withSNPs.html",sp_write,name,type_tag)), selfcontained = F, libdir = "lib")

  }
}

print(df_meta)
df_meta_sub <- unique(dplyr::select(df_meta, Directory, Name))

for (i in 1:nrow(df_meta_sub)) {
  dir <- df_meta_sub$Directory[i]
  name <- df_meta_sub$Name[i]

  infile_mapped <- file.path(indir, dir, sprintf("Results/mapped_gc_content_threegenomes%s.csv", type_tag))
  outdir <- dirname(infile_mapped)

  df_mapped <- read_csv(infile_mapped)
  df_mapped$infile <- factor(df_mapped$infile)

  df_mapped <- filter(df_mapped, window_size == window_size_use)
  gap_size <- df_mapped$pos[2] - df_mapped$pos[1]

  min_pos <- df_mapped$pos %>% min
  max_pos <- df_mapped$pos %>% max
  df_mapped_nonoverlapping <- df_mapped %>% filter(pos%in%seq(min_pos, to=max_pos, by = window_size_use))
  
#  do_snp_plot(df_mapped_nonoverlapping, df_meta, name)
  do_line_plot(df_mapped_nonoverlapping, name)
}


