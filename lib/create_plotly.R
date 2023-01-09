library('htmlwidgets')
library('plotly')
library('tidyverse')
library('scales')
library('Biostrings')

args <- commandArgs(trailingOnly=TRUE)

infile_mapped <- args[1]
df_mapped <- read_csv(infile_mapped)
date <- args[2]
type_tag <- args[3]
outdir <- dirname(infile_mapped)
window_size_use <- 1000
min_qual <- 20

translation_dict <- list()
translation_dict[["EcoliO157EDL933"]] <- "s__Escherichia_coli"
translation_dict[["Lmonocytogenes_4b"]] <- "s__Listeria_monocytogenes"
translation_dict[["Senterica_ATCC14028"]] <- "s__Salmonella_enterica"

df_mapped$infile <- factor(df_mapped$infile)

df_mapped <- filter(df_mapped, window_size == window_size_use)
gap_size <- df_mapped$pos[2] - df_mapped$pos[1]

min_pos <- df_mapped$pos %>% min
max_pos <- df_mapped$pos %>% max
df_mapped_nonoverlapping <- df_mapped %>% filter(pos%in%seq(min_pos, to=max_pos, by = window_size_use))

do_line_plot <- function(df_mapped_nonoverlapping) {
  levels(df_mapped_nonoverlapping$infile)  <- c("Powerblade","Standard")

  for (sp in unique(df_mapped_nonoverlapping$species)) {
    if (sp%in%names(translation_dict)) {
      sp_write <- translation_dict[[sp]]
    } else {
      sp_write <- sp
    }
    sp_use <- str_replace_all(sp,"s__","") %>% str_replace_all("_"," ")
    p <- ggplot(filter(df_mapped_nonoverlapping, (species == sp)), aes(x = pos)) + geom_line(alpha=0.5, aes(y=covg, color=infile, group=infile)) + scale_color_manual(values = c("#4682b4", "#D74B4B"))  + theme_bw(16) + scale_x_continuous(labels=comma) + theme(panel.grid=element_blank(), axis.text = element_text(color="black"), legend.position=c(.9,.75), legend.title=element_blank()) + xlab("") + ylab("Coverage") + ggtitle(sp_use)

      if (FALSE) {
      if (length(contig_lengths) > 0) { 
      s <- 0
      for (i in 1:length(contig_lengths)) {
        s <- s + contig_lengths[[i]]
        if (contig_lengths[[i]] < 10000)
          break
        p <- p + geom_vline(xintercept=s, linetype="dashed", size=0.1, alpha=0.5)
      }
    }
    }

  
    pp <- ggplotly(p) %>% layout(legend=list(x=0.015,y=0.97,title=list(text=""),bgcolor="#FFFFFFFF",borderwidth=1, bordercolor="grey"))
    saveWidget(pp, file.path(outdir,sprintf("%s_%s_%s.html",sp_write,date,type_tag)), selfcontained = F, libdir = "lib")
  }
}


#####
# Do the SNP Plot
####
do_snp_plot <- function(df_mapped_nonoverlapping) {
  df_meta <- read_csv("metadata.txt")

  indir <- dirname(dirname(infile_mapped))
  for (sp in unique(df_mapped_nonoverlapping$species)) {

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
    stopifnot(length(unique(c(pos_one, pos_two))) == length(pos_one) + length(pos_two))
    
    df_snps <- data.frame(pos=c(pos_one,pos_two), sample=c(rep(basename(names(df_list)[1]), times=length(pos_one)), rep(basename(names(df_list)[2]), times=length(pos_two))))
    #df_snps_full <- bind_rows(df_snps_full, df_snps) %>% left_join(df_meta)
    df_snps_full <- df_snps %>% left_join(df_meta)
    print(filter(df_snps_full, (pos > 2787000) & (pos < 2788000)))

  # Count the unique SNPs in each bin
  df_summary_all <- data.frame()
  for (i in seq(from=0,to=max(df_snps_full$pos),by=window_size_use)) {
    df_snps_sub <- filter(df_snps_full, (Date == date) & (pos > i) & (pos <= i + window_size_use))
    if (nrow(df_snps_sub) == 0)
      next
    df_summary <- df_snps_sub %>% group_by(sample) %>% dplyr::summarize(n_snps=length(pos))
    df_summary$pos <- (i + window_size_use/2)
    df_summary_all <- bind_rows(df_summary_all, df_summary)
  }
  
  # Create the two dataframes for plotting
  df_plot_line <- filter(df_mapped_nonoverlapping, species==sp)
  df_plot_line$sample <- str_replace(df_plot_line$infile, sprintf("_threegenomes%s_mapped.pileup",type_tag),"")
  df_plot_point <- df_summary_all %>% left_join(df_plot_line)

  p <- NULL
  pl <- NULL
  p <- ggplot() + geom_line(data=df_plot_line, aes(x = pos, y = covg, color=sample), alpha=0.5) + geom_point(data=df_plot_point, aes(x=pos, y=covg, color=sample, size=n_snps), pch=19, show.legend=FALSE) + 
          theme_bw(16) + scale_color_manual(values = c("#4682b4", "#D74B4B")) + scale_fill_manual(values = c("#4682b4", "#D74B4B")) + scale_x_continuous(labels=comma) + theme(panel.grid=element_blank(), axis.text = element_text(color="black"), legend.position=c(.9,.75), legend.title=element_blank()) + xlab("") + ylab("Coverage") + ggtitle(sp_use)
  pl <- ggplotly(p, tooltip=c("n_snps","pos")) %>% style(hoverinfo = "skip", traces = c(1,2))

  pl <- style(pl, name = "Powerblade", traces = 1) %>% style(pl, name = "Standard", traces = 2)
  pl <- pl %>% layout(legend=list(x=0.015,y=0.97,title=list(text=""),bgcolor="#FFFFFFFF",borderwidth=1, bordercolor="grey"))
    saveWidget(pl, file.path(outdir,sprintf("%s_%s_%s_withSNPs.html",sp_write,date,type_tag)), selfcontained = F, libdir = "lib")

  }
}

do_snp_plot(df_mapped_nonoverlapping)
#do_line_plot(df_mapped_nonoverlapping)
