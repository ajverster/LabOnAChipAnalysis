library('htmlwidgets')
library('plotly')
library('tidyverse')
library('scales')
library('Biostrings')

args <- commandArgs(trailingOnly=TRUE)

df_mapped <- read_csv(args[1])
date <- args[2]
window_size_use <- 1000
type_tag <- args[3]

dir <- dirname(args[1])
outdir <- file.path(dir, "Results")


# Load up the genoome contigs
contig_lengths_all <- list()
if (type_tag != "reference") {
    genome <- readDNAStringSet("Data/GCA_009990925.1_PDT000129940.2_genomic.fna")
    contig_lengths <- lapply(genome, length)
    contig_lengths_all[["s__Listeria_monocytogenes"]] <- contig_lengths
    genome <- readDNAStringSet("Data/GCF_000732965.1_ASM73296v1_genomic.fna")
    contig_lengths <- lapply(genome, length)
    contig_lengths_all[["s__Escherichia_coli"]] <- contig_lengths
    genome <- readDNAStringSet("Data/GCF_016864495.1_ASM1686449v1_genomic.fna")
    contig_lengths <- lapply(genome, length)
    contig_lengths_all[["s__Salmonella_enterica"]] <- contig_lengths
} else {
    contig_lengths_all[["s__Listeria_monocytogenes"]] <- list()
    contig_lengths_all[["s__Escherichia_coli"]] <- list()
    contig_lengths_all[["s__Salmonella_enterica"]] <- list()

}

translation_dict <- list()
translation_dict[["EcoliO157EDL933"]] <- "s__Escherichia_coli"
translation_dict[["Lmonocytogenes_4b"]] <- "s__Listeria_monocytogenes"
translation_dict[["Senterica_ATCC14028"]] <- "s__Salmonella_enterica"

df_mapped$infile <- factor(df_mapped$infile)
levels(df_mapped$infile)  <- c("Powerblade","Standard")

df_mapped <- filter(df_mapped, window_size == window_size_use)
gap_size <- df_mapped$pos[2] - df_mapped$pos[1]
df_mapped_nonoverlapping <- df_mapped[seq(from=1,to=nrow(df_mapped), by = window_size_use / gap_size),] 

for (sp in unique(df_mapped_nonoverlapping$species)) {
  if (sp%in%names(translation_dict)) {
    sp_write <- translation_dict[[sp]]
  } else {
    sp_write <- sp
  }
  #contig_lengths <- contig_lengths_all[[sp_write]]
  sp_use <- str_replace_all(sp,"s__","") %>% str_replace_all("_"," ")
  p <- ggplot(filter(df_mapped_nonoverlapping, (species == sp)), aes(x = pos)) + geom_line(alpha=0.5, aes(y=covg, color=infile, group=infile)) + scale_color_manual(values = c("#D74B4B","#4682b4"))  + theme_bw(16) + scale_x_continuous(labels=comma) + theme(panel.grid=element_blank(), axis.text = element_text(color="black"), legend.position=c(.9,.75), legend.title=element_blank()) + xlab("") + ylab("Coverage") + ggtitle(sp_use)

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


