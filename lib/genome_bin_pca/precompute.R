
library('tidyverse')
library('Hmisc')


filter_bins <- function(df_gc_content, min_cov=10, min_pos=0, max_pos=2.75e6) {
  return(filter(df_gc_content, (covg > min_cov) & (pos >= min_pos) & (pos <= max_pos)))
}

nonoverlapping_bins <- function(df_gc_content, bin_size=1000) {
  df_gc_content_nonoverlapping <- data.frame()
  for (bin_size in unique(df_gc_content$window_size)) {
    df_gc_content_sub <- dplyr::filter(df_gc_content, window_size==bin_size)
    pos_all <- df_gc_content_sub$pos %>% unique %>% sort
    step_size <- pos_all[2] - pos_all[1]
    pos_sub_use <- pos_all[seq(from=1,to=length(pos_all), by = bin_size/step_size)]
    
    df_gc_content_nonoverlapping_bin <- filter(df_gc_content_sub, pos%in%pos_sub_use)
    df_gc_content_nonoverlapping <- bind_rows(df_gc_content_nonoverlapping, df_gc_content_nonoverlapping_bin)
  }
  return(df_gc_content_nonoverlapping)
}

remove_poor_samples <- function(df_gc_content, min_count=2000) {
  # Remove bad genomes where we see few regions with decent coverage
  infile_drop <- df_gc_content %>% group_by(infile) %>% dplyr::summarize(count=length(covg)) %>% filter(count < min_count) %>% pull(infile)
  return(filter(df_gc_content, !infile%in%infile_drop))
}

normalize_coverage <- function(df_gc_content) {
  return(df_gc_content %>% group_by(infile, window_size) %>% dplyr::summarize(pos=pos,gc=gc,covg=covg,covg_z=(covg - mean(covg)) / sd(covg)))
}

correlation_data <- function(df_gc_content, x=0.48, y=4) {
  # Correlation between GC and coverage
  df_text <- data.frame()
  for (infile_sub in unique(df_gc_content$infile)) {
    df_gc_content_sub <- filter(df_gc_content,  (infile == infile_sub))
    df_text <- bind_rows(df_text, data.frame(infile=infile_sub,x=x, y=y, text=sprintf("rho=%.2f",cor(df_gc_content_sub$gc, df_gc_content_sub$covg_z, method = "spearman"))))
    print(sprintf("%s: %.2f, mean covg %.2f", infile_sub, cor(df_gc_content_sub$gc, df_gc_content_sub$covg_z, method = "spearman"), mean(df_gc_content_sub$covg)))
  }
  return(df_text)
}


process_gc_bins <- function(infile, species_oi="s__Listeria_monocytogenes") {
  
  df_gc_content <- read_csv(infile)
  df_gc_content <- filter(df_gc_content, species == species_oi)
  df_gc_content_sub <- filter_bins(df_gc_content, min_cov=30)
  df_gc_content_sub <- nonoverlapping_bins(df_gc_content_sub)
  df_gc_content_sub <- remove_poor_samples(df_gc_content_sub)
  df_gc_content_sub <- normalize_coverage(df_gc_content_sub)
  df_gc_content_sub <- tibble(df_gc_content_sub)
  return(df_gc_content_sub)
}

pca_multi <- function(infiles, bin_size=1000, species_oi="s__Listeria_monocytogenes") {
  
  df_pca <- data.frame()
  for (inf in infiles) {
    df_gc_content_sub <- process_gc_bins(inf, species_oi)
    df_gc_content_sub <- filter(df_gc_content_sub, window_size == bin_size)
    df_pca_single <- df_gc_content_sub %>% select(-gc,-covg,-window_size) %>% pivot_wider(names_from=pos,values_from=covg_z)
    df_pca <- bind_rows(df_pca, df_pca_single)
  }
  pca <- df_pca %>% ungroup %>% select(-infile) %>% data.matrix %>% impute %>% prcomp

  df_plot <- pca$x %>% data.frame
  df_plot$infile <- df_pca %>% pull(infile)
  return(df_plot)
  
}



args <- commandArgs(trailingOnly=TRUE)
infile_metadata <- args[1]
genome <- args[2]

df_meta <- read_csv(infile_metadata)

indir <- "/media/Data/LabOnAChip/"
outdir <- "/media/Data/LabOnAChip/Pipeline/lib/genome_bin_pca/"

expt_dirs <- unique(df_meta$Directory)

infiles <- sapply(expt_dirs, function(x) file.path(indir, x, sprintf("Results/mapped_gc_content_threegenomes%s.csv", genome)))

if (genome == "strain") {
  strain_names <- c("Lmonocytogenes_4b","EcoliO157EDL933","Senterica_ATCC14028")
} else if (genome == "Ref") {
  strain_names <- c("s__Listeria_monocytogenes","s__Escherichia_coli","s__Salmonella_enterica")
}

df_list <- list()

for (sp_oi in strain_names) {
    for (bin_size in c(100,1000,10000)) {
      df_pca <- pca_multi(infiles, bin_size, sp_oi)
      df_pca$Samples <- sapply(df_pca$infile, function(x) str_split(x,"_")[[1]][1])
      df_list[[paste0(bin_size,sp_oi)]] <- df_pca %>% select(-infile) %>% left_join(df_meta, on="Samples")
    }
}


save.image(file.path(outdir,"pca_data.RData"))

