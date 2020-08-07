################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo2017_preimplantation/settings.R")
} else {
  stop("Computer not recognised")
}

# Folders with the differential analysis results
io$diff.acc <- paste0(io$basedir,"/acc/differential/feature_level")

# Output directory
io$outdir <- paste0(io$basedir,"/acc/boxplots")

####################
## Define options ##
####################

# Define stages
opts$stages <- c(
  "Zygote",
  "2-cell",
  "4-cell",
  "8-cell",
  "Morula",
  "ICM",
  "TE"
)

# Define genomic contexts
opts$acc.annos <- c(
  # "2cell_H3K27me3",
  # "4cell_3PN_H3K4me3",
  # "4cell_H3K27me3",
  # "8cell_3PN_H3K27me3",
  # "8cell_H3K27ac",
  # "CGI",
  # "H1_H3K27ac",
  # "H1_H3K27me3",
  # "H1_H3K4me1",
  # "H1_H3K4me3",
  # "H1_distal_H3K27ac",
  # "ICM_H3K27ac",
  # "ICM_H3K4me3",
  # "TE_H3K27me3",
  # "atac_peaks_2cell_2pn",
  # "atac_peaks_2cell_3pn",
  # "atac_peaks_4cell_3pn",
  # "atac_peaks_8cell_2pn",
  # "atac_peaks_8cell_3pn",
  # "atac_peaks_icm_2pn",
  # "distal_4cell_3PN_H3K4me3",
  # "distal_8cell_H3K27ac",
  # "distal_ICM_H3K27ac",
  # "distal_ICM_H3K4me3",
  # "genebody",
  # "prom_2000_2000",
  # "prom_2000_2000_cgi",
  # "prom_2000_2000_noncgi",
  # "prom_200_200",
  "prom_200_200_cgi",
  "prom_200_200_noncgi"
  # "proximal_4cell_3PN_H3K4me3",
  # "proximal_8cell_H3K27ac",
  # "proximal_ICM_H3K27ac",
  # "proximal_ICM_H3K4me3"
)

# Options for selecting differential hits
opts$min.fdr <- 0.10
opts$min.acc.diff <- 25


# Update metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages] %>%
  .[!is.na(id_acc) & pass_accQC==T]

