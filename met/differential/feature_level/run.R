###################################################################################
## Script to compute (in parallel) differential methylation at the feature level ##
###################################################################################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$script <- "/Users/ricard/Guo_2017/met/differential/feature_level/diffmet_feature_level.R"
  io$outdir <- "/Users/ricard/data/Guo_2017/met/differential/feature_level"
  ricard/Guo_2017/met/differential/feature_level
} else if(grepl("yoda",Sys.info()['nodename'])){
  io$script <- "/homes/ricard/Guo_2017/met/differential/feature_level/diffmet_feature_level.R"
  io$outdir <- "/hps/nobackup2/stegle/users/ricard/Guo_2017/met/differential/feature_level"
  io$tmpdir <- "/hps/nobackup2/stegle/users/ricard/Guo_2017/met/differential/feature_level/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
dir.create(io$outdir, showWarnings=F, recursive = T)


## Options ##
opts <- list()

opts$groups <- list(
  "ICM_vs_TE" = list(c("ICM"), c("TE")),
  "Morula_vs_ICM" = list(c("Morula"), c("ICM")),
  "Morula_vs_TE" = list(c("Morula"), c("TE"))
)

# Minimum number of cells per group
opts$min.cells <- 5

# Genomic contexts
opts$anno <- c(
  "CGI",
  # "E3.5_H3K27ac_distal",
  "ESC_CTCF",
  "ESC_DHS",
  "ESC_p300",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "LINE",
  "LTR",
  "prom_2000_2000",
  "prom_2000_2000_cgi",
  "prom_2000_2000_noncgi"
  # "prom_200_200",
  # "prom_200_200_cgi",
  # "prom_200_200_noncgi"
  # "window2000_step250"
)
# opts$anno <- c("prom_2000_2000")

for (i in names(opts$groups)) {
  groupA <- opts$groups[[i]][[1]]
  groupB <- opts$groups[[i]][[2]]
  for (j in opts$anno) {
    outfile <- sprintf("%s/%s_%s.txt.gz", io$outdir, i, j)
    lsf <- sprintf("bsub -M 4096 -n 1 -q standard -o %s/%s_%s.txt", io$tmpdir, i, j)
    # lsf <- ""
    cmd <- sprintf("%s Rscript %s --anno %s --group1 %s --group2 %s --min.cells %d --outfile %s", 
                   lsf, io$script, j, paste(groupA, collapse=" "), paste(groupB, collapse=" "), opts$min.cells, outfile)
    print(cmd)
    system(cmd)
  }
}
