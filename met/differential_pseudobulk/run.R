###########################################################################
## Script to compute differential methylation using pseudobulked samples ##
###########################################################################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$script <- "/Users/ricard/Guo2017_preimplantation/met/differential_pseudobulk/diffmet_pseudobulk.R"
  io$outdir <- "/Users/ricard/data/Guo2017_preimplantation/met/differential/pseudobulk"
} else if (grepl("ricard",Sys.info()['nodename'])) {
  # io$script <- "/homes/ricard/Guo2017_preimplantation/met/differential/feature_level/diffmet_supervised.R"
  # io$outdir <- "/hps/nobackup/stegle/users/ricard/Guo2017_preimplantation/met/differential/feature_level/test"
  # io$tmpdir <- "/hps/nobackup/stegle/users/ricard/Guo2017_preimplantation/met/differential/feature_level/tmp"
} else {
  stop("Computer not recognised")
}
dir.create(io$outdir, showWarnings=F)


## Options ##
opts <- list()

opts$groups <- list(
  "ICM_vs_TE" = list(c("ICM"), c("TE")),
  "Morula_vs_ICM" = list(c("Morula"), c("ICM")),
  "Morula_vs_TE" = list(c("Morula"), c("TE"))
)

# Genomic contexts
opts$anno <- c(
  "CGI"
  # "E3.5_H3K27ac_distal",
  # "ESC_CTCF",
  # "ESC_DHS",
  # "ESC_p300",
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  # "H3K27ac_distal_E7.5_End_intersect12",
  # "H3K27ac_distal_E7.5_Mes_intersect12",
  # "LINE",
  # "LTR",
  # "prom_2000_2000",
  # "prom_2000_2000_cgi",
  # "prom_2000_2000_noncgi",
  # "prom_200_200",
  # "prom_200_200_cgi",
  # "prom_200_200_noncgi"
  # "window2000_step250"
)

for (group in names(opts$groups)) {
  stage1 <- opts$groups[[group]][[1]]
  stage2 <- opts$groups[[group]][[2]]
  for (anno in opts$anno) {
    outfile <- sprintf("%s/%s_%s.txt", io$outdir, group, anno)
    # lsf <- sprintf("bsub -M 2048 -n 1 -q standard -o %s/%s_%s.txt", io$tmpdir, group, anno)
    lsf <- ""
    cmd <- sprintf("%s Rscript %s --anno %s --stage1 %s --stage2 %s --outfile %s", 
                   lsf, io$script, anno, paste(stage1, collapse=" "), paste(stage2, collapse=" "), outfile)
    system(cmd)
  }
}
