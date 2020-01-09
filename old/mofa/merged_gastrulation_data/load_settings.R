library(data.table)
library(purrr)
library(ggplot2)
library(scater)

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir1 <- "/Users/ricard/data/gastrulation"
  io$basedir2 <- "/Users/ricard/data/Guo_2017"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/Users/ricard/data/Guo_2017/mofa"
} else {
  # io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  # io$basedir2 <- "/hps/nobackup/stegle/users/ricard/Guo_2017"
  # io$gene_metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  # io$outdir <- "/homes/ricard/gastrulation/metaccrna/mofa/E7.5/out"
}
io$sample.metadata1 <- paste0(io$basedir1,"/sample_metadata_scNMT.txt")
io$met.dir1 <- paste0(io$basedir1,"/met/parsed")
io$acc.dir1 <- paste0(io$basedir1,"/acc/parsed")
io$rna.file1 <- paste0(io$basedir1,"/rna/parsed/SingleCellExperiment.rds")

io$sample.metadata2 <- paste0(io$basedir2,"/sample_metadata.txt")
io$met.dir2 <- paste0(io$basedir2,"/met/parsed")
io$acc.dir2 <- paste0(io$basedir2,"/acc/parsed")

io$annos_dir  <- paste0(io$basedir1, "/features/filt")

## Define options ##
opts <- list()

# Define which annotations to look at
opts$met.annos <- c(
  "prom_2000_2000"
  # "H3K27ac_distal_E7.5_Mes_intersect12_500",
  # "H3K27ac_distal_E7.5_Ect_intersect12_500",
  # "H3K27ac_distal_E7.5_End_intersect12_500"
)

opts$acc.annos <- c(
  "prom_200_200"
  # "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  # "H3K27ac_distal_E7.5_End_intersect12"
)


# Define which stage and lineages to look at 
opts$stage_lineage1 <- c(
  # "E4.5_EPI","E4.5_PE",
  # "E5.5_EPI","E5.5_PE"
  # "E6.5_EPI", "E6.5_VE","E6.5_PS",
  "E7.5_Ectoderm","E7.5_Mesoderm","E7.5_Endoderm"
)

opts$stage_lineage2 <- c(
  # "ESC",
  # "TE",
  "ICM",
  "Morula",
  "16-cell",
  "8-cell",
  "4-cell",
  "2-cell",
  "Zygote"
)

# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 5      # minimum number of cells per feature (and per stage_lineage)
opts$met_nfeatures <- 2500    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.GpCs <- 3        # minimum number of GpC sites per feature
opts$acc_min.cells <- 5      # minimum number of cells per feature (and per stage_lineage)
opts$acc_nfeatures <- 2500    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 2500       # maximum number of genes (filter based on variance)

# Use only samples that passed QC for all omics (opts$include_all=FALSE)?
opts$include_all <- TRUE

# window length for the overlap between genes and features
opts$overlapGenes  <- FALSE
opts$gene_window  <- 5e4

# Define which cells to use
tmp1 <- fread(io$sample.metadata1) %>% 
  .[,stage_lineage:=paste(stage,lineage,sep="_")] %>%
  .[!is.na(id_met) & !is.na(id_acc)]
if (opts$include_all) { 
  opts$met_cells1 <- tmp1 %>% .[pass_metQC==T & outlier==F & stage_lineage%in%opts$stage_lineage1,id_met]
  opts$rna_cells1 <- tmp1 %>% .[pass_rnaQC==T & outlier==F & stage_lineage%in%opts$stage_lineage1,id_rna]
  opts$acc_cells1 <- tmp1 %>% .[pass_accQC==T & outlier==F & stage_lineage%in%opts$stage_lineage1,id_acc]
} else {
  opts$met_cells1 <- tmp1 %>% 
    .[pass_metQC==T & pass_rnaQC==T & pass_accQC==T & outlier==F & stage_lineage%in%opts$stage_lineage1,id_met]
  opts$rna_cells1 <- tmp1 %>% 
    .[pass_metQC==T & pass_rnaQC==T & pass_accQC==T & outlier==F & stage_lineage%in%opts$stage_lineage1,id_rna]
  opts$acc_cells1 <- tmp1 %>% 
    .[pass_metQC1==T & pass_rnaQC==T & pass_accQC==T & outlier==F & stage_lineage%in%opts$stage_lineage1,id_acc]
}

tmp2 <- fread(io$sample.metadata2)
opts$met_cells2 <- tmp2 %>% .[stage%in%opts$stage_lineage2,id_met]
opts$acc_cells2 <- tmp2 %>% .[stage%in%opts$stage_lineage2,id_acc]
