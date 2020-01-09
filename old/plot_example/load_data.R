suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(scater))
suppressMessages(library(purrr))

# Initialize argument parser
# p <- ArgumentParser(description='')
# p$add_argument('-rna.id','--rna.id', type="character", help='Gene. RNA expression')
# p$add_argument('-met.id','--met.id', type="character", help='Feature. DNA methylation')
# p$add_argument('-met.anno','--met.anno', type="character", help='Annotation. DNA methylation')
# p$add_argument('-acc.id','--acc.id', type="character", help='Feature. Chromatin accessibility')
# p$add_argument('-acc.anno','--acc.anno', type="character", help='Annotation. Chromatin accessibility')

# Read arguments
# args <- p$parse_args(commandArgs(TRUE))

load_data <- function(rna.id, met.id, met.anno, acc.id, acc.anno, min.cpg=1, min.gpc=1) {
  
  ################
  ## Define I/O ##
  ################
  
  io <- list()
  # io$sample.metadata <- "/Users/ricard/data/gastrulation/sample_metadata_scNMT.txt"
  io$met.dir <- "/Users/ricard/data/gastrulation/met/parsed"
  io$acc.dir <- "/Users/ricard/data/gastrulation/acc/parsed"
  io$rna.file <- "/Users/ricard/data/gastrulation/rna/parsed/SingleCellExperiment.rds"
  # io$annos_dir <- "/Users/ricard/data/gastrulation/features/filt"
  # io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  
  ##########
  ## Load ##
  ##########
  
  # Load DNA methylation data
  met_dt <- fread(sprintf("zcat < %s/%s.tsv.gz",io$met.dir,met.anno), stringsAsFactors=F, quote="")
  colnames(met_dt) <- c("id_met","id","anno","Nmet","N","rate")
  
  # Load DNA accessibility data
  acc_dt <- fread(sprintf("zcat < %s/%s.tsv.gz",io$acc.dir,acc.anno), stringsAsFactors=F, quote="")
  colnames(acc_dt) <- c("id_acc","id","anno","Nmet","N","rate")
  
  # Load RNA data
  sce <- readRDS(io$rna.file)
  rna_dt <- exprs(sce) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
    melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
    merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))
  
  ############
  ## Filter ##
  ############
  
  # Select ID
  met_dt <- met_dt[id%in%met.id] %>% setnames("rate","value")
  acc_dt <- acc_dt[id%in%acc.id] %>% setnames("rate","value")
  rna_dt <- rna_dt[ens_id%in%rna.id] %>% setnames("expr","value")
  
  # Filter by coverage
  met_dt <- met_dt[N>=min.cpg]
  acc_dt <- acc_dt[N>=min.gpc]
  
  return(list("met"=met_dt, "acc"=acc_dt, "rna"=rna_dt))
}