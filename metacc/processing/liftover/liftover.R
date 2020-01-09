library(data.table)
library(purrr)
library(liftOver)
library(GenomicRanges)

io <- list()
# io$indir <- "/Users/ricard/data/Guo_2017/acc/raw"
# io$outdir <- "/Users/ricard/data/Guo_2017/acc/raw/liftover"; dir.create(io$outdir)

io$indir <- "/hps/nobackup/stegle/users/ricard/Guo_2017/acc/not_used/ESC"
io$outdir <- "/hps/nobackup/stegle/users/ricard/Guo_2017/acc/not_used/ESC/liftover"; dir.create(io$outdir)


# Liftover from mm9 to mm10
# chain = import.chain("/Users/ricard/data/mm10_sequence/mm9ToMm10.over.chain")
chain = import.chain("/hps/nobackup/stegle/users/ricard/mm10/sequence/mm9ToMm10.over.chain")

files <- list.files(io$indir, pattern=".tsv.gz$", full.names=FALSE)
files <- tools::file_path_sans_ext(files)

for (i in files) {
  print(i)
  data <- fread(cmd=sprintf("zcat < %s/%s.gz",io$indir, i)) %>%
    # setnames(c("chr","start","end","rate")) %>%
    setnames(c("chr","start","rate")) %>% .[,end:=start] %>%

    makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE) %>%
    liftOver(chain) %>%
    as.data.frame() %>%
    setDT()  %>% 
    .[, .(seqnames, start, rate)] %>% setnames("seqnames","chr") %>%
    .[, chr := gsub("chr", "", chr)] %>%
    setkey(chr,start)
  fwrite(data, sprintf("%s/%s",io$outdir,i), sep="\t", quote=FALSE, col.names=FALSE, na = "NA")
}

system(sprintf("pigz -p 4 %s/*.tsv",io$outdir))