#####################
## Define settings ##
#####################

source("/Users/ricard/Guo2017_preimplantation/settings.R")

# Define I/O
# io$mm10.genome <- "/Users/ricard/data/mm10_sequence/mm10.genome"
# io$stats <- paste0(io$basedir,"/met/stats/stats_per_chromosome.txt.gz")
io$outdir <- paste0(io$basedir,"/met/stats/pdf")

# Define options
# opts$chr <- c(paste0("chr",1:19),"X")
opts$chr <- c(1:19,"X")

# Update sample metadata
sample_metadata <- sample_metadata %>% 
  .[!is.na(id_met)] %>%
  droplevels

############################
## Load precomputed stats ##
############################

stats <- fread(io$met.stats_per_chr) %>%
  .[chr%in%opts$chr] %>%
  .[,chr:=factor(chr,levels=opts$chr)] %>%
  # .[,mean_coverage:=mean(coverage),by="id_met"] %>%
  # .[,relative_coverage:=coverage/mean_coverage,by="id_met"] %>%
  merge(sample_metadata, by="id_met")

# Sanity check
# all(sort(unique(stats$id_met)) == sort(sample_metadata$id_met))

# Load chromosome length
mm10.genome <- fread(io$mm10.genome) %>%
  setnames(c("chr","chr_length")) %>%
  .[chr%in%opts$chr] %>% .[,chr:=factor(chr,levels=opts$chr)]

####################################
## Plot relative coverage per chr ##
####################################

to.plot <- stats %>%
  .[,.(coverage=as.double(sum(coverage))), by=c("embryo","chr")] %>%
  merge(mm10.genome, by="chr") %>% 
  .[,coverage:=coverage/as.double(chr_length), by="embryo"]# %>%
  # .[,norm_coverage:=coverage/mean(coverage),by="embryo"]

ggscatter(to.plot, x="embryo", y="coverage") +
  # geom_hline(yintercept=1, linetype="dashed") +
  facet_wrap(~chr, scales="fixed") +
  labs(x="", y="normalised DNAm coverage") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

###################################
## Plot methylation rate per chr ##
###################################

to.plot <- stats %>%
  .[sex%in%c("Female","Male")] %>% droplevels %>% 
  .[,.(rate=mean(mean)), by=c("embryo","chr","sex","stage")]

sample_metadata[,.N,by=c("embryo","stage")]

for (i in unique(to.plot$stage)) {
  # p <- ggbarplot(to.plot[stage==i], x="sex", y="rate", fill="sex", stat="identity") +
  p <- ggplot(to.plot[stage==i], aes_string(x="embryo", group="embryo", y="rate", fill="sex")) +
    geom_bar(stat="identity", color="black") +
    # geom_hline(yintercept=1, linetype="dashed") +
    facet_wrap(~chr, scales="fixed") +
    labs(x="", y="Global DNA methylation (%)") +
    theme_classic() +
    theme(
      axis.text.y = element_text(size=rel(0.75)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  pdf(sprintf("%s/methylation_per_chr_%s.pdf",io$outdir,i), width=8, height=10)
  # png(sprintf("%s/barplots_%s.png",io$outdir,i), width=6, height=4, units="in", res=400)
  print(p)
  dev.off()
}
