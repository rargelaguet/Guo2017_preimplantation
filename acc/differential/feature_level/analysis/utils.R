gg_barplot <- function(tmp, title = "", ylim=NULL) {
  
  if (is.null(ylim)) {
    ylim <- c(min(tmp$value, na.rm=T), max(tmp$value, na.rm=T))
  }
  
  p <- ggplot(tmp, aes(x=anno, y=value, group=anno)) +
    geom_bar(aes(fill=anno), color="black", stat="identity", position="dodge") +
    # scale_fill_manual(values=opts$colors) +
    geom_hline(yintercept=0, color="black") +
    scale_y_continuous(limits=c(ylim[1],ylim[2])) +
    labs(title="", x="", y="Number of hits") +
    theme_classic() +
    theme(
      axis.text = element_text(size=rel(1.0), color='black'),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size=rel(1.0), color='black'),
      axis.line = element_line(color="black"),
      legend.position = "right",
    )
  
  return(p)
}

gg_volcano_plot <- function(tmp, top_genes=10, xlim=NULL, ylim=NULL) {
  
  negative_hits <- tmp[sig==TRUE & diff<0,id]
  positive_hits <- tmp[sig==TRUE & diff>0,id]
  all <- nrow(tmp)
  
  if (is.null(xlim))
    xlim <- max(abs(tmp$diff), na.rm=T)
  if (is.null(ylim))
    ylim <- max(-log10(tmp$p.value), na.rm=T)
  
  p <- ggplot(tmp, aes(x=diff, y=-log10(p.value))) +
    labs(title="", x="Differential levels (%)", y=expression(paste("-log"[10],"(p.value)"))) +
    # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
    geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange") +
    ggrastr::geom_point_rast(aes(color=sig), size=1) +
    scale_color_manual(values=c("black","red")) +
    scale_x_continuous(limits=c(-xlim-2,xlim+2)) +
    scale_y_continuous(limits=c(0,ylim+1)) +
    annotate("text", x=0, y=ylim+1, size=7, label=sprintf("(%d)", all)) +
    annotate("text", x=-75, y=ylim+1, size=7, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=75, y=ylim+1, size=7, label=sprintf("%d (+)",length(positive_hits))) +
    ggrepel::geom_text_repel(data=head(tmp[sig==T],n=top_genes), aes(x=diff, y=-log10(p.value), label=id), size=4) +
    theme_classic() +
    theme(
      axis.text = element_text(size=rel(1.25), color='black'),
      axis.title = element_text(size=rel(1.5), color='black'),
      # axis.title.y = element_text(margin=margin(0,10,0,0)),
      # axis.title.x = element_text(margin=margin(10,0,0,0)),
      legend.position = "none"
    )
  
  return(p)
}

