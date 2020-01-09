to.plot <-stats %>% melt(
  id.vars = c("sample","chr"), 
  measure.vars = list("rate"=c("rate_forward", "rate_reverse") , "coverage"=c("coverage_forward", "coverage_reverse")),
  variable.name = c("strand"),
  variable.factor = TRUE
)

ggplot(to.plot, aes(x=strand, y=rate)) +
  facet_wrap(~chr) +
  # geom_bar(stat="identity", position="dodge") +
  geom_boxplot() +
  theme_classic()


