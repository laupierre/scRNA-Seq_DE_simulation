library (ggplot2)

res <- read.delim ("simulation_results.txt")

p <- ggplot (res, aes(x=method, y=percentage, fill=type)) +
     geom_bar (stat="identity",  position=position_dodge()) + facet_wrap (~cells + dropout) + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))
p <- p + coord_flip () + ggtitle ("Detection of DEG in simulated single-cell RNA-Seq") + geom_hline(yintercept=100, linetype = "dashed")
p
ggsave ("SCRIP simulation.pdf")

