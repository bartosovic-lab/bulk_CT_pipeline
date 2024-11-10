library(ggplot2)

# Description: This script reads the Picard MarkDuplicates output file and returns the number of read pairs examined and the percent duplication.
options <- commandArgs(trailingOnly = TRUE)
if (length(options) == 0) {
  stop("Usage: Rscript get_percent_duplication.R <Picard MarkDuplicates output file>")
}

get_picard_metrics <- function(file){
  d <- read.table(file=file,skip=6,nrows=1,sep='\t',header=TRUE)
  return(d)
}

dupes.ls <- lapply(options, get_picard_metrics)
dupes.df <- data.frame(do.call(rbind, dupes.ls))
colnames(dupes.df) <- names(dupes.ls[[1]])

ggplot(data=dupes.df, aes(x=READ_PAIRS_EXAMINED, y=PERCENT_DUPLICATION)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Percent Duplication vs. Read Pairs Examined", x="Read Pairs Examined", y="Percent Duplication") +
  theme_minimal()
ggsave("duplication_plot.png")

ggplot(data=dupes.df, aes(x=READ_PAIRS_EXAMINED, y=ESTIMATED_LIBRARY_SIZE)) +
  geom_point() +
  labs(title="Estimated Library Size vs. Read Pairs Examined", x="Read Pairs Examined", y="Estimated Library Size") +
  theme_minimal() + coord_cartesian(ylim=c(0, max(dupes.df$ESTIMATED_LIBRARY_SIZE)*1.1))
ggsave("library_size_plot.png")