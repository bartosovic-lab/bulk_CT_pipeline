library(ggplot2)
library(ggthemes)

# Argument parsing
library(argparse)
parser <- ArgumentParser(description = "Aggregate Picard MarkDuplicates outputs into a single table and plot the data using R. The output files are saved as duplicates_table.tsv, sequencing_saturation.png, and library_size_plot.png.", usage = "Example usage: Rscript get_percent_duplication.R -i /path/to/picard_metrics/sample1_metrics.txt /path/to/picard_metrics/sample2_metrics.txt -o /path/to/output/")
parser$add_argument("-i", "--input", help="Picard MarkDuplicates output file", required=TRUE, nargs="+")
parser$add_argument("-o", "--output", help="Output folder", required=TRUE)
args <- parser$parse_args()


get_picard_metrics <- function(file){
  d <- read.table(file=file,skip=6,nrows=1,sep='\t',header=TRUE)
  return(d)
}

dupes.ls <- lapply(args$input, get_picard_metrics)
dupes.df <- data.frame(do.call(rbind, dupes.ls))
colnames(dupes.df) <- names(dupes.ls[[1]])

write.table(dupes.df, file=paste0(args$output, "duplicates_table.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

ggplot(data=dupes.df, aes(x=READ_PAIRS_EXAMINED, y=PERCENT_DUPLICATION)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Percent Duplication vs. Read Pairs Examined", x="Read Pairs Examined", y="Percent Duplication") +
  theme_base()
ggsave(paste0(args$output, "sequencing_saturation.png"))

ggplot(data=dupes.df, aes(x=READ_PAIRS_EXAMINED, y=ESTIMATED_LIBRARY_SIZE)) +
  geom_point() +
  labs(title="Estimated Library Size vs. Read Pairs Examined", x="Read Pairs Examined", y="Estimated Library Size") +
  theme_base() + coord_cartesian(ylim=c(0, max(dupes.df$ESTIMATED_LIBRARY_SIZE)*1.1))
ggsave(paste0(args$output,"library_size_plot.png"))