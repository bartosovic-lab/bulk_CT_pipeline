
args <- commandArgs(trailingOnly=TRUE)

out.file <- args[length(args)]
reports  <- args[-length(args)]


parse_picard <- function(report){
  cat(paste0("*** Reading file ",basename(report),"\n"))
  d <- read.table(file = report,
                  nrows = 2,
                  stringsAsFactors=FALSE,
                  sep = "\t")
  colnames(d) <- d[1,]
  d <- d[-1,]
  d$SAMPLE <- basename(report)
  return(d)
  
}

r.ls <- lapply(reports,parse_picard)


percent.duplication.ls <- lapply(r.ls,function(x){
   return(x[,c('SAMPLE','PERCENT_DUPLICATION','UNPAIRED_READS_EXAMINED','READ_PAIRS_EXAMINED')])
})


percent.duplication.df <- do.call('rbind',percent.duplication.ls)
write.csv(x=percent.duplication.df,
          file=out.file,
          quote=FALSE)
