# Aggregate model choices


# logfile
LOG <- snakemake@log$log1
cat("creating new log file\n", file = LOG, append = F)


# read data
bayes_factors <- list()
for (infileid in 1:length(snakemake@input$bayes_factors)) {
  infile.name <- snakemake@input$bayes_factors[infileid]
  bayes_factors[[infileid]] <- readRDS(infile.name)
  split <-
    strsplit(infile.name, split = "_|\\.|/")[[1]]  # read wc from filename
  bayes_factors[[infileid]]$statcomposition <-
    rep(as.numeric(split[which(split == "statcomp") + 1]), length(bayes_factors[[infileid]]$podid))
  bayes_factors[[infileid]]$pls <-
    rep(as.numeric(split[which(split == "pls") + 1]), length(bayes_factors[[infileid]]$podid))
  bayes_factors[[infileid]]$tolid <-
    rep(as.numeric(split[which(split == "tolid") + 1]), length(bayes_factors[[infileid]]$podid))
}
bayes_factors <- do.call(rbind.data.frame, bayes_factors)


# log
cat("aggregated data\n",
    file = LOG,
    append = T)


# save file
saveRDS(bayes_factors, snakemake@output$bayes_factors)


# log
cat(
  "saved RDS",
  snakemake@output$bayes_factors,
  "\n",
  file = LOG,
  append = T
)
