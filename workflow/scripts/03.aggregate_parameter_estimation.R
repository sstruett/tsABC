# Aggregate parameter estimation
#
# Here, we aggregate the estimated parameters; we will provide the final output
# as a list which includes the dataframes of interest and the corresponding
# information: posteriors rej/adj; priors; put together with all wildcards


library(tidyverse)


# logfile
LOG <- snakemake@log$log1
cat("creating new log file\n", file = LOG, append = F)


# read in files
max_tolid <- 0
parameter.estimations <- list()
for (infileid in 1:length(snakemake@input$estims)) {
  infile.name <- snakemake@input$estims[infileid]
  parameter.estimations[[infileid]] <- readRDS(infile.name)
  parameter.estimations[[infileid]]$filename <- infile.name
  split <-
    strsplit(infile.name, split = "_|\\.|/")[[1]]  # read wc from filename
  parameter.estimations[[infileid]]$plsid <-
    as.numeric(split[which(split == "pls") + 1])
  parameter.estimations[[infileid]]$tolid <-
    as.numeric(split[which(split == "tolid") + 1])
  max_tolid <-
    max(c(max_tolid, as.numeric(split[which(split == "tolid") + 1])))
  parameter.estimations[[infileid]]$statcomposition <-
    as.numeric(split[which(split == "statcomp") + 1])
}


# log
cat("read data\n",
    file = LOG,
    append = T)

# loop through all estimates
df_collector <- list()
for (parestid in 1:length(parameter.estimations)) {
  parest <- parameter.estimations[[parestid]]
  
  # get values that are not posteriors; afterwards we can loop through the
  # posteriors
  parest.filename <- parest$filename
  parest.plsid <- parest$plsid
  parest.tolid <- parest$tolid
  parest.statcomposition <- parest$statcomposition
  
  # remove from list; afterwards we can loop through
  parest$filename <- NULL
  parest$plsid <- NULL
  parest$tolid <- NULL
  parest$statcomposition <- NULL
  parest$prior <- NULL
  
  
  # extract true params
  for (parest.id in 1:length(parest)) {
    true_params <- parest[[parest.id]]$true_params
    
    
    # loop through parameters and extract the posterior for each parameter
    for (param.name in colnames(true_params)) {
      post.rej <- parest[[parest.id]]$rej[, param.name]
      post.adj <- parest[[parest.id]]$adj[, param.name]
      
      stopifnot(parest.tolid == length(post.rej))
      stopifnot(parest.tolid == length(post.adj))
      
      
      df_0 <- data.frame(
        param = param.name,
        true_value = true_params[[param.name]],
        statcomposition = parest.statcomposition,
        pls = parest.plsid,
        tol = parest.tolid,
        regression = c("rej", "adj")
      )
      df_1 <-
        rbind.data.frame(c(post.rej, rep(NA, max_tolid - parest.tolid)),
                         c(post.adj, rep(NA, max_tolid - parest.tolid)))
      colnames(df_1) <- paste0("acc_", 1:parest.tolid)
      
      # collect df
      df_collector[[length(df_collector) + 1]] <-
        cbind.data.frame(df_0, df_1)
      
      # log
      cat(
        "collected data: ",
        parestid,
        "of",
        length(parameter.estimations),
        ";",
        parest.id,
        "of",
        length(parest),
        ";",
        param.name,
        "of",
        length(colnames(true_params)),
        ";",
        "\n",
        file = LOG,
        append = T
      )
      
      rm(df_0, df_1)
    }
  }
}
rm(parest, parameter.estimations, true_params)


# put into a single tibble
df <- do.call(rbind.data.frame, df_collector)

# log
cat("extracted and pasted data\n",
    file = LOG,
    append = T)


# add mode, mean, median
posteriors <- df %>%
  select(starts_with("acc_"))

df$mean <- apply(posteriors, 1, mean)
df$mode <- apply(posteriors, 1, function(x) {
  return(quantile(as.numeric(x), 0.5))
})
df$median <- apply(posteriors, 1, median)

# log
cat(
  "calculated mean, mode, median of each posterior\n",
  file = LOG,
  append = T
)


# rearrange columns
df <- df %>% select(!starts_with("acc_"))
df <- cbind.data.frame(df,
                       posteriors)


# save to file
saveRDS(object = df, file = snakemake@output$estims)


# log
cat("saved RDS",
    snakemake@output$estims,
    "\n",
    file = LOG,
    append = T)
