library(locfit)
library("MASS")
library(abc)


# Obtain parameters for estimate
{
  OUTFILE_ESTIMS  <- snakemake@output$estims
  INFILE_ABC  <- snakemake@input$sumstats
  INFILE_POD  <- snakemake@input$podstats
  INFILE_PODID <- snakemake@input$podid
  NTOLERATED <-
    as.numeric(snakemake@params$ntol)
  STATCOMP_NO <- as.numeric(snakemake@params$statcomposition_no)
  STATCOMP <- snakemake@params$statcompositions[STATCOMP_NO]
  PLS_COMP <- as.numeric(snakemake@params$pls)
  PLOT_DIR <- snakemake@output$postplots
  REGRESSION <- snakemake@params$regression
  LOG <- snakemake@log$log1
}

cat("creating new log file\n", file = LOG, append = F)

# Read abc sim data and subset to pls
{
  # read in table
  df_abc <- data.frame(read.table(INFILE_ABC, header = T))
  
  
  # subset stats (remove parameters)
  all_sumstats_abc <-
    df_abc[, grep("LinearCombination", names(df_abc))]
  if (ncol(all_sumstats_abc) <  PLS_COMP) {
    cat(
      "  changed PLS set from ",
      PLS_COMP,
      " to ",
      ncol(all_sumstats_abc),
      "\n",
      file = LOG,
      append = T
    )
    PLS_COMP <- ncol(all_sumstats_abc)
  }
  sumstats_abc <-
    all_sumstats_abc[, paste("LinearCombination", 0:(PLS_COMP - 1), sep = "_")]
  
  
  # Get params and remove the columns that do not have variance
  params_abc_all <- df_abc[, grep("param", names(df_abc))]
  params_abc <- params_abc_all[, apply(params_abc_all, 2, var) != 0]
  
  cat("  read and prepared simulations",
      "\n",
      file = LOG,
      append = T)
}


# Read pod data
{
  # read in table
  df_pod <- data.frame(read.table(INFILE_POD, header = T))
  
  
  # read pod index
  podid <- read.csv(INFILE_PODID, header = FALSE)[[1]] + 1  # podid are 0-based
  
  
  # subset stats (remove parameters)
  all_sumstats_pods <-
    df_pod[, grep("LinearCombination", names(df_pod))]
  sumstats_pods <-
    all_sumstats_pods[paste("LinearCombination", 0:(PLS_COMP - 1), sep = "_")]
  
  # read pods_params
  pod_params <- df_pod[, grep("param_", names(df_pod))]
  
  cat("  read and prepared pods and pod_params",
      "\n",
      file = LOG,
      append = T)
}


# Produce logit boundaries from actual priors
{
  logit_boundaries <- t(apply(params_abc, 2, range))
  colnames(logit_boundaries) <- c("minimal", "maximal")
}

# results lists
{
  prior_and_posteriors <- list()
}


{
  # parameter preparation
  my_tolerance <- NTOLERATED / nrow(sumstats_abc)
}


# make estimates with library abc
{
  # make an estimate for each set of sumstats in the table (i. e. the samples)
  for (i in 1:nrow(sumstats_pods)) {
    # extract summary statistic to estimate from
    my_target <- as.numeric(sumstats_pods[i, ])
    
    # extract true parameters
    true_params <- pod_params[i, ]
    
    
    # create estimate
    flag <- TRUE
    tryCatch(
      expr = {
        abc_result <- abc(
          target = my_target,
          param = params_abc,
          sumstat = sumstats_abc,
          tol = my_tolerance,
          method = REGRESSION,
          MaxNWts = 5000,
          transf = "logit",
          logit.bounds = logit_boundaries
        )
      },
      error = function(e) {
        cat(
          "* Caught a plotting error on itertion ",
          i,
          "\n",
          file = LOG,
          append = T
        )
        flag <- FALSE
      }
    )
    if (!exists("abc_result"))
      next
    
    
    # clean the data from na's
    rejection_values <- na.omit(abc_result$unadj.values)
    adjusted_values <- na.omit(abc_result$adj.values)
    
    
    # save results into list
    prior_and_posteriors[[i]] <-
      list(rej = rejection_values,
           adj = adjusted_values,
           true_params = true_params,
           podid = podid[i]
           )
    
    
    
    cat(
      "  saved priors and posteriors to list: ",
      i,
      " of ",
      nrow(sumstats_pods),
      "\n",
      file = LOG,
      append = T
    )
    
    # plot diagnostics
    tryCatch(
      expr = {
        dir.create(PLOT_DIR, showWarnings = FALSE)
        pdf(paste0(PLOT_DIR, "/diagnostic_plots_", i, ".pdf", collapse = ""))
        plot(abc_result, param = params_abc, ask = F)
        dev.off()
      },
      error = function(e) {
        cat(
          "* Caught a plotting error on itertion ",
          i,
          "\n",
          file = LOG,
          append = T
        )
      }
    )
  }
}


# save results list
{
  # add prior to the result list
  prior_and_posteriors[["prior"]] <- params_abc
  
  saveRDS(prior_and_posteriors, file = OUTFILE_ESTIMS)
}
