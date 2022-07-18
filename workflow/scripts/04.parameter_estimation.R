# # for development
# save.image("rdev.RData")
# stop(
#   "saved rdev.RData; ########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################"
# )
# # setwd(
# #   "/Volumes/netscratch/dep_tsiantis/grp_laurent/struett/git_tsABC/tsABC/"
# # )
# # # setwd("/Users/abgushtdizi/Dropbox/professional/phd/git_wolbachia/wolbachia_abc/")
# # load(file = "rdev.RData")
# # snakemake@log$log1 = ""


library(locfit)
library("MASS")
library(abc)
library(arrow)


# Obtain parameters for estimate
{
  OUTFILE_ESTIMS  <- snakemake@output$estims
  INFILE_ABC  <- snakemake@input$sumstats
  INFILE_OBS  <- snakemake@input$observed
  INFILE_ID <- snakemake@input$identifier
  NTOLERATED <-
    as.numeric(snakemake@params$ntol)
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


# Read obs data
{
  # read in table
  df_obs <- data.frame(read.table(INFILE_OBS, header = T))
  
  
  # read obs index
  obsid <- read_feather(INFILE_ID)$population
  
  
  # subset stats (remove parameters)
  all_sumstats_obs <-
    df_obs[, grep("LinearCombination", names(df_obs))]
  sumstats_obs <-
    all_sumstats_obs[paste("LinearCombination", 0:(PLS_COMP - 1), sep = "_")]
  
  # read obs_params
  obs_params <- df_obs[, grep("param_", names(df_obs))]
  
  cat("  read and prepared obs and obs_params",
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
  for (i in 1:nrow(sumstats_obs)) {
    # extract summary statistic to estimate from
    my_target <- as.numeric(sumstats_obs[i, ])
    
    # extract true parameters
    true_params <- obs_params[i, ]
    
    
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
          "* Caught an estimation error on itertion ",
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
           obsid = obsid[i]
           )
    
    
    cat(
      "  saved priors and posteriors to list: ",
      i,
      " of ",
      nrow(sumstats_obs),
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

