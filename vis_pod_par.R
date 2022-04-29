library(arrow)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)


setwd(
  "/Users/struett/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/git_tsABC/tsABC/"
)


df.nomask <- read_feather("results/abc/pods/podstats.feather")
df.masked <-
  read_feather("results/abc/pods/podstats.masked.feather")



# start a function here
plot_list_per_stat <- function(df) {
  pods <- unique(interaction(df[, grep("^param_", colnames(df))]))
  npods <-
    length(unique(interaction(df[, grep("^param_", colnames(df))])))
  
  pod_plot_list <- list()
  
  for (podid in 1:npods) {
    df.pod <-
      df[pods[podid] == interaction(df[, grep("^param_", colnames(df))]), ]
    params <- df.pod[, grep("^param_", colnames(df))]
    stopifnot(all(apply(params, 2, function(x) {
      length(unique(x)) == 1
    })))
    
    # print params:
    for (parid in 1:ncol(params)) {
      cat(colnames(params)[parid], ": ", as.numeric(params[1, parid]), "\n")
    }
    
    # visualize the distribution of stats
    stats <- df.pod[, grep("^param_", colnames(df), invert = T)]
    stats <- stats[, grep("podid", colnames(stats), invert = T)]
    stats.names <- colnames(stats)
    stats.names <- lapply(strsplit(stats.names, "_"), function(x) {
      l <- length(x)
      return(paste(x[1:(l - 1)], collapse = "_"))
    }) %>% unlist() %>% unique()
    
    # start function/ for loop
    
    plot_list <- list()
    for (statid in 1:length(stats.names)) {
      stats.set <- stats[, grep(stats.names[statid], colnames(stats))]
      stats.set <-
        stats.set[, sapply(colnames(stats.set), function(x) {
          var(stats.set[[x]]) != 0
        })]
      stats.set.long <- stats.set %>%
        pivot_longer(cols = colnames(.),
                     names_to = "statid",
                     values_to = "count")
      
      stats.set.long$statid <-
        sapply(stats.set.long$statid, function(x) {
          return(as.numeric(rev(strsplit(x, "_")[[1]])[1]))
        })
      
      
      plot_list[[statid]] <- stats.set.long %>%
        ggplot(aes(as.factor(statid), count)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(
          position = position_jitter(width = 0.3),
          col = "gray",
          shape = 16,
          size = 0.5
        )
    }
    
    pod_plot_list[[podid]] <- plot_list
    
  }
  return(pod_plot_list)
}


pod_plot_list <- plot_list_per_stat(df.nomask)

pod_plot_list_masked <- plot_list_per_stat(df.masked)

which_stat <- 2; {
  new_list <- list()
  for (i in 1:length(pod_plot_list)) {
    new_list[[i]] <- pod_plot_list[[i]][[which_stat]]
  }
  plot_grid(plotlist = new_list, align = T, nrow = 3)
}





