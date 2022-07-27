# # for development
# save.image("rdev.RData")
# stop(
#   "saved rdev.RData; ########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################"
# )
# # setwd("/Users/struett/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/git_tsABC/tsABC/")
# # # setwd("/Users/abgushtdizi/Dropbox/professional/phd/git_wolbachia/wolbachia_abc/")
# # load(file="rdev.RData")
# # snakemake@log$log1=""


# log
LOG <- snakemake@log$log1
cat("creating new log file\n", file = LOG, append = F)


library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(wesanderson)


data <- snakemake@input$estims %>%
  readRDS %>%
  as_tibble %>%
  ungroup


# for parameters in the table get the prior boundaries
param_prior <- (data$param %>%
                  unique %>%
                  sort %>% strsplit("_") %>%
                  bind_cols())[2, ] %>%
  t %>%
  as_tibble %>%
  rename(paramid = V1)
boundaries <- apply(param_prior, 1, function(x) {
  # this function is maybe a bit ugly, but it translates the list of list, which
  # is provided only as a vector inside R from the snakemake object
  paramid <- as.numeric(x)
  priors <-
    matrix(
      snakemake@config$ABC$athaliana$priors,
      nrow = length(snakemake@config$ABC$athaliana$priors) / nrow(param_prior)
    ) %>% t %>% as_tibble
  
  
  return(priors[paramid + 1, 1:2] %>%
           rename(lower = V1,
                  upper = V2))
}) %>%
  bind_rows()
param_prior <- cbind(param_prior,
                     boundaries) %>%
  as_tibble %>%
  mutate(
    paramid = as.numeric(paramid),
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  )


plot.list <- list()
for (paramid in param_prior$paramid) {
  plot.list[[1 + length(plot.list)]] <- data %>%
    filter(param == paste0(c("param", paramid), collapse = "_")) %>%
    select(!starts_with("acc_")) %>%
    ggplot(aes(x = "regression", y = mode, fill = podid)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), shape = 16, size = 0.15)+
    facet_grid(.~regression) +
    scale_y_continuous(limits = param_prior[param_prior$paramid == paramid, 2:3] %>% as_vector %>% unname) +
    scale_fill_manual(values = wes_palette("Darjeeling1")) +
    theme(
      aspect.ratio = 1/0.707,
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "bottom"
      
      )
}


# print to file
pdf(snakemake@output$plot)
for (my_plot in 1:length(plot.list)) {
  show(plot.list[[my_plot]])
  
  # log
  cat(
    "successfully printed plot to file:",
    my_plot,
    "of",
    length(plot.list),
    "\n",
    file = LOG,
    append = T
  )
}
dev.off()

# log
cat("successfully printed plots to file\n",
    file = LOG,
    append = T)


# print csv
write_csv(data, file = snakemake@output$table)
