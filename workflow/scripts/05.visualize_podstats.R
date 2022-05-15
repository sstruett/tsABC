# This will plot the parameters of the pods and the correlation to the
# statistics; we will maximally consider the first 20 PLS components


library(ggplot2)
library(cowplot)
library(wesanderson)
library(tidyverse)


# theme set
theme_set(theme_cowplot())
logbreak <- sapply(10 ** (-10:10), function(x) {
  x * (1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10 ** (-10:10))
    return(as.character(log10(x)))
  else
    return("")
})


# small functions
find_max_pls  <- function(my_cols) {
  colslist <- strsplit(my_cols, "LinearCombination_")
  plsno <- integer(length = length(colslist))
  for (colid in 1:length(colslist)) {
    plsno[colid] <- as.numeric(colslist[[colid]][2])
  }
  return(max(plsno))
}


# log
LOG <- snakemake@log$log1
cat("creating new log file\n", file = LOG, append = F)


# read data; pseudo-observed data sets
pod.stats <- list()
for (infileid in 1:length(snakemake@input$transformed)) {
  infile.name <- snakemake@input$transformed[infileid]
  df_raw <- read.csv(infile.name, sep = "\t")
  
  
  # select all cols that arent stats and add identifiers
  nostat <- df_raw %>%
    select(!starts_with("LinearCombination_"))
  split <-
    strsplit(infile.name, split = "_|\\.|/")[[1]]  # read wc from filename
  nostat$statcomposition <-
    as.numeric(split[which(split == "statcomp") + 1])
  
  # add pod index
  nostat$podid <-
    as.numeric(interaction(nostat[, grep("^param_", colnames(nostat))]))
  
  
  # collect all cols that are stats and reduce or extend to 20 columns
  stats <- df_raw %>%
    select(starts_with("LinearCombination_"))
  if (ncol(stats) > 20) {
    wanted_cols <- paste("LinearCombination", 1:20 - 1, sep = "_")
    colindexes <- integer(length = 20)
    for (mycolid in 1:20) {
      mycolname <- wanted_cols[mycolid]
      colindexes[mycolid] <- which(mycolname == colnames(stats))
    }
    stats <- stats[, colindexes]
    
    # log
    cat("read data and reduced stats to 20 columns\n",
        file = LOG,
        append = T)
    
  } else if (ncol(stats) < 20) {
    max_pls <- find_max_pls(colnames(stats))
    extender <-
      setNames(data.frame(# LinearCombinations are zero-based
        matrix(
          ncol = 20 - max_pls - 1, nrow = nrow(stats)
        )),
        paste("LinearCombination", (max_pls + 1):(20 - 1), sep = "_"))
    stats <- cbind.data.frame(stats, extender)
    
    # log
    cat("read data and extended stats to 20 columns\n",
        file = LOG,
        append = T)
  }
  stopifnot(ncol(stats) == 20)
  
  
  # read into list
  pod.stats[[infileid]] <- cbind.data.frame(nostat, stats)
}


# log
cat("read data\n", file = LOG, append = T)
cat("note, pod indexes may differ from config file\n",
    file = LOG,
    append = T)


# put into a single tibble
df <- do.call(rbind.data.frame, pod.stats) %>% tibble()



# log
cat("start plotting\n", file = LOG, append = T)


# collect plots in list
plot_list <- list()
plot.id <- 0

# parameters
plot.id <- plot.id + 1
plot_list[[plot.id]] <- df %>%
  select(starts_with("param_"),
         starts_with("podid")) %>%
  pivot_longer(cols = -c("podid"),
               names_to = "param",
               values_to = "value") %>%
  ggplot(aes(podid, value, fill = param)) +
  geom_point(shape = 23) +
  facet_grid(param ~ ., scales = "free") +
  scale_y_continuous(trans = "log",
                     breaks = logbreak,
                     labels = round(log10(logbreak), 1)) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  theme(
    aspect.ratio = 0.707 ,
    legend.position = "none",
    panel.background = element_rect(fill = "gray99"),
    strip.background = element_blank()
  )

# log
cat("plotted parameters\n", file = LOG, append = T)


# param sumstat correlation
plot_list[[plot.id]] <-
  for (pls_name in colnames(df)[grep("^LinearCombination_", colnames(df))]) {
    plot.id <- plot.id + 1
    plot_list[[plot.id]] <- df %>%
      select(!starts_with("LinearCombination_"),
             starts_with(pls_name)) %>%
      rename(LinearCombination = pls_name) %>%
      pivot_longer(cols = colnames(df)[grep("^param_", colnames(df))],
                   names_to = "param",
                   values_to = "value") %>%
      ggplot(aes(podid, LinearCombination, fill = as.factor(statcomposition))) +
      geom_point(position = "jitter", shape = 23) +
      facet_grid(statcomposition ~ ., scales = "free") +
      scale_fill_manual(values = wes_palette("Darjeeling2")) +
      labs(title = paste(pls_name)) +
      theme(aspect.ratio = 0.707 / 2,
            legend.position = "none")
    
    
    # log
    cat("plotted",
        pls_name,
        "\n",
        file = LOG,
        append = T)
  }


# log
cat("start printing to file\n", file = LOG, append = T)


# print to file
pdf(snakemake@output$pdf)
for (my_plot in 1:length(plot_list)) {
  show(plot_list[[my_plot]])
  
  # log
  cat(
    "successfully printed plot to file:",
    my_plot,
    "of",
    length(plot_list),
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
