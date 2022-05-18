# visualize discretized bayes factors

library(ggplot2)
library(cowplot)
library(tidyverse)
library(wesanderson)


# set theme
theme_set(theme_cowplot())


# helper functions and values
discretize_bayes_factor <- function(bf, do.factor = T) {
  if (bf <= 10 ** 0)
    dbf = "Negative"
  else if (bf <= 10 ** 0.5)
    dbf = "Barely worth mentioning"
  else if (bf <= 10 ** 1)
    dbf = "Substantial"
  else if (bf <= 10 ** 1.5)
    dbf = "Strong"
  else if (bf <= 10 ** 2)
    dbf = "Very strong"
  else if (bf > 10 ** 2)
    dbf =  "Decisive"
  else
    stop(paste0(c("unknown bayes: ", as.character(bf)), collapse = ""))
  
  dbf <- factor(
    dbf,
    levels = c(
      "Negative",
      "Barely worth mentioning",
      "Substantial",
      "Strong",
      "Very strong",
      "Decisive"
    )
  )
  
  if (do.factor) {
    return(dbf)
  } else {
    return(as.numeric(dbf) - 1)
  }
}
logbreak <- sapply(10 ** (-10:10), function(x) {
  x * (1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10 ** (-10:10))
    return(as.character(log10(x)))
  else
    return("")
})


# log
LOG <- snakemake@log$log1
cat("creating new log file\n", file = LOG, append = F)


# create color palette
mcol <- wes_palette("Cavalcanti1")[c(5, 2)]
mcol <- c(mcol[1], "gray90", mcol[2])
colfunc <- colorRampPalette(mcol)
mcol <- colfunc(9)[c(4:9)]
mcol <- colfunc(9)[c(1, 5:9)]


# log
cat("creating color palette\n", file = LOG, append = T)


# read in data
df <- readRDS(snakemake@input$model_choice) %>%
  tibble() %>%
  pivot_longer(
    cols = c("rejection", "mnlogistic"),
    names_to = "regression_method",
    values_to = "bf"
  )

# log
cat("read data\n", file = LOG, append = T)


# find tsigma per podid
tsigma <- snakemake@params$tsigma_per_podid %>% as.numeric()
df$tsigma <- sapply(df$podid, function(x)
  return(tsigma[x]))

# discretize bf
df$bf_discrete <- sapply(df$bf, discretize_bayes_factor)

# log
cat("discretized bayes factor\n",
    file = LOG,
    append = T)

# create empty plot list
plot_list <- list()
plot_index <- 0


# make as bar plot per podid
plot_index <- plot_index + 1
plot_list[[plot_index]] <- df %>%
  ggplot(aes(podid, fill = bf_discrete)) +
  geom_bar(position = "fill") +
  facet_grid(statcomposition ~ regression_method + pls) +
  scale_fill_manual(values = mcol) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


# calculate percentage of discrete bayes factor
df <- df %>%
  group_by(podid, statcomposition, pls, tolid, regression_method, tsigma) %>%
  mutate(N = n()) %>%
  group_by(podid,
           statcomposition,
           pls,
           tolid,
           regression_method,
           tsigma,
           bf_discrete) %>%
  mutate(n = n()) %>%
  mutate(bf_proportion = n / N) %>%
  ungroup() %>%
  select(-c(bf, N, n)) %>%
  distinct()


# check if they sum up to 1
a <- df %>%
  group_by(podid, statcomposition, pls, tolid, regression_method, tsigma) %>%
  summarise(total_proportion = sum(bf_proportion))
stopifnot(all(a$total_proportion == 1))


# remove column
df <- df %>% select(-podid)  # no podid needed as we use tsigma


# add zeros for proper plotting
zcounter = 0
for (a in unique(df$statcomposition)) {
  for (b in unique(df$regression_method)) {
    for (c in unique(df$tsigma)) {
      for (d in unique(df$bf_discrete)) {
        for (e in unique(df$tolid)) {
          for (f in unique(df$pls)) {
            sdf <- df %>%
              subset(statcomposition == a) %>%
              subset(regression_method == b) %>%
              subset(tsigma == c) %>%
              subset(bf_discrete == d)
            n <-  sdf %>%
              nrow()
            if (n == 0) {
              zcounter = zcounter + 1
              df = rbind.data.frame(
                df,
                data.frame(
                  statcomposition = a,
                  pls = f,
                  tolid = e,
                  regression_method = b,
                  tsigma = c,
                  bf_discrete = d,
                  bf_proportion = 0
                )
              )
            }
            
          }
        }
      }
    }
  }
}


# log
cat(
  "calculated proportion per group and added zeros\n",
  file = LOG,
  append = T
)
cat("there were",
    zcounter,
    " zeros to be added\n",
    file = LOG,
    append = T)


# create plot as proportional area plot over time
for (a in unique(df$statcomposition)) {
  for (b in unique(df$pls)) {
    for (c in unique(df$tolid)) {
      for (d in unique(df$regression_method)) {
        sdf <- df %>%
          subset(statcomposition == a) %>%
          subset(pls == b) %>%
          subset(tolid == c) %>%
          subset(regression_method == d)
        
        if (nrow(sdf) > 0) {
          plot_index <- plot_index + 1
          plot_list[[plot_index]] <- sdf %>%
            ggplot(aes(
              x = tsigma,
              y = bf_proportion * 100,
              fill = bf_discrete
            )) +
            geom_area() +
            geom_hline(yintercept = 5,
                       col = "white",
                       size = 0.1) +
            geom_hline(yintercept = 20,
                       col = "white",
                       size = 0.1) +
            geom_hline(yintercept = 50,
                       col = "white",
                       size = 0.1) +
            geom_hline(yintercept = 80,
                       col = "white",
                       size = 0.1) +
            geom_hline(yintercept = 95,
                       col = "white",
                       size = 0.1) +
            geom_vline(xintercept = 200000,
                       col = "white",
                       size = 0.1) +
            # facet_grid(regression_method+pls~statcomposition)+
            scale_x_continuous(trans = "log",
                               breaks = logbreak,
                               labels = loglabel) +
            scale_y_continuous(breaks = c(5, 95),
                               labels = c("0%", "100%")) +
            scale_fill_manual(values = mcol) +
            theme(
              # axis.text.x = element_text(angle = 60, hjust = 1),
              aspect.ratio = 1,
              strip.background = element_blank(),
              axis.title.y = element_blank(),
              # axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks.y = element_blank(),
              axis.ticks.length.x = unit(-0.1, "lines"),
              # panel.spacing = unit(-0.9, "lines"),
              # panel.border = element_rect(colour = "black")
            ) +
            labs(
              x = "time after transition [log10]",
              fill = "model support",
              title = paste0(
                "statcomposition: ",
                a,
                "\npls: ",
                b,
                "\ntolid: ",
                c,
                "\nregression_method: ",
                d,
                collapse = ""
              )
            )
        }
      }
    }
  }
}




# log
cat("created plots\n",
    file = LOG,
    append = T)


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
