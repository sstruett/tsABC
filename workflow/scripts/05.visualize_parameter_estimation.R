# visualize parameter estimations using average quantiles


library(ggplot2)
library(cowplot)
library(tidyverse)
library(wesanderson)


# set theme
theme_set(theme_cowplot())


# log
LOG <- snakemake@log$log1
cat("creating new log file\n", file = LOG, append = F)


# read parameter estimation file
data <- readRDS(snakemake@input$estimations)


# get the variable parameter
# log
cat(
  "for this plotting we make the assumption that only param_3 (t_sigma) is varying\n",
  file = LOG,
  append = T
)
cat(
  "script will run also otherwise, but plots won't be useful\n",
  file = LOG,
  append = T
)


# define interquantile ranges
interquantile.ranges <-
  sort(unique(as.numeric(as.character(
    snakemake@params$iqr
  ))))


# create color palette
mcol <- wes_palette("Darjeeling1")[1:2]
mcol <- c(mcol[1], mcol[2])
colfunc <- colorRampPalette(mcol)
mcol <- colfunc(length(interquantile.ranges))
if (0 %in% interquantile.ranges) {
  mcol <- c(mcol, rev(mcol)[2:length(mcol)])
} else {
  mcol <- c(mcol, rev(mcol))
}



# helper functions
calc_quants <-
  function(posterior.data.frame,
           interquantile.ranges = c(0.99, 0.95, 0.9, 0.8, 0.5, 0.25, 0.1, 0)) {
    # Computes the quantiles row wise
    #
    # Args:
    #   posterior.data.frame: data.frame; should contain a single explicit
    #     distribution per row
    #   interquantile.ranges: numeric vector; which interquantiles to calculate,
    #     e. g. 0.99 calculates the 0.005-th and 0.995-th quantiles; e. g. 0
    #     calculates the mode
    #
    # Returns:
    #   data.frame with quantiles
    iqr <- sort(interquantile.ranges)
    stopifnot(all(iqr >= 0))
    
    qr <- sort(unlist(sapply(iqr, function(x) {
      return(unique(c(0.5 - x / 2, 0.5 + x / 2)))
    })))
    
    apply.f <- function(x) {
      quantile(x, probs = qr, na.rm = T)
    }
    
    return(t(apply(posterior.data.frame, 1, apply.f)))
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


# calculate quantiles of posteriors; remove accepted posteriors
data.quantiles <-
  calc_quants(data %>% select(starts_with("acc_")), interquantile.ranges = interquantile.ranges)
data.index <- data %>% select(!starts_with("acc_"))
data <- cbind.data.frame(data.index,
                         data.quantiles) %>% tibble()


# find tsigma per podid
tsigma <- snakemake@params$tsigma_per_podid %>% as.numeric()
data$tsigma <- sapply(data$podid, function(x)
  return(tsigma[x]))


# create empty plot list
plot_list <- list()
plot_index <- 0


# create plots
msize <- 2


for (a in unique(data$param)) {
  for (b in unique(data$statcomposition)) {
    for (c in unique(data$pls)) {
      for (d in unique(data$tol)) {
        for (e in unique(data$regression)) {
          sdf <- data %>%
            subset(param == a) %>%
            subset(statcomposition == b) %>%
            subset(pls == c) %>%
            subset(tol == d) %>%
            subset(regression == e)
          
          if (nrow(sdf) > 0) {
            sdf <- sdf %>%
              pivot_longer(
                -c(
                  param,
                  true_value,
                  statcomposition,
                  pls,
                  tol,
                  regression,
                  podid,
                  mean,
                  mode,
                  median,
                  tsigma
                ),
                names_to = "quantile",
                values_to = "value"
              ) %>%
              group_by(
                param,
                true_value,
                statcomposition,
                pls,
                tol,
                regression,
                podid,
                quantile,
                tsigma
              ) %>%
              summarise(mean.iqr = mean(value))
            
            sdf$quantile <- as.numeric(sub("%", "", sdf$quantile))
            mlevels <- sort(unique(sdf$quantile))
            sdf$quantile <- factor(sdf$quantile, levels = mlevels)
            
            mxlim <- range(sdf$tsigma)
            mylim <- c(NA, NA)
            
            plot_index <- plot_index + 1
            plot_list[[plot_index]] <- sdf %>%
              ggplot(aes(
                x = tsigma,
                y = mean.iqr,
                col = quantile
              )) +
              geom_line(aes(tsigma, true_value),
                        col = "black",
                        size = msize) +
              geom_line(show.legend = F, size = msize) +
              # geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.02, show.legend = FALSE)+
              # geom_text() %>%
              # facet_grid(plsComp ~ sscomp) +
              scale_x_continuous(
                trans = "log10",
                limits = mxlim,
                breaks = logbreak,
                labels = loglabel
              ) +
              scale_color_manual(values = mcol) +
              theme(
                legend.position = "none",
                aspect.ratio = 1,
                panel.border = element_rect(colour = "black", size = msize),
                # text = element_text(size = 12),
                #   strip.background = element_blank(),
                # axis.title.y = element_blank(),
                # axis.text = element_blank(),
                # axis.text.y = element_blank(),
                axis.line = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_line(size = msize * 0.8),
                # axis.ticks.y = element_blank(),
                # axis.ticks.x = element_line(size = 1),
                # axis.ticks.length.x = unit(-1, "lines"),
                # panel.spacing = unit(-0.9, "lines")
              ) +
              labs(
                title = paste0(
                  "param: ",
                  a,
                  "\nstatcomposition: ",
                  b,
                  "\npls: ",
                  c,
                  "\ntol: ",
                  d,
                  "\nregression_method: ",
                  e,
                  collapse = ""
                )
              )
            
            # make log scale for some parameters
            if (a %in% c("param_0", "param_3"))
              plot_list[[plot_index]] <- plot_list[[plot_index]] +
              scale_y_continuous(
                trans = "log10",
                limits = mylim,
                breaks = logbreak,
                labels = loglabel
              )
          }
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
