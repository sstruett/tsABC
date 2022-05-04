library(arrow)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyverse)


theme_set(
  theme_cowplot() +
    theme(axis.text.x = element_blank(),
          legend.position = "none")
)



setwd(
  "/Users/struett/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/git_tsABC/tsABC/"
)


df.nomask <- read_feather("results/athal/observed/population_Relicts.masked.feather")
# df.nomask <- read_feather("results/abc/pods/podstats.masked.feather")


df <- df.nomask[,which(apply(df.nomask, 2, var) != 0)]
df$index <- 1:nrow(df)


df %>%
  select(starts_with("sfs_") | starts_with("index")) %>%
  pivot_longer(., -index) %>%
  ggplot(aes(name, value, col=as.factor(index))) +
  geom_jitter(width = 0.2)


df %>%
  select(starts_with("ld_") | starts_with("index")) %>%
  pivot_longer(., -index) %>%
  ggplot(aes(name, value, col=as.factor(index))) +
  geom_jitter(width = 0.2)

tm_win <- df %>% select(starts_with("tm_win"))
tm_win_mean <- apply(tm_win, 2, mean)
tm_win <- tm_win[, order(tm_win_mean, decreasing = T)]
tm_win$index <- df$index

tm_win_long <- tm_win %>%
  pivot_longer(., -index)
tm_win_long$name <- factor(tm_win_long$name, levels = colnames(tm_win))

tm_win_long %>%
  ggplot(aes(name, value, col=as.factor(index))) +
  geom_jitter(width = 0.2)


tm_win <- df %>% select(starts_with("tm_win"))
matrix(apply(tm_win, 2, mean), nrow = sqrt(ncol(tm_win))) %>%
  image(zlim=c(0.0, 0.2))


