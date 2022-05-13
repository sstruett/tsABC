library(arrow)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyverse)
library(wesanderson)
library(reshape)

# plot stats from observed athal; files are hardcoded


theme_set(
  theme_cowplot() +
    theme(
      aspect.ratio = 0.707,
      legend.position = "right",
      panel.background = element_rect(fill = "gray99")
    )
)

logbreak <- sapply(10 ** (-10:10), function(x) {
  x * (1:10)
}) %>% unique() %>% sort()
loglabel <- sapply(logbreak, function(x) {
  if (x %in% 10 ** (-10:10))
    return(as.character(log10(x)))
  else
    return("")
})


# setwd(
#   "/Users/struett/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/git_tsABC/tsABC/"
# )


# prepare files
pattern_from_vector <- function(my.list, my.pattern, reverse = T) {
  new.list <- c()
  for (my.item in my.list) {
    if (reverse) {
      if (!grepl(pattern = my.pattern,
                 x = my.item,
                 fixed = TRUE)) {
        new.list <- c(new.list, my.item)
      }
    }
    else {
      if (grepl(pattern = my.pattern,
                x = my.item,
                fixed = TRUE)) {
        new.list <- c(new.list, my.item)
      }
    }
  }
  return (new.list)
}

my_path <- "results/athal/observed/"
file.list <- list.files(path = my_path)
for (pattern in c("CEU", "IBnr", "Relicts"))
  file.list <- pattern_from_vector(file.list, pattern)
file.list.map <- pattern_from_vector(file.list, "map", reverse = F)


df_list <- list()
for (i in 1:length(file.list.map)) {
  df_file <-
    str_replace(str_remove(file.list.map[i], "map."), "dataset", "sumstats")
  
  map <- read_feather(paste0(my_path, file.list.map[i]))
  df <- read_feather(paste0(my_path, df_file))
  
  df$population <- map$population
  df$index <- 1:nrow(df)
  df$region <- !grepl("all_regions", df_file)
  df$masked <- grepl("masked", df_file)
  
  df_list[[i]] <- df
}

# clean workspace
rm(map, i, pattern, file.list, file.list.map, df_file)

# make long format
df <- do.call(rbind.data.frame, df_list)
df <- df %>%
  pivot_longer(-c(index, population, region, masked), names_to = "stat") %>%
  extract(
    col = "stat",
    into = c("statname", "statid"),
    regex = "(.+)_([[:alnum:]]+)",
    remove = T
  )

df$statid <- as.numeric(df$statid)
df$masked[df$masked] <-
  "masked"
df$masked[df$masked != "masked"] <- "nomask"
df$region[df$region] <-
  "region"
df$region[df$region != "region"] <- "non-peri"
df$masked <- factor(df$masked, levels = c("nomask", "masked"))

message("prepared data")
plot_list <- list()

plot_list[[1]] <- df %>%
  subset(statname == "sfs") %>%
  subset(value > 0) %>%
  ggplot(aes(statid, value, fill = population)) +
  geom_point(
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
    shape = 16,
    col = "gray50",
    alpha = 0.2
  ) +
  stat_summary(
    fun = mean,
    geom = 'point',
    size = 3,
    shape = 23,
    position = position_dodge(width = 0.3)
  ) +
  facet_grid(region ~ masked, scales = "free") +
  scale_y_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  scale_fill_manual(values = wes_palette("Darjeeling1"))+
  labs(x="SFS")


plot_list[[2]] <- df %>%
  subset(statname == "ld") %>%
  subset(value > 0) %>%
  ggplot(aes(statid, value, fill = population)) +
  geom_point(
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
    shape = 16,
    col = "gray50",
    alpha = 0.2
  ) +
  stat_summary(
    fun = mean,
    geom = 'point',
    size = 3,
    shape = 23,
    position = position_dodge(width = 0.5)
  ) +
  facet_grid(region ~ masked, scales = "free") +
  # scale_y_continuous(trans = "log10",
  #                    breaks = logbreak,
  #                    labels = loglabel) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  scale_fill_manual(values = wes_palette("Darjeeling1"))+
  labs(x="LD")



plot_list[[3]] <- df %>%
  subset(statname == "tm_win") %>%
  subset(value > 0) %>%
  ggplot(aes(statid, value, fill = population)) +
  # geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
  #            shape = 16, col="gray50",
  #            alpha=0.2) +
  stat_summary(
    fun = mean,
    geom = 'point',
    size = 3,
    shape = 23,
    position = position_dodge(width = 0.5),
    alpha = 0.6
  ) +
  facet_grid(
    region ~ masked,
    # scales = "free"
    ) +
  scale_y_continuous(trans = "log10",
                     breaks = logbreak,
                     labels = loglabel) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  scale_fill_manual(values = wes_palette("Darjeeling1"))+
  labs(x="TM_WIN")




# explore mean stats of tm_win
mean_df <-
  df %>%
  subset(statname == "tm_win") %>%
  group_by(population, region, masked, statname, statid) %>%
  summarize(value=mean(value, na.rm=F))

mean_df$index <- interaction(
  mean_df$population,
  mean_df$region,
  mean_df$masked,
  mean_df$statname
  )

for (group in unique(mean_df$index)) {
  this_df <- mean_df %>%
    subset(index == group)
  
  this_m <- matrix(data = this_df$value[this_df$statid + 1], nrow = sqrt(nrow(this_df)))
  
  colnames(this_m) <- 1:ncol(this_m)
  rownames(this_m) <- 1:nrow(this_m)
  
  this_df <- melt(this_m)
  colnames(this_df) <- c("x", "y", "p")
  
  plot_list[[length(plot_list) + 1]] <- ggplot(this_df, aes(x = x, y = y, fill = p)) +
    geom_tile()+
    scale_fill_gradient(low = "white", high = "hotpink4") +
    coord_fixed()+
    labs(
      title = group,
      x="n-th segment",
      y="(n+1)-th segment"
      ) +
    theme(
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )
}



message("prepared plots")
pdf(paste0(my_path, "vis_athal_obs.pdf"), useDingbats = F)
for (i in 1:length(plot_list)) {
  show(plot_list[[i]])
}
dev.off()


message("printed plots")
