library(abc)
library(tidyr)
library(pbapply)

save.image(file = "rdev.RData")
stop("saved rdev.RData")
setwd(
  "/Users/struett/MPIPZ/netscratch-2/dep_tsiantis/grp_laurent/struett/git_tsABC/tsABC/"
)
load("rdev.RData")

# perform a model choice between two proposed models

# read data
sumstats <- read.csv(snakemake@input$sumstats, sep = "\t")
sumstats <-
  sumstats[, grep("LinearCombination", colnames(sumstats))]
alternative_sumstats <-
  read.csv(snakemake@input$alternative_sumstats, sep = "\t")
alternative_sumstats <-
  alternative_sumstats[, grep("LinearCombination", colnames(alternative_sumstats))]
podstats <- read.csv(snakemake@input$podstats, sep = "\t")
podparam <- podstats[, grep("param", colnames(podstats))]
podstats <-
  podstats[, grep("LinearCombination", colnames(podstats))]
podindex <- interaction(podparam, sep = "_")


# make podindex being same as in config file
podconfig <- matrix(data=as.numeric(snakemake@config$ABC$performance$pods), nrow = 4)
podconfig

# read paremeter for inference
ntolerated <- as.numeric(snakemake@wildcards$tolid)
npls <- as.numeric(snakemake@wildcards$plsid)


# prepare data
sumstats$model = "transition"
alternative_sumstats$model = "constant_selfing"
sumstats_model_choice <-
  rbind.data.frame(sumstats, alternative_sumstats)
model_indices <- sumstats_model_choice$model
sumstats_model_choice$model <- NULL
sumstats_model_choice <-
  sumstats_model_choice[, which(colnames(sumstats_model_choice) %in% paste("LinearCombination", 1:npls -
                                                                             1, sep = "_"))]
podstats <-
  podstats[, which(colnames(podstats) %in% paste("LinearCombination", 1:npls - 1, sep = "_"))]
                                                                        

model_choice_result <- pbapply(podstats, 1, function(x) {
  a <- postpr(
    x,
    model_indices,
    sumstats_model_choice,
    tol = ntolerated / nrow(sumstats),
    method = "mnlogistic",
    corr = TRUE   # corr seems not to work
  )
  return(summary(a, rejection=T, print=F))
})

bayes_rejection <- numeric(length = length(model_choice_result))
bayes_mnlogistic <- numeric(length = length(model_choice_result))
for (i in 1:length(model_choice_result)) {
  res <- model_choice_result[[i]]
  bayes_rejection[i] <- model_choice_result[[i]]$rejection$BayesF["transition", "constant_selfing"]
  bayes_mnlogistic[i] <- model_choice_result[[i]]$mnlogistic$BayesF["transition", "constant_selfing"]
}

# save Bayes Factors RDS



# plot
pdf(snakemake@output$bayes_plot)
plot(podindex, bayes_rejection, log="y")
plot(podindex, bayes_mnlogistic, log="y")
dev.off()
