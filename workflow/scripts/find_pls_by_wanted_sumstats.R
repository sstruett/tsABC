# old file

# this script is a replacement of the more redundant version of the performance
# analysis

# input: sumstat table
# output:
#   transformer file (definition of pls)
#   reduced untransformed sumstat table as txt


library(arrow)

source("workflow/scripts/rfunctions.R")
message("loaded rfunctions")

# read arguments from command line
args = commandArgs(trailingOnly=TRUE)

# explicify output files
outfile_transformer <- args[1]
outfile_reduced_sumstats_abc <- args[2]

# extract max number of pls to be determined
max_num_pls <- args[4]

# extract the sumstat names composition to use
sscomp_id <- as.character(args[5])
sscomp <- strsplit(args[6], "/")[[1]]

# read in data frames from pickle
df_abc = read_feather(args[3])

# remove columns with no variance from abc simulations
df_abc <- remove_columns_with_no_variance(df_abc)

# extract the sumstat composition
sumstats_abc <- obtain_columns_with_wanted_sumstats(df_abc, tolower(sscomp))

# extract the parameters
params_abc <- obtain_columns_with_wanted_sumstats(df_abc, "param")

# print head of params and raw sumstats
cat("____________", "variant parameters", "\n")
print(head(params_abc))
cat("____________", "variant sumstats", "\n")
if (ncol(sumstats_abc) < ncol(params_abc)) {
  print(head(sumstats_abc,))
} else {
  print(head(sumstats_abc[,c(1:ncol(params_abc), ncol(sumstats_abc))]))
}

# define number of pls components to be determined
num_sumstats <- ncol(sumstats_abc)

num_pls <- get_number_of_pls_to_define(num_sumstats, max_num_pls)

# don't do pls if only one column of sumstats
print(num_pls)
stopifnot(as.logical(num_pls))

#open File
dir_prefix <- paste(strsplit(args[1], "/")[[1]][1:(length(strsplit(args[1], "/")[[1]])-1)], collapse = "/")
directory <- paste0(dir_prefix, "/", collapse = "");
dir.create(file.path(directory), showWarnings = F, recursive = T);
message("created: ", directory)
filename<-paste(sscomp_id, collapse = "_");
numComp<-num_pls;

# reassign the stats
stats<-sumstats_abc; params<-params_abc;
message("dimensions of summary statistics data frame: ", paste(dim(stats), collapse = " "))

# use only maximal 100,000 rows to calculate the pls
if (nrow(stats)>100000) {
  stats<-stats[1:100000,]; params<-params[1:100000,]
  stats<-remove_columns_with_no_variance(stats)
}

params <- as.data.frame(params)
stats <- as.data.frame(stats)


#standardize the params


for(i in 1:length(params)) {
  params[,i]<-(params[,i]-mean(params[,i]))/sd(params[,i]);
}
message("done standardize the params")

#force stats in [1,2]
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:length(stats)){
  myMax<-c(myMax, max(stats[,i]));
  myMin<-c(myMin, min(stats[,i]));
  stats[,i]<-1+(stats[,i]-myMin[i])/(myMax[i]-myMin[i]);
}
message("done force stats in [1,2]")

#transform statistics via boxcox  
library("MASS");	
for(i in 1:length(stats)){		
  d<-cbind(stats[,i], params);
  mylm<-lm(as.formula(d), data=d)
  myboxcox<-boxcox(mylm, lambda=seq(-50, 80, 1/10), plotit=F, interp=T, eps=1/50);
  lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);			
  print(paste(names(stats)[i], myboxcox$x[myboxcox$y==max(myboxcox$y)]));
  myGM<-c(myGM, exp(mean(log(stats[,i]))));			
}
message("done transform statistics via boxcox")

#standardize the BC-stats
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:length(stats)){
  stats[,i]<-(stats[,i]^lambda[i] - 1)/(lambda[i]*myGM[i]^(lambda[i]-1));	
  myBCSDs<-c(myBCSDs, sd(stats[,i]));
  myBCMeans<-c(myBCMeans, mean(stats[,i]));		
  stats[,i]<-(stats[,i]-myBCMeans[i])/myBCSDs[i];
}
message("done standardize the BC-stats")

cat("____________", "variant sumstats", "\n")
if (ncol(sumstats_abc) < ncol(params_abc)) {
  print(head(sumstats_abc,))
} else {
  print(head(sumstats_abc[,c(1:ncol(params_abc), ncol(sumstats_abc))]))
}
message("done prepare summary stats for pls")

#perform pls
library("pls");
myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp=numComp);
message("done perform pls")

#write pls to a file
myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
for(i in 2:numComp) { myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i]); }
write.table(cbind(names(stats), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs, myPlsrDataFrame), file=paste(directory, "Routput_", filename, sep=""), col.names=F, row.names=F, sep="\t", quote=F);

#make RMSE plot
pdf(paste(directory, "RMSE_", filename, ".pdf", sep=""));
plot(RMSEP(myPlsr));
dev.off();

message("finished finding pls")

