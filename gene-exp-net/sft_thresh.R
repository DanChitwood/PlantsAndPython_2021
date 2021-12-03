### script for determining soft thresholding power for WGCNA network building
### tutorial says to pick first number to cross threshold of 0.9 for scale independence
### some other resources say to pick number that represents inflection point since the graph isn't always clean..
### mean connectivity for power chosen should also be below 100
### this will output a 2 page pdf, first page scale independence & second page mean connectivity
### you will need to either change the output_file variable for different data or rename the file sft_thresh.pdf after you make it so it isn't overwritten
### use xdg-open sft_thresh.pdf command to open graph from terminal

library(WGCNA)
library(tidyr)

###setup
output_file='sft_thresh.pdf'
#setwd('/home/beebo_bebop/HRT841 - Plants & Python') #disables bc working from my home directory on hpcc
options(stringsAsFactors = F) #WGCNA req
data=read.csv('root_consensus.csv') #one type of tissue
enableWGCNAThreads() #need this in cli R but supposedly not studio, enables multi-threading

###convert data from long to wide
#may not be necessary if using subsets provided by other groups, but mine are made as a raw column subset of the original data rather than the computationally heavy method used by the data cleaning group
gene_names <- data[,1] # get gene names
data_t <- as.data.frame(t(data[,-c(1)])) # remove gene name column and transpose dataframe
names(data_t) <- gene_names # add gene names

###make soft threshold graph (per WGCNA tutorial I step 2)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(data_t, powerVector = powers, verbose = 5) #network topology analysis

pdf(output_file)
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main= paste("Scale independence"));
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
# line to choose best value   manually
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()
