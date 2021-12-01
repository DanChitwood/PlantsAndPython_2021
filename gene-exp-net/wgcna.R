### single sample wgcna

library(WGCNA)
library(tidyr)

# setup
file_name = "root_consensus" # with no extension
sft_thresh=8 #initialize soft thresholding power variable
minModuleSize=30 #initialize minimum module size variable (>x genes, 30 recommended by WGCNA tutorial)
eigendiss_thresh=0.25 #merge modules where expression correlation of module eigengene is >1-x (WGCNA tutorial I part 2.b.5)
#setwd('/home/beebo_bebop/HRT841 - Plants & Python') #disabled bc working from my home directory on hpcc
options(stringsAsFactors = F) #WGCNA req
data=read.csv(paste(file_name,'.csv',sep='')) #one type of tissue

enableWGCNAThreads() #need this in cli R but supposedly not studio

###convert data from long to wide
#may not be necessary if using subsets provided by other groups, but mine are made as a raw column subset of the original data rather than the computationally heavy method used by the data cleaning group
#for first subset
gene_names <- data[,1] # get gene names
data_t <- as.data.frame(t(data[,-c(1)])) # remove gene name column and transpose dataframe
names(data_t) <- gene_names # add gene name

### make network
# should do soft thresh determination with sft_thresh.R before using this script

# build network from first subset and identify modules 
# soft thresh & minimum module size are initialized in setup section
adj=adjacency(data_t, power=sft_thresh) #defines adjacency
tom=TOMsimilarity(adj) #turns adjacency into topological overlap matrix
diss_tom=1-tom #dissimilarity from similarity
dend=hclust(as.dist(diss_tom),method="average") #make dendrogram
modules=cutreeDynamic(dendro = dend, distM = diss_tom, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize) #determine modules
modcolors=labels2colors(modules) #makes the colors for the dendrogram
merge=mergeCloseModules(data_t,modcolors,cutHeight=eigendiss_thresh,verbose=3)
merge_modcolors=merge$colors
merge_eigen=merge$newMEs

# prepare data for heatmap (just plot a few random genes)
nSelect = 600
set.seed(20)
select = sample(ncol(data_t),size=nSelect)
select_tom = diss_tom[select,select]
select_tree = hclust(as.dist(select_tom),method='average')
selectcolors = modcolors[select];
plot_diss = select_tom^10
diag(plot_diss) = NA

# plot dendro and heatmap
pdf(paste('dendro_and_heatmap-',file_name,'.pdf',sep=''))
plotDendroAndColors(dend, merge_modcolors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, main = "Gene clustering on TOM-based dissimilarity with modules")
TOMplot(plot_diss,select_tree,selectcolors,main=paste('Network heatmap plot (',nSelect,' random selected genes)',sep = ''))
dev.off()
