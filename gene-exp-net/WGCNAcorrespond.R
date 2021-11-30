### WGCNA consensus modules & how they correspond to individual modules
library(dplyr)
library(WGCNA)
library(ggplot2)
library(tidyr)
library(data.table)

# setup
setwd('/home/beebo_bebop/HRT841 - Plants & Python')
options(stringsAsFactors = F) #WGCNA req
root_data=read.csv('root_consensus.csv')
leaf_data=read.csv('leaf_consensus.csv')
#enableWGCNAThreads() #need this in terminal R but using studio so disabled

#convert root data from long to wide
root_gene_names <- root_data[,1] # get gene names
root_data_t <- as.data.frame(t(root_data[,-c(1)])) # remove gene name column and transpose dataframe
names(root_data_t) <- root_gene_names # add gene names
#convert leaf data from long to wide
leaf_gene_names <- leaf_data[,1] # get gene names
leaf_data_t <- as.data.frame(t(leaf_data[,-c(1)])) # remove gene name column and transpose dataframe
names(leaf_data_t) <- leaf_gene_names # add gene names


### prep by doing individual networks
# build root network and identify modules automatically, soft thresh = 8
root_net <- blockwiseModules(root_data_t,power=8,
                             TOMType = "unsigned", minModuleSize = 30,
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = TRUE, maxBlockSize = 16000,
                             saveTOMFileBase = "rootTOM",
                             verbose = 3)
# plot dendrogram with colors of assigned modules
root_mergedColors <- labels2colors(root_net$colors) # convert module labels to colors
plotDendroAndColors(root_net$dendrograms[[1]],root_mergedColors[root_net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang=0.03, 
                    addGuide = TRUE, guideHang = 0.05)
table(root_net$colors)
#save the module assignment and module eigengene information necessary for subsequent analysis
rootLabels = root_net$colors
rootColors = labels2colors(root_net$colors)
rootMEs = root_net$MEs
rootTree = root_net$dendrograms[[1]]

# build leaf network and identify modules automatically, soft thresh = 8
leaf_net <- blockwiseModules(leaf_data_t,power=8,
                             TOMType = "unsigned", minModuleSize = 30,
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = TRUE, maxBlockSize = 16000,
                             saveTOMFileBase = "leafTOM",
                             verbose = 3)

# plot dendrogram with colors of assigned modules
leaf_mergedColors <- labels2colors(leaf_net$colors) # convert module labels to colors
plotDendroAndColors(leaf_net$dendrograms[[1]],leaf_mergedColors[leaf_net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang=0.03, 
                    addGuide = TRUE, guideHang = 0.05)
table(leaf_net$colors)
#save the module assignment and module eigengene information necessary for subsequent analysis
leafLabels = leaf_net$colors
leafColors = labels2colors(leaf_net$colors)
leafMEs = leaf_net$MEs
leafTree = leaf_net$dendrograms[[1]]

### step 1 of tutorial
nSets = 2 #two data sets
labels = c("root", "leaf")
multiExpr = vector(mode = "list", length = nSets) # sets up vector for both datasets to be put in
multiExpr[[1]] = list(data = root_data_t) # puts data into multiExpr as a list
names(multiExpr[[1]]$data) = root_data$Sample # provides geneID as names
rownames(multiExpr[[1]]$data) = colnames(root_data[,-1]) #provides rownames from colnames in original data (minus first col bc that's sample)
multiExpr[[2]] = list(data = leaf_data_t)
names(multiExpr[[2]]$data) = leaf_data$Sample
rownames(multiExpr[[2]]$data) = colnames(leaf_data[,-1])
exprSize=checkSets(multiExpr) # verify data format
nGenes=exprSize$nGenes # initialize variables for later
nSamples=exprSize$nSamples # initialize variables for later
#they say to output various things to a file but unnecessary since using single session

###step 2 of tutorial
powers = c(c(1:10), seq(from = 12, to=20, by=2)) # soft threshold powers
powerTables = vector(mode = 'list', length = nSets) #list to hold scale-free analysis results
for(set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage() ## for loop to call topology analysis for each data set

# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points plot.window
ylim = matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
#pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col])
    addGrid()
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set])
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set])
  if (col==1)
  {
    legend("bottomright", legend = labels, col = colors, pch = 20)
  } else
    legend("topright", legend = labels, col = colors, pch = 20)
}

#build network.. choosing power 8 but idk what's best.. maybe 7 maybe 12
#blockwise is necessary on laptop.. single step uses over 11Gb memory
net = blockwiseConsensusModules(
  multiExpr, power = 8, minModuleSize = 30, deepSplit = 2,
  pamRespectsDendro = FALSE, maxBlockSize = 16000,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)
consMEs = net$multiMEs; #guessing these are the conserved modules
moduleLabels = net$colors
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]];
#dendrogram
sizeGrWindow(8,6);
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")

### step 4 of tutorial - relating consensus to set-specific modules
#Isolate the module labels in the order they appear in ordered module eigengenes
rootModuleLabels = substring(names(rootMEs), 3)
leafModuleLabels = substring(names(leafMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
rootModules = labels2colors(as.numeric(rootModuleLabels))
leafModules = labels2colors(as.numeric(leafModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of female and consensus modules
nRootMods = length(rootModules)
nLeafMods = length(leafModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTableRoot = matrix(0, nrow = nRootMods, ncol = nConsMods)
pTableLeaf = matrix(0, nrow = nLeafMods, ncol = nConsMods)
rootCountTbl = matrix(0, nrow = nRootMods, ncol = nConsMods)
leafCountTbl = matrix(0, nrow = nLeafMods, ncol = nConsMods)
# Execute all pairwaise comparisons
for (fmod in 1:nRootMods)
  for (cmod in 1:nConsMods)
  {
    rootMembers = (rootColors == rootModules[fmod]);
    consMembers = (moduleColors == consModules[cmod]);
    pTableRoot[fmod, cmod] = -log10(fisher.test(rootMembers, consMembers, alternative = "greater")$p.value);
    rootCountTbl[fmod, cmod] = sum(rootColors == rootModules[fmod] & moduleColors ==
                                     consModules[cmod])
  }

# Truncate p values smaller than 10^{-50} to 10^{-50}
pTableRoot[is.infinite(pTableRoot)] = 1.3*max(pTableRoot[is.finite(pTableRoot)])
pTableRoot[pTableRoot>50 ] = 50 
# Marginal counts (really module sizes)
rootModTotals = apply(rootCountTbl, 1, sum)
consModTotals = apply(rootCountTbl, 2, sum)
# Actual plotting
sizeGrWindow(10,7 )
#pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7)
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTableRoot,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", rootModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("root ", rootModules, ": ", rootModTotals, sep=""),
               textMatrix = rootCountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of root set-specific and root-leaf consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)

# again with leaves..
for (fmod in 1:nLeafMods)
  for (cmod in 1:nConsMods)
  {
    leafMembers = (leafColors == leafModules[fmod]);
    consMembers = (moduleColors == consModules[cmod]);
    pTableLeaf[fmod, cmod] = -log10(fisher.test(leafMembers, consMembers, alternative = "greater")$p.value);
    leafCountTbl[fmod, cmod] = sum(leafColors == leafModules[fmod] & moduleColors ==
                                     consModules[cmod])
  }
# Truncate p values smaller than 10^{-50} to 10^{-50}
pTableLeaf[is.infinite(pTableLeaf)] = 1.3*max(pTableLeaf[is.finite(pTableLeaf)])
pTableLeaf[pTableLeaf>50 ] = 50 
# Marginal counts (really module sizes)
leafModTotals = apply(leafCountTbl, 1, sum)
consModTotals = apply(leafCountTbl, 2, sum)
# Actual plotting
sizeGrWindow(10,7 )
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTableLeaf,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", leafModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("leaf ", leafModules, ": ", leafModTotals, sep=""),
               textMatrix = leafCountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of leaf set-specific and root-leaf consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
