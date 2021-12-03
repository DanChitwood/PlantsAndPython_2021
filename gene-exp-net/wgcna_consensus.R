
### WGCNA 
### individual networks for each subset & consensus networks & how they correspond to individual modules
### should output 'dendro_a.pdf', 'dendro_b.pdf', 'dendro_cons.pdf', 'merging.pdf', &  'consensus_vs_subset.pdf'
### the first three are dendrograms for each network (sample a, sample b, consensus)
### 'merging.pdf' shows the module merging for each network
### 'consensus_vs_subset.pdf' shows heatmap of module discongruity between each subset & consensus (# of genes from each module not found in consensus)

### note: to run WGCNA on hpcc you need to use R/4.1 && install BiocManager & WGCNA to local lib directories 

library(WGCNA)
library(tidyr)


### setup
sft_thresh_a= 8   #initialize soft thresholding power variables
sft_thresh_b= 8
sft_thresh_cons= 8
minModuleSize= 30   # initialize minimum module size variable (>x genes, 30 recommended by WGCNA tutorial)
#colors= c("#0f0f12", "#505359", "#b6bfbc", "#f2fbff", “#5ee7ff", "#00a1db", "#1d5bb8", "#1f2c66", "#1b5245", "#2e8f46", "#58d92e", "#cbff70", "#ffff8f", "#ffdf2b", "#f0771a", "#e32239", "#851540", "#401a24", "#9c3b30", "#3929cc", "#8a5cff", "#ffbca6", "#eb75be", "#77388c”)   # 24x custom color init, trying to make it more readable
options(stringsAsFactors=  F)   # WGCNA req
#setwd('/home/beebo_bebop/HRT841 - Plants & Python')   # disabled bc working from my home directory on hpcc due to limited write size issue
data_a_long= read.csv('root_consensus.csv')   # read in Subset A
data_b_long= read.csv('leaf_consensus.csv')   # &Subset B
enableWGCNAThreads()   # need this in Rscript but supposedly not functional in Rstudio


### convert data from long to wide
# may not be necessary if using subsets provided by other groups, but mine are made as a raw column subset of the original data rather than the computationally heavy method used by the data cleaning group

# for Subset A
a_gene_names=   data_a_long[,1]   # get gene names
data_a=   as.data.frame(t(data_a_long[,-c(1)]))   # remove gene name column and transpose dataframe
names(data_a)=   a_gene_names   # add gene names

# again for Subset B
b_gene_names=   data_b_long[,1]
data_b=   as.data.frame(t(data_b_long[,-c(1)]))
names(data_b)=   b_gene_names


### prep by doing individual networks
# should do soft thresh determination with sft_thresh.R for each subset before using this script
# here i chose to split the difference between best for a & best for b. not sure what best practice is..
# maybe pick different threshold value for each?

# build network from Subset A & identify modules 
# soft thresh & minimum module size are initialized in setup section
adj_a= adjacency(data_a, power= sft_thresh_a)   # defines adjacency
tom_a= TOMsimilarity(adj_a)   # turns adjacency into topological overlap matrix
tom_a_diss= 1-tom_a   # dissimilarity from similarity
dend_a= hclust(as.dist(tom_a_diss),method= "average")   # make dendrogram
modlabels_a= cutreeDynamic(dendro=  dend_a, distM=  tom_a_diss, deepSplit=  1, pamRespectsDendro=  FALSE, minClusterSize=  minModuleSize)   # determine modules & list module label for each gene
modcolors_a= labels2colors(modlabels_a)   # makes the colors for the dendrogram #convert label list to color list
cutheight_a= dynamicMergeCut(length(names(data_a_long))-1,0.85,2.35)   # function to select cut height using sample size, merge correlation & z statistic // tutorial uses set values chosen based on visual inspection of the dendrograms but this function is available & seems to have good output but my stats knowledge is a bit lacking for determining whether it’s used appropriately
merge_a= mergeCloseModules(data_a, modlabels_a, cutHeight= cutheight_a, verbose= 3)   # function to merge modules with similar expression profiles
modlabels_merged_a= merge_a$colors   # module labels for each gene after merge
modcolors_merged_a=labels2colors(modlabels_merged_a)   # convert to list of colors
eigens_merged_a= merge_a$newMEs   # list of new module eigengenes

# re-order modules & initialize some variables for later
a_eigens=  orderMEs(eigens_merged_a, greyName=  "ME0") #reorders modules by eigengene correlation // modules as columns, eigengenes as rows
a_labs=  substring(names(a_eigens),3) #list nimeric labels for each *module*, in order passed by previous line
a_mods=  labels2colors(as.numeric(a_labs)) #list of strings containing color value for each *gene*
a_colors= modcolors_merged_a


# output dendrogram with modules for Subset A
pdf('dendro_a.pdf', wi=9, he=5)
sizeGrWindow(12,9)
plotDendroAndColors(dend_a, modcolors_merged_a, xlab= "Dynamic Tree Cut", dendroLabels= FALSE, hang= 0.03, addGuide= TRUE, guideHang= 0.05, main=  "Gene clustering on TOM-based dissimilarity with modules for Subset A")
dev.off()


### build network from Subset B & identify modules 
# soft thresh & minimum module size are initialized in setup section
adj_b= adjacency(data_b, power= sft_thresh_b)   # defines adjacency
tom_b= TOMsimilarity(adj_b)   # turns adjacency into topological overlap matrix
tom_b_diss= 1-tom_b   # dissimilarity from similarity
dend_b= hclust(as.dist(tom_b_diss),method= "average")   # make dendrogram
modlabels_b= cutreeDynamic(dendro=  dend_b, distM=  tom_b_diss, deepSplit=  1, pamRespectsDendro=  FALSE, minClusterSize=  minModuleSize)   # determine modules & list module label for each gene
modcolors_b= labels2colors(modlabels_b)   # makes the colors for the dendrogram #convert label list to color list
cutheight_b= dynamicMergeCut(length(names(data_b_long))-1,0.85,2.35)   # function to select cut height using merge correlation & z statistic // tutorial uses set values chosen based on visual inspection of the dendrograms but this function is available & seems to have good output but my stats knowledge is a bit lacking for determining whether it’s used appropriately
merge_b= mergeCloseModules(data_b, modlabels_b, cutHeight= cutheight_b, verbose= 3)   # function to merge modules with similar expression profiles
modlabels_merged_b= merge_b$colors   # module labels for each gene after merge
modcolors_merged_b=labels2colors(modlabels_merged_b)   # convert to list of colors
eigens_merged_b= merge_b$newMEs   # list of new module eigengenes

# re-order modules & initialize some variables for later
b_eigens=  orderMEs(eigens_merged_b, greyName=  "ME0") #reorders modules by eigengene correlation // modules as columns, eigengenes as rows
b_labs=  substring(names(b_eigens),3) #list nimeric labels for each *module*, in order passed by previous line
b_mods=  labels2colors(as.numeric(b_labs)) #list of strings containing color value for each *gene*
b_colors= modcolors_merged_b

# output dendrogram with modules for Subset B
pdf('dendro_b.pdf', wi=9, he=6)
sizeGrWindow(12,9)
plotDendroAndColors(dend_b, modcolors_merged_b, xlab= "Dynamic Tree Cut", dendroLabels= FALSE, hang= 0.03, addGuide= TRUE, guideHang= 0.05, main=  "Gene clustering on TOM-based dissimilarity with modules for Subset B")
dev.off()


### step 1 of tutorial II // initialize vectors & variables for loops that put the two TOM vectors together
# multiExpr needs to be set up this way for the module merging function
nSets= 2   #two data sets
labels= c("Subset A", "Subset B")
multiExpr=  vector(mode=  "list", length=  nSets)   #sets up vector for both datasets to be placed in
multiExpr[[1]]=  list(data=  data_a)   #puts data into multiExpr as a list
names(multiExpr[[1]]$data)=  data_a_long$Sample   #provides geneID as names as column names
rownames(multiExpr[[1]]$data)=  colnames(data_a_long[,-1])   #provides rownames from colnames in original data (minus first col bc that's sample)
multiExpr[[2]]=  list(data=  data_b)
names(multiExpr[[2]]$data)=  data_b_long$Sample
rownames(multiExpr[[2]]$data)=  colnames(data_b_long[,-1])
exprSize= checkSets(multiExpr)   #verify data format
nGenes= exprSize$nGenes   #initialize variables for later
nSamples= exprSize$nSamples   #initialize variables for later


### step 2 of tutorial II // put TOMs together

### build consensus network
tom_cons= pmin(tom_a,tom_b)   #calculate consensus TOM
tom_cons_diss= 1-tom_cons   #dissimilarity from similarity
dend_cons= hclust(as.dist(tom_cons_diss),method= "average")   #make dendrogram
modlabels_cons= cutreeDynamic(dendro=dend_cons, distM=tom_cons_diss, deepSplit=  1, pamRespectsDendro=  FALSE, minClusterSize=  minModuleSize)    #construct modules & list module label for each gene
modcolors_cons= labels2colors(modlabels_cons)
shorter=ifelse(length(names(data_a_long))>length(names(data_b_long)),length(names(data_b_long))-1,length(names(data_a_long))-1)   #use smaller of subset sample size for cutheight
cutheight_cons= dynamicMergeCut(shorter,0.8,2.35)
merge_cons= mergeCloseModules(multiExpr,modlabels_cons,cutHeight= cutheight_cons,verbose= 3)
modlabels_merged_cons= merge_cons$colors
modcolors_merged_cons= labels2colors(modlabels_merged_cons)
eigens_merged_cons= merge_cons$newMEs

# re-order modules & initialize some variables for later
cons_eigens=  consensusOrderMEs(eigens_merged_cons, greyName=  "ME0") #reorders modules by eigengene correlation // modules as columns, eigengenes as rows
cons_labs=  substring(names(cons_eigens[[1]]$data),3) #list nimeric labels for each *module*, in order passed by previous line
cons_mods=  labels2colors(as.numeric(cons_labs)) #list of strings containing color value for each *gene*
cons_colors=modcolors_merged_cons

pdf('dendro_cons.pdf', wi=9, he=6)
sizeGrWindow(12,9)
plotDendroAndColors(dend_cons, modcolors_merged_cons, xlab= "Dynamic Tree Cut", dendroLabels= FALSE, hang= 0.03, addGuide= TRUE, guideHang= 0.05, main=  "Gene clustering on TOM-based dissimilarity with modules for Consensus Module")
dev.off()


### if we want to show effect of filtering by eigengene similarity
eigenlist_a= moduleEigengenes(data_a, colors= modcolors_a) #creates list of eigengenes (imaginary representative of module)
eigen_a= eigenlist_a$eigengenes
eigendiss_a= 1-cor(eigen_a) #eigengene dissimilarity as function of correlation of expression
eigentree_a= hclust(as.dist(eigendiss_a), method=  "average")
eigenlist_b= moduleEigengenes(data_b, colors= modcolors_b)
eigen_b= eigenlist_b$eigengenes
eigendiss_b= 1-cor(eigen_b)
eigentree_b= hclust(as.dist(eigendiss_b), method=  "average")
eigen_cons= multiSetMEs(multiExpr, colors=  NULL, universalColors=  modcolors_cons)
eigendiss_cons= consensusMEDissimilarity(eigen_cons)
eigentree_cons= hclust(as.dist(eigendiss_cons), method=  "average")

pdf('merging.pdf')   #output pdf showing module merging & resulting dendrogram/module graph
sizeGrWindow(12,9)
plot(eigentree_a, main= "Clustering of Module Eigengenes for Subset A", xlab=  "", sub=  "")
abline(h= cutheight_a, col= "red")
plotDendroAndColors(dend_a, cbind(modcolors_a,modcolors_merged_a), c("Dynamic Tree Cut","Merged Dynamic"), dendroLabels= FALSE, hang= 0.03, addGuide= TRUE, guideHang= 0.05, main=  "Gene clustering on TOM-based dissimilarity with modules for Subset A")
plot(eigentree_b, main= "Clustering of Module Eigengenes for Subset B", xlab=  "", sub=  "")
abline(h= cutheight_b, col= "red")
plotDendroAndColors(dend_b, cbind(modcolors_b,modcolors_merged_b), c("Dynamic Tree Cut","Merged Dynamic"), dendroLabels= FALSE, hang= 0.03, addGuide= TRUE, guideHang= 0.05, main=  "Gene clustering on TOM-based dissimilarity with modules for Subset B")
plot(eigentree_cons, main= "Clustering of Module Eigengenes for Consensus Model", xlab=  "", sub=  "")
abline(h= cutheight_cons, col= "red")
plotDendroAndColors(dend_cons, cbind(modcolors_cons,modcolors_merged_cons), c("Dynamic Tree Cut","Merged Dynamic"), dendroLabels= FALSE, hang= 0.03, addGuide= TRUE, guideHang= 0.05, main=  "Gene clustering on TOM-based dissimilarity with modules for Consensus Model")
dev.off()


### step 4 of tutorial - relating consensus to set-specific modules

# count modules for each
a_nMods=  length(a_mods)
b_nMods=  length(b_mods)
cons_nMods=  length(cons_mods)

# Initialize tables of p-values and of the corresponding counts
a_pTable=  matrix(0, nrow=  a_nMods, ncol=  cons_nMods)
b_pTable=  matrix(0, nrow=  b_nMods, ncol=  cons_nMods)
a_CountTbl=  matrix(0, nrow=  a_nMods, ncol=  cons_nMods)
b_CountTbl=  matrix(0, nrow=  b_nMods, ncol=  cons_nMods)

# Execute all pairwaise comparisons for sample a
for (amod in 1:a_nMods)
  for (cmod in 1:cons_nMods)   #loops through consensus modules for each module in subset a
  {
    a_mbrs=  (a_colors == a_mods[amod])   #creates list of bool values for each gene (true if in [amod] module)
    cons_mbrs=  (cons_colors == cons_mods[cmod])   #creates list of bool values for each gene (true if in [cmod] module)
    a_pTable[amod, cmod]=  -log10(fisher.test(a_mbrs, cons_mbrs, alternative=  "greater")$p.value) #determines p value of by fisher test and writes -log10(p) to vector // p=  liklihood that odds ratio is >1 for ~x gene in [amod] is present in [cmod]~
    a_CountTbl[amod, cmod]=  sum(a_colors == a_mods[amod] & cons_colors == cons_mods[cmod]) #determines number of genes belonging to both [amod] & [cmod], returns value to vector
  }

# Truncate p values smaller than 10^{-50} to 10^{-50} // not really sure why this is necessary but tutorial says to
a_pTable[is.infinite(a_pTable)] = 1.3*max(a_pTable[is.finite(a_pTable)])
a_pTable[a_pTable>50 ] = 50 

# Marginal counts (module sizes?)
a_ModTotals=  apply(a_CountTbl, 1, sum) #creates list with total number of genes in module for each module
consa_ModTotals=  apply(a_CountTbl, 2, sum)

# pairwise again for sample b
for (bmod in 1:b_nMods)
  for (cmod in 1:cons_nMods)
  {
    b_mbrs=  (b_colors == b_mods[bmod])
    cons_mbrs=  (cons_colors == cons_mods[cmod])
    b_pTable[bmod, cmod]=  -log10(fisher.test(b_mbrs, cons_mbrs, alternative=  "greater")$p.value)
    b_CountTbl[bmod, cmod]=  sum(b_colors==   b_mods[bmod] & cons_colors==   cons_mods[cmod])
  }
  
# Truncate p values smaller than 10^{-50} to 10^{-50}
b_pTable[is.infinite(b_pTable)]=  1.3*max(b_pTable[is.finite(b_pTable)])
b_pTable[b_pTable>50 ]=  50 

# Marginal counts (really module sizes)
b_ModTotals=  apply(b_CountTbl, 1, sum)
consb_ModTotals=  apply(b_CountTbl, 2, sum)

### consensus plotting, set to output a single pdf with two pages
# refer to the caption of tutorial version of this, but basically is a heatmap of module overlap where red indicates 'high' overlap between subset & consensus
pdf('consensus_vs_subset.pdf', wi=10, he=7)
par(mfrow= c(1,1))
par(cex=  1.0)
par(mar= c(8, 10.4, 2.7, 1)+0.3)
# Use function labeledHeatmap to produce the color-coded table
labeledHeatmap(Matrix=  a_pTable,
               xLabels=  paste(" ", cons_mods),
               yLabels=  paste(" ", a_mods),
               colorLabels=  TRUE,
               xSymbols=  paste("Consensus ", cons_mods, ": ", consa_ModTotals, sep= ""),
               ySymbols=  paste("Subest A ", a_mods, ": ", a_ModTotals, sep= ""),
               #textMatrix=  a_CountTbl,
               colors=  greenWhiteRed(100)[50:100],
               main=  "Correspondence of Subset A modules vs Subset A/B Consensus modules",
               cex.text=  1.0, cex.lab=  1.0, setStdMargins=  FALSE)
labeledHeatmap(Matrix=  b_pTable,
               xLabels=  paste(" ", cons_mods),
               yLabels=  paste(" ", b_mods),
               colorLabels=  TRUE,
               xSymbols=  paste("Consensus ", cons_mods, ": ", consb_ModTotals, sep= ""),
               ySymbols=  paste("Subset B ", b_mods, ": ", b_ModTotals, sep= ""),
               #textMatrix=  b_CountTbl,
               colors=  greenWhiteRed(100)[50:100],
               main=  "Correspondence of Subset B modules vs Subset A/B Consensus modules",
               cex.text=  1.0, cex.lab=  1.0, setStdMargins=  FALSE)
dev.off()


#### very rough export to cytoscape

# Read in the annotation file // we don't have one currently but maybe later? leaving in but disabled
#annot=  read.csv(file=  "GeneAnnotation.csv")

# Select modules // okay but can i select all the modules?
#a_modules=  c("brown", "red")

# Select module probes
#probes=  names(data_a)
#a_inModule=  is.finite(match(a_colors, a_modules))
#a_modProbes=  probes[a_inModule]
#modGenes=  annot$gene_symbol[match(a_modProbes, annot$some_column_idkwhy)] #enable if using annotation, not sure what the second argument of match() should be

# Select the corresponding Topological Overlap
#a_modTOM=  TOM[a_inModule, a_inModule]
#dimnames(a_modTOM)=  list(a_modProbes, a_modProbes)

# Export the network into edge and node list files Cytoscape can read
#cyt=  exportNetworkToCytoscape(modTOM,
#edgeFile=  paste("CytoscapeInput-edges-", paste(modules, collapse= "-"), ".txt", sep= ""),
#nodeFile=  paste("CytoscapeInput-nodes-", paste(modules, collapse= "-"), ".txt", sep= ""),
#weighted=  TRUE,
#threshold=  0.02,
#nodeNames=  modProbes,
#altNodeNames=  modGenes,
#nodeAttr=  moduleColors[inModule])
