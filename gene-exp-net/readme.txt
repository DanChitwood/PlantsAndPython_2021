files for gene coexpression network group

### workflow for WCGNA should be: subset.py > obtain subset from gene_FPKM_200501.csv via bash script > filter.py > sft_thresh.R > wgcna.R or wgcna_consensus.R
### can then do export to cytoscape stuff
### alt workflow with subsets provided by other groups: filter.py > sft_thresh.R > wgcna.R or wgcna_consensus.R
### note that making a consensus network requires an additional processing step to keep only genes that survived filter.py for both group a & group b 


subset.ipynb is my method of making a subset from the metadata & should output a functional bash line for retreiving your subset from the original data. I spent way too much time on this but it should be very easy to use & making it helped me get more comfy with python & pandas.

subset.py is a stripped down script from the subset.ipynb notebook with parameter variables declared at the beginning for convenience
please use it (remember to adjust parameters & output file name)

fpkm_filtering.ipynb shows how I filtered data for WGCNA. the ~filtering by low expression across samples~ might be useful for other groups, just don't use the log transform part 

filter.py can be run on a developmental node on smaller subsets, I suggest making a subset then filering but it might run okay on the full set on a processing node. also worth noting that this does a log2(x+1) transform so if you're adapting it for a different purpose make sure to remove that

sft_thresh.R is a script for graphical determination of soft threshold power for use in subsequent WGCNA scripts (wgcna.R & wgcna_consensus.R)

wgcna.R is a script for running single sample WGCNA (tutorial I), needs a single filtered input dataframe (ie output from filter.py)

wgcna_consensus.R is a script for running consensus analysis in WGCNA (tutorial II), needs two filtered input dataframes; will output correlation networks between each subset network & a consensus network
//I think this will be useful for determining modules of interest in each subsample 
//logic is: compare each network to the consensus network & look at what's not included in the consensus
//this might be able to be played with a bit to make a dis-consensus network

metafix.ipynb is just documentation of how I generated the indexed metadata for subsetting. you should be fine to use the indexed_meta.csv file located at /mnt/research/HRT841_F21/Class_project/subsets as input for the subset.ipynb workflow



- SMH
