# Gene coexpression network group

## Workflow

- workflow for WCGNA should be: `subset.py` > obtain subset from `gene_FPKM_200501.csv` via bash script > `filter.py` > `sft_thresh.R` > `wgcna.R` or `wgcna_consensus.R`
- can then do export to cytoscape stuff
- alt workflow with subsets provided by other groups: `filter.py` > `sft_thresh.R` > `wgcna.R` or `wgcna_consensus.R`
- note that making a consensus network requires an additional processing step to keep only genes that survived `filter.py` for both group a & group b 

## Files

- `subset.ipynb` is my method of making a subset from the metadata & should output a functional bash line for retreiving your subset from the original data. I spent way too much time on this but it should be very easy to use & making it helped me get more comfy with python & pandas.

- `subset.py` is a stripped down script from the subset.ipynb notebook with parameter variables declared at the beginning for convenience
please use it (remember to adjust parameters & output file name)
	- I would like to tweak this some to take CLI inputs for filter initialization
	- ideally this would make it possible to filter by any metadata column (or all, or just some), & input filtering parameters / desired output without editing the script
	- also considering making an option to load in the full data & output the subset.csv directly rather than just providing a bash sommand
	- as is, it works well to run locally with low resources, & for small subsets the bash cmd is pretty quick at retrieving subsets
	- I think this is ideal usage for small samples bc it minimizes computational resources needed
	- BUT, for large subsets it would likely be more efficient if this were to load the data & output the subsets itself

- `fpkm_filtering.ipynb` shows how I filtered data for WGCNA. the *filtering by low expression across samples* might be useful for other groups, just don't use the log transform part 

- `filter.py` can be run on a developmental node on smaller subsets, I suggest making a subset then filering but it might run okay on the full set on a processing node. also worth noting that this does a log2(x+1) transform so if you're adapting it for a different purpose make sure to remove that

- `sft_thresh.R` is a script for graphical determination of soft threshold power for use in subsequent WGCNA scripts (`wgcna.R` & `wgcna_consensus.R`)
	- the 'deluxe' version comes from tutorial II for consensus modules & provides more metrics but doesn't currently function

- `wgcna.R` is a script for running single sample WGCNA (tutorial I), needs a single filtered input dataframe (ie output from `filter.py`)

- `same_genes.py` is a script needed to prepare subsets for WGCNA consensus analysis after filtering since different subsets will have some variation in number of genes after filtering
	- might roll this into filtering if I get around to making the python scripts interactive

- `wgcna_consensus.R` is a script for running consensus analysis in WGCNA (tutorial II), needs two filtered input dataframes; will output correlation networks between each subset network & a consensus network
	- I think this will be useful for determining modules of interest in each subsample 
	- logic is: compare each network to the consensus network & look at what's not included in the consensus
	- this might be able to be played with a bit to make a dis-consensus network
	- this takes >10Gb free RAM, probably more, even with small subsets & will likely require much more with larger subsets, but runs relatively quickly

- `WGCNAcorrespond.R` is the RStudio version of wgcna_consensus.R, it's kinda messy but works to run small subsets locally
	- this uses a different command to create dendrograms which runs better on low resource workstations but has some limitations

- `subsets.R` is the initial RStudio version of subset.py
	- kinda sloppy & not necessary to keep, it's just how I initially made the subset code bc i have a background in R
	- converting was super helpful in getting more comfotrable with python though
	- seriously don't use this, the python version is much better

- `metafix.ipynb` is just documentation of how I generated the indexed metadata for subsetting. 
	- you should be fine to use the `indexed_meta.csv` file located at `/mnt/research/HRT841_F21/Class_project/subsets` as input for the `subset.ipynb` workflow


## Running R scripts on hpcc

Module setups:

```
module purge
module load GCC/8.3.0 OpenMPI/3.1.4
module load R/4.1
```

open R 

```
R
```

Then install requiered packages

```
intstall.packages("BiocManager") #type yes so it installs to a local directory 
BiocManager::install("WGCNA") # works best than install.packages("WGCNA") 
q() # to exit R
```
Run script

```
Rscript SCRIPT_NAME.R
```


-- SMH
