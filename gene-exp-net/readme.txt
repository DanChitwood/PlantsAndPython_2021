files for gene coexpression network group

subset.ipynb is my method of making a subset from the metadata & should output a functional bash line for retreiving your subset from the original data. I spent way too much time on this but it should be very easy to use & making it helped me get more comfy with python & pandas.

subset.py is a stripped down script from the subset.ipynb notebook with parameter variables declared at the beginning for convenience
please use it (remember to adjust parameters & output file name)

fpkm_filtering.ipynb shows how I filtered data for WGCNA. the ~filtering by low expression across samples~ might be useful for other groups, just don't use the log transform part 

filter.py can be run on a developmental node on smaller subsets, I suggest making a subset then filering but it might run okay on the full set on a processing node. also worth noting that this does a log2(x+1) transform so if you're adapting it for a different purpose make sure to remove that

metafix.ipynb is just documentation of how I generated the indexed metadata for subsetting. you should be fine to use the indexed_meta.csv file located at /mnt/research/HRT841_F21/Class_project/subsets as input for the subset.ipynb workflow



- SMH
