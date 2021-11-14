### this is written to filter all samples by UMR, then create subsets by ecotype, genotype, & tissue. 
### these subsets are then aggregated into a single subset filtered by all 3 criteria, but this behavior can be altered by changing the merge functions
### please remember to adjust the parameters

umr_thresh=90 # unique mapped rate cutoff, 90 cuts 75% of samples while 80 only cuts 40%
### variables below should hold strings seperated by | operators without spaces. 
### eco & geno are set to str.fullmatch() by default while tissue is set to str.contains()
eco_to_keep="Col-0|Col-4|Landsberg erecta"
geno_to_keep="wild type|Wild type|wild-type|Wild Type"
tissue_to_keep="lea" 


import pandas as pd
import numpy as np

meta=pd.read_csv('indexed_meta.csv', header=0, index_col=0)
meta['UniqueMappedRate'] = meta['UniqueMappedRate'].str.rstrip('%').astype('float') #UMR was held as string with % appended so get rid of that & convert to float

### filter by unique mapped rate. change umr_thresh to your preferred cutoff
meta_by_umr=meta[meta["UniqueMappedRate"] > umr_thresh]
#just for reference, how many samples made it past that cut?
complete=len(meta)
by_umr=len(meta_by_umr)
print("\n")
print(complete, "samples reduced to", by_umr, "by removing samples with unique mapped read rate below",umr_thresh, "%\nTotal number of samples reduced by ",((1-(by_umr/complete))*100),"%.\n")

### filter by ecotype
meta_by_eco=meta_by_umr[meta_by_umr["Ecotype"].astype(str).str.fullmatch(eco_to_keep)]
eco_kept=meta_by_eco["Ecotype"].value_counts().to_frame()
eco_kept=eco_kept.reset_index()
eco_kept.columns = ['ecotype', 'samples']
print("ecotypes kept:\n", eco_kept, "\n") #print what's in the resulting df to check that it worked
### samples left
by_eco=len(meta_by_eco)
print(by_umr, "samples reduced to", by_eco, "after filtering by ecotype\nnumber of samples reduced by ",((1-(by_eco/by_umr))*100),"%.\n\n")

### filter by genotype
meta_by_geno=meta_by_umr[meta_by_umr["Genotype"].astype(str).str.fullmatch(geno_to_keep)]
geno_kept=meta_by_geno["Genotype"].value_counts().to_frame()
geno_kept=geno_kept.reset_index()
geno_kept.columns = ['genotype', 'samples']
print("genotypes kept:\n", geno_kept, "\n") #print what's in the resulting df to check that it worked
### samples left
by_geno=len(meta_by_geno)
print(by_umr, "samples reduced to", by_geno, "after filtering by genotype\nnumber of samples reduced by ",((1-(by_geno/by_umr))*100),"%.\n\n")

### filter by tissue
meta_by_tissue=meta_by_umr[meta_by_umr["Tissue"].astype(str).str.contains(tissue_to_keep)]
tissue_kept=meta_by_tissue["Tissue"].value_counts().to_frame()
tissue_kept=tissue_kept.reset_index()
tissue_kept.columns = ['tissue', 'samples']
print("tissue types kept:\n", tissue_kept,"\n") #print what's in the resulting df to check that it worked
### samples left
by_tissue=len(meta_by_tissue)
print(by_umr, "samples reduced to", by_tissue, "after filtering by tissue\nnumber of samples reduced by ",((1-(by_tissue/by_umr))*100),"%.\n\n")

### aggregate subsets
meta_eco_geno=meta_by_eco.merge(meta_by_geno, how='inner') # merge result is only those shared between the two df
meta_eco_geno_tissue=meta_eco_geno.merge(meta_by_tissue, how='inner')
### samples left
by_ecoxgeno=len(meta_eco_geno)
by_all=len(meta_eco_geno_tissue)
print(by_umr, "samples reduced to", by_all, "after filtering by ecotype, genotype, & tissue\nfiltering by ecotype removed", (by_umr-by_eco), "samples\nfiltering by genotype removed an additional", (by_eco-by_geno), "samples\nfiltering by tissue removed an additional", (by_eco-by_geno-by_ecoxgeno), "samples\nnumber of samples reduced by ",((1-(by_all/by_umr))*100),"%.\n")

#### alternate aggregation code for filtering by 2 criteria rather than 3
#meta_eco_geno=meta_by_eco.merge(meta_by_geno, how='inner') # merge result is only those shared between the two df
### samples left
#by_ecoxgeno=len(meta_eco_geno)
#print(by_umr, "samples reduced to", by_ecoxgeno, "after filtering by ecotype & genotype\nfiltering by ecotype removed", (by_umr-by_eco), "samples\nfiltering by genotype removed an additional", (by_eco-by_geno), "samples\nnumber of samples reduced by ",((1-(by_ecoxgeno/by_umr))*100),"%.\n")

### generate bash cmd
subset=meta_eco_geno_tissue #chonky name is annoying

subset_index=subset['index'] # list the index numbers of subset
subindx2=subset_index+1 # +1 bc awk starts at 1 & first column is geneIDs
indx=subindx2.astype(str).tolist() #make it a list of strings so we can format it how awk likes
awk_index = ["$"+ n for n in indx] #add $ to index for awk (i think this tells it to pick columns? not sure tbh)

print("\n bash code for this subset is:\nawk -F, \'{OFS=\",\";print $1,", ', '.join(awk_index), "}\' gene_FPKM_200501.csv > subset.csv\n") # print the whole bash cmd




