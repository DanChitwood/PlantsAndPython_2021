library(readr)
library(readxl)
library(dplyr)

setwd('/home/beebo_bebop/HRT841 - Plants & Python')

badmeta=read_excel('Arabidopsis_metadata.xlsx')
massive=read_excel('meta_index_corrected.xlsx') # everything
huge=massive[massive$UniqueMappedRate > 90, ] # 90% umr
large=huge[huge$Ecotype == 'Col-0', ] # just Col-0
big=large[large$Genotype == "wild type", ] # just wild type Col-0
bigger=large[large$Genotype == 'wild type' | large$Genotype == '/', ] # wild type or balnk Col-0
#tissues=unique(big$Tissue) # prints unique tissue factors
root=big[big$Tissue == 'root' |
           big$Tissue == 'whole root' |
           big$Tissue == '10 day old roots' |
           big$Tissue == 'seedling root' |
           big$Tissue == 'Root', ] # most roots except tip
leaf=big[big$Tissue == 'leaves' |
           big$Tissue == 'rosette leaf' |
           big$Tissue == 'leaf' |
           big$Tissue == 'rosette' |
           big$Tissue == 'Rosette leaves' |
           big$Tissue == '4-week old leaves' |
           big$Tissue == 'leave' |
           big$Tissue == 'Leaf' |
           big$Tissue == 'Leaves' |
           big$Tissue == 'Rosette leaf' |
           big$Tissue == 'adult vascular leaf', ] # most leaves
# rows for root
rootsub=root$`index`
leafsub=leaf$`index`
rootsub=rootsub+1 #awk index starts at 1, first column of data is gene name
leafsub=leafsub+1
rootugh=paste(as.character(rootsub), collapse=", $") # mostly copypastable to awk function
leafugh=paste(as.character(leafsub), collapse=", $") # mostly copypastable to awk function
rootbash=print(paste0(";print $1,$",rootugh,"}' gene_FPKM_200501.csv > rootsub.csv")) #most of bash
leafbash=print(paste0(";print $1,$",leafugh,"}' gene_FPKM_200501.csv > leafsub.csv"))
cat('{OFS=","',leafbash) #adds more of bash, just needs ~awk -F, '~ added by hand at beginning to function




#lets fix that indexing
# metafix=read_excel('metafix.xlsx') #i struggle with read.csv for some reason & am too lazy to figure out the issue so i just convert to xl
#colnames(metafix)[1]='index' #corrects the index column
#colnames(massive)[1]='indexfkd'
#indexedmeta=merge(massive,metafix,by.x='Sample',by.y='Sample',)
#propermeta=indexedmeta[order(indexedmeta$index), ]
#meta=propermeta%>%relocate(index, .before=indexfkd)
#meta_index=meta[-c(3)]
#meta=meta_index%>%relocate(index, .before=Sample)
#write.csv(meta,file='meta_index_corrected.csv')

#metacheck = massive[seq(1, nrow(massive), 100), ]#creating a subsample of whole dataset to doublecheck meta fixing
#metacheckindex=metacheck$index+1
#metaugh=paste(as.character(metacheckindex), collapse=", $")
#metabash=print(paste0(";print $1,$",metaugh,"}' gene_FPKM_200501.csv > metacheck.csv"))

#dblchk=propermeta[c(1:4,12)] #check back to the merged meta
#dblchk=dblchk[order(dblchk$index), ]
#dbl = dblchk[seq(1, nrow(dblchk), 100), ]
#metacheck$Sample == dbl$Sample
#metacheck=merge(dbl,badmeta,by='Sample')
#metacheck=metacheck[c(1,2,6)]
