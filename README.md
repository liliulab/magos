# MAGOS
## Model-based Adaptive Grouping of SNVs. 
 MAGOS performs Model Based Hierarchical Clustering on single nucleotide polymorphism from cancer sequencing samples based on the allele frequency. MAGOS can cluster SNVs and deconvolute (sub)clones from read counts from targeted sequencing, whole-exom sequencing or whole-genome sequencing. MAGOS analysis can be performed on data from single sample or multiple samples. 


## Install

```
library(devtools)
devtools::install_github('liliulab/MAGOS')
```


## Input Data
MAGOS only requires the reference and the alternate read counts from the VCF file. In the the standard VCF format this information can be extracted from the 'AD' (Allelic Depths) section in the FORMAT column. 


### Single Sample
The input data for single sample analysis should be a dataframe or a matrix with two columns. Each row corresponds to a point mutation. The reference reads in the first column and the alterante reads in the second column. 

| | Ref Counts | Alt Counts|
|----|---------|--------|
|mut 1  | ref 1      | alt 1     |
|mut 2 | ref 2      | alt 2     |
|mut 3 | ref 3      | alt 3     |
|...|...|...|


### Multiple Sample
The input data for multiple sample is similar to sigle sample. The input for multiple sample should include the read counts for all the samples. Each row correspond to a SNV and the columns correspond to reference and alternate reads for that SNV in all the samples. Each sample data should be in a pair of columns, eg. data from the first sample in columns 1 and 2, data from the second sample in columns 3 and 4 and .... 

| | Ref Sample 1 | Alt Sample 1| Ref Sample 2| Alt Sample 2| Ref Sample 3 | Alt Sample 3|...|...|
|----|---------|--------|---|---|---|---|---|---|
|mut 1  | ref 1 s1     | alt 1 s1    | ref 1 s2| alt 1 s2| ref 1 s3 | alt 1 s3|...|...|
|mut 2 | ref 2 s1     | alt 2 s1     |ref 2 s2| alt 2 s2| ref 2 s3 | alt 2 s3|...|...|
|mut 3 | ref 3 s1    | alt 3 s1    |ref 3 s2| alt 3 s2| ref 3 s3 | alt 3 s3|...|...|
|...|...|...|...|...|...|...|...|...|


## Usage
```
run = mag.single.run(input.data)


# the results will be in run$results
# the results can be visualized using the following script: 

plot(run$results$vaf.1, run$results$depth.1, col= run$results$colors, xlab= 'VAF', ylab='Depth') 
```

Examples and test data and explanations on different functionalities can be found in the ACE_workshop folder. 




## Reference

The paper in under review in Bioinformatics. 




## Contributors
Algorithm of MAGOS was developed by Li Liu and Navid Ahmadinejad. Please contact liliu at asu.edu for any questions or suggestions. 
