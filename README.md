# phylosuite: a package to perform phylogenetically-informed analysis of genome-influenced traits

**Description**

The R package, phylosuite, is a collection of phylogenetically-informed algorithms for analyzing genome-influenced traits, such as antibiotic resistance and infection severity. 

This tool has three main workhorses:
1. **asr_clustering_algorithm()**: An algorithm that leverages ancestral state reconstruction to evaluate the clustering of features across a phylogenetic tree
2. **nearest_neighbor_algorithm()**: An algorithm that permits the identification of an isolate's nearest neighbor on the phylogenetic tree
3. **permute_burden_algorithm()**: An algorithm that implements a permutation test to identify genes with elevated mutation frequency. 

This package's suite of algorithms and workflows was developed to address numerous biological questions, including:  
1. How often does a genome-informed trait emerge and spread across a healthcare network?
2. Do hitchhiking or compensatory mutations follow a genome-informed trait?
3. Do subpopulations have unique or shared mutational pathways?  

**Installation**
```
# Install devtools
install.packages("devtools", dep=TRUE)

# Install phyloaware from github
devtools::install_github("kylegontjes/phylosuite", force=TRUE, build_vignettes=TRUE)
```

**Questions**
| Question | Method | Function(s) | Inputs | Output |
|---|---|---|---|---|
| How often does a genome-informed trait emerge and spread across a healthcare network? | Ancestral state reconstruction and tracing of estimates| asr_clustering_algorithm() | Dataframe with genome-influenced trait and a phylogenetic tree | ancestral reconstruction and calls for clusters and singletons | 
| Who is the nearest neighbor of my isolate | Evaluation of phylogenetic distance | nearest_neighbor_algorithm() | Dataframe with genome-influenced trait and a phylogenetic tree | An isolate's nearest neighbor |  
| What genes are mutated more often than expected? | Permutation testing of observed mutations | permute_burden_algorithm() | VCF file and dictionary with details on gene length | Statistically significant genes with more mutations than expected | 
