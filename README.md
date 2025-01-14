# phylosuite: a package to perform phylogenetically-informed analysis of genome-influenced traits

**Description**

The R package, phylosuite, is a collection of phylogenetically informed algorithms for analyzing genome-influenced traits, such as antibiotic resistance and infection severity. 

This tool's main workhorse is the **asr_clustering_algorithm()**, which leverages ancestral state reconstruction to characterize the gain, loss, and continuuation of features across a phylogenetic tree. 

This information can be leveraged, using complementary downstream algorithms, to addresss numerous biological questions, including: 
1. How often does a genome-informed trait emerge and spread across a healthcare network?
2. Do hitchhiking or compensatory mutations follow a genome-informed trait?
3. Do subpopulations have unique or shared mutational pathways?  

Additional phylogenetically-informed algorithms have been developed:
1. **nearest_neighbor_algorithm()**: An algorithm that permits the identification of an isolate's nearest neighbor on the phylogenetic tree
2. **permute_burden_algorithm()**: An algorithm implementing a permutation test to identify genes with elevated mutation frequency. 
 
**Installation**
```
# Install devtools
install.packages("devtools", dep=TRUE)

# Install phyloaware from github
devtools::install_github("kylegontjes/phylosuite", force=TRUE, build_vignettes=TRUE)
```

**Questions that this tool can address**
| Question | Method | Function(s) | Inputs | Output |
|---|---|---|---|---|
| Do isolates belong to episodes of trait emergence or spread? | Ancestral state reconstruction and tracing of estimates | asr_clustering_algorithm() | Dataframe with genome-influenced trait and a phylogenetic tree | ancestral reconstruction and calls for clusters and singletons | 
| How often does a genome-informed trait emerge and spread across a healthcare network? | Ancestral state reconstruction and tracing of estimates| asr_clustering_algorithm() + asr_clustering_analysis() | Dataframe with genome-influenced trait and a phylogenetic tree | Descriptive statistics on cluster calls and phylogenetic singletons | 
| Genetic features associated with trait/lineage emergence/spread? | Ancestral state reconstruction of phenotype and genotype | asr_clustering_algorithm() + genotypic_phenotypic_convergence() **WIP** | Dataframe with genome-influenced trait + genotypes of interest and a phylogenetic tree | Genotypes with convergent episodes of gain or loss 
| What descriptive characteristics associated with phenotypic emergence and spread? | Association testing for characteristics and phenotypic emergence (singletons) and clusters | phyloaware_regression() **WIP** | Cluster calls from asr_clustering_algorithm() and dataframe with characteristics of interest | Statistical association testing results for characteristics of interest | 
| What genes are mutated more often than expected in my population of interest? | Permutation testing of observed mutations | permute_burden_algorithm() **WIP** | VCF file and dictionary with details on gene length | Statistically significant genes with more mutations than expected | 
| Who is the nearest neighbor of my isolate | Evaluation of phylogenetic distance | nearest_neighbor_algorithm() | Dataframe with genome-influenced trait and a phylogenetic tree | An isolate's nearest neighbor |  
