# phylosuite: a package to perform phylogenetically-informed analysis of genome-influenced traits

**Description**

This R package, phylosuite, is a collection of phylogenetically informed algorithms for analyzing genome-influenced traits, such as antibiotic resistance, infection severity, or source of isolation. 

This tool's main workhorse is the **asr()**, which leverages corHMM's ancestral state reconstruction algorithm to characterize genome-influenced features' gain, loss, and continuation across a phylogenetic tree. 

This information can be leveraged, using our downstream algorithms, to address numerous biological questions, including: 
1. How often does a regime with a genome-informed trait emerge and spread across a healthcare network?
2. Is there evidence of synchronous gain/loss of traits on the phylogenetic tree?
3. Do hitchhiking or compensatory mutations follow a genome-informed trait?
4. Do subpopulations have unique or shared mutational pathways?  

Additional phylogenetically-informed algorithms, that do not use ancestral state reconstruction, have been developed:
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
| What are the ancestral states for our genome-influenced trait? | Ancestral state reconstruction using corHMM | asr() | Dataframe with genome-influenced trait and a phylogenetic tree | ancestral reconstruction states | 
| Do isolates belong to episodes of trait emergence or spread? | Ancestral state reconstruction and tracing of estimates | asr() + asr_cluster_detection() | Dataframe with genome-influenced trait and a phylogenetic tree | ancestral reconstruction and calls for clusters and singletons | 
| How often does a genome-informed trait emerge and spread across a healthcare network? | Ancestral state reconstruction and tracing of estimates| asr() + asr_cluster_detection() + asr_cluster_analysis() | Dataframe with genome-influenced trait and a phylogenetic tree | Descriptive statistics on cluster calls and phylogenetic singletons | 
| Genetic features associated with trait/lineage emergence/spread? | Ancestral state reconstruction of phenotype and genotype | asr() + convergence_detection() + convergence_permutation() **WIP** | Dataframe with genome-influenced trait + genotypes of interest and a phylogenetic tree | Genotypes with convergent episodes of gain or loss 
| What genotypes are jointly gained/lost (hitchhikers) or subsequently gained/lost upon acquisition of a trait (compensatory mutation) | Ancestral state reconstruction of phenotype and genotypes | asr() + downstream_gain_loss()  + convergence_permutation() **WIP** | Dataframe with genome-influenced trait + genotypes of interest and a phylogenetic tree | Genotypes classified as hitchhikers or compensatory mutations |  
| What descriptive characteristics are associated with phenotypic emergence and spread? | Association testing for characteristics and phenotypic emergence (singletons) and clusters | phyloaware_regression() **WIP** | Cluster calls from asr_cluster_detection() algorithm and a dataframe with characteristics of interest | Statistical association testing results for characteristics of interest | 
| What genes are mutated more often than expected in my population of interest? | Permutation testing of observed mutations | permute_burden_algorithm() **WIP** | VCF file and dictionary with details on gene length | Statistically significant genes with more mutations than expected | 
| Who is the nearest neighbor of my isolate | Evaluation of phylogenetic distance | nearest_neighbor_algorithm() | Dataframe with genome-influenced trait and a phylogenetic tree | An isolate's nearest neighbor |  