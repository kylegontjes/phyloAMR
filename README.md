# phylosuite: a package to perform phylogenetically-informed analysis of genome-influenced traits

**Description**

This R package, phylosuite, is a collection of phylogenetically informed algorithms for analyzing genome-influenced traits, such as antibiotic resistance, infection severity, or source of isolation. 

We were inspired to build this package for numerous reasons. Notably, with an increased appreciation of how bacterial phenotypes emerge and spread across populations, we recognized a noticeable gap in technology. That is the absence of easy-to-use tools that can be used to characterize how traits emerge and spread across a population, alongside how genotypes and phenotypes interact. While genome-wide association studies (GWAS), such as hogwash and treeWAS, can directly or indirectly address these questions, their output does not permit easy downstream manipulations. These tools often rely upon ancestral state reconstruction algorithms, notably ape's ace function or corHMM. Recognizing the power of these algorithms, we wanted to provide a direct-to-consumer algorithm that produces easy-to-use output and addresses several questions that traditional GWAS studies are unable to answer in their current state.
 
This tool's main workhorse is the **asr()**, which leverages corHMM's ancestral state reconstruction algorithm to characterize genome-influenced features' gain, loss, and continuation across a phylogenetic tree. This function requires a phylogenetic tree and dataframe with variables corresponding to the trait of interest and the phylogeny tip labels. After performing joint or ancestral state reconstruction using corHMM, the resultant state predictions of ancestral nodes are parsed to generate a parent-child dataframe. By traversing the phylogenetic tree from the tips to the root, the episodes of trait gain, loss, and continuation are added to the phylogenetic tree's edge matrix.

This information can be leveraged, using our downstream algorithms, to address numerous biological questions, including: 
1. How often does a regime with a genome-informed trait emerge and spread across a healthcare network?
2. Is there evidence of synchronous gain/loss of traits on the phylogenetic tree?
3. Do compensatory or revertant mutations follow a genome-informed trait?
4. Do subpopulations have unique or shared mutational pathways?  

Additional phylogenetically-informed algorithms, that do not use ancestral state reconstruction, have been developed:
1. **nearest_neighbor_algorithm()**: An algorithm that permits the identification of an isolate's nearest neighbor on the phylogenetic tree
2. **permute_burden_algorithm()**: An algorithm implementing a permutation test to identify genes with elevated mutation frequency. 
 
**How to install phylosuite**
```
# Install devtools
install.packages("devtools", dep=TRUE)

# Install phyloaware from GitHub
devtools::install_github("kylegontjes/phylosuite", force=TRUE, build_vignettes=TRUE)
```

**Questions that this tool can address**
| Question | Method | Function(s) | Inputs | Output |
|---|---|---|---|---|
| What are the ancestral states for our genome-influenced trait? | Ancestral state reconstruction using corHMM | asr() | Dataframe with genome-influenced trait and a phylogenetic tree | ancestral reconstruction states | 
| Do isolates belong to episodes of trait emergence or spread? | Ancestral state reconstruction and tracing of estimates | asr() + asr_cluster_detection() | Dataframe with genome-influenced trait and a phylogenetic tree | ancestral reconstruction and calls for clusters and singletons | 
| How often does a genome-informed trait emerge and spread across a healthcare network? | Ancestral state reconstruction and tracing of estimates| asr() + asr_cluster_detection() + asr_cluster_analysis() | Dataframe with genome-influenced trait and a phylogenetic tree | Descriptive statistics on cluster calls and phylogenetic singletons | 
| Genetic features associated with trait/lineage emergence/spread? | Ancestral state reconstruction of phenotype and genotype | asr() + synchronous_detection() + synchronous_permutation_test()  **WIP** | Dataframe with genome-influenced trait + genotypes of interest and a phylogenetic tree | Two traits with synchronous episodes of gain or loss 
| What genotypes are subsequently gained/lost upon acquisition of a trait (i.e., potential compensatory or revertant mutations) | Ancestral state reconstruction of phenotype and genotypes | asr() + downstream_gain_loss()  + downstream_permutation_test()  **WIP** | Dataframe with genome-influenced trait + genotypes of interest and a phylogenetic tree | Genotypes classified as downstream mutations from a traits gain event |  
| What descriptive characteristics are associated with phenotypic emergence and spread? | Association testing for characteristics and phenotypic emergence (singletons) and clusters | phyloaware_regression() | Cluster calls from asr_cluster_detection() algorithm and a dataframe with characteristics of interest | Statistical association testing results for characteristics of interest | 
| What genes are mutated more often than expected in my population of interest? | Permutation testing of observed mutations | permute_burden_algorithm() **WIP** | VCF file and dictionary with details on gene length | Statistically significant genes with more mutations than expected | 
| Who is the nearest neighbor of my isolate? | Evaluation of phylogenetic distance | nearest_neighbor_algorithm() | Dataframe with genome-influenced trait and a phylogenetic tree | An isolate's nearest neighbor |  

**Package Vignettes**

Numerous [vignettes](https://github.com/kylegontjes/phyloaware/tree/master/vignettes) have been constructed to illustrate how to use this package:

1. [Ancestral state reconstruction of a trait](https://github.com/kylegontjes/phylosuite/blob/master/vignettes/ancestral_state_reconstruction_of_a_trait.Rmd)

2. [Detection and characterization of phylogenetic emergence and spread of a trait using ancestral state reconstruction](https://github.com/kylegontjes/phylosuite/blob/master/vignettes/ancestral_state_reconstruction_cluster_detection.Rmd) 

3. [Investigation of genotype and phenotype associations using ancestral state reconstruction](https://github.com/kylegontjes/phylosuite/blob/master/vignettes/ancestral_state_reconstruction_phenotype_genotype_investigations.Rmd)
