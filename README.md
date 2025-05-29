# phyloAMR: an R package to perform PHYLOgenetic analysis of AntiMicrobial Resistance 

[![CI](https://github.com/kylegontjes/phyloAMR/actions/workflows/ci.yml/badge.svg)](https://github.com/kylegontjes/phyloAMR/actions/workflows/ci.yml) 
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![Visit Website](https://img.shields.io/badge/Website-Visit-blue)]([\[https://example.com](https://kylegontjes.github.io/phyloAMR/)])


**Description**

phyloAMR is an R package that characterizes the evolutionary history of genome-influenced traits and investigate phenotype-genotype associations. 

While built to evaluate the phylogenetics of antimicrobial resistance, it can be applied to genome-influenced traits, including infection severity and source of infection.

Leveraging corHMM's ancestral state reconstruction, this package characterizes a trait's gain, loss, and continuation across a phylogeny. These inferred ancestral states can adress numerous questions, such as:  
1. How often does an antibiotic-resistant bacterium emerge and spread across a healthcare network?
2. Are there genotypes associated with the gain or loss of antibiotic resistance?  
3. Do compensatory or revertant mutations follow a genome-informed trait? 

Additional phylogenetically-informed algorithms, that do not use ancestral state reconstruction, have been developed:
1. **nearest_neighbor_algorithm()**: Identification of an isolate's nearest neighbor on the phylogenetic tree. Can be used to identify an isolate's nearest neighbor without a feature (i.e., susceptible neighbor of a resistant isolate).  

**How to install phyloAMR**

```
# Install devtools
install.packages("devtools", dep = TRUE)

# Install corHMM & phyloAMR from GitHub
## corHMM - installed because the GitHub version is preferred
devtools::install_github("thej022214/corHMM", force = TRUE, build_vignettes = FALSE)

## phyloAMR
devtools::install_github("kylegontjes/phyloAMR", force = TRUE, build_vignettes = TRUE)
```

**Questions that this tool can address**

| Question | Method | Function(s) | Inputs | Output |
|---|---|---|---|---|
| What are the ancestral states for our trait? | Ancestral state reconstruction using corHMM | asr() | Dataframe with a trait and a phylogenetic tree | Ancestral reconstruction states | 
| Do isolates belong to episodes of trait emergence or spread? | Ancestral state reconstruction and tree traversal of estimates | asr() + asr_cluster_detection() | Dataframe with genome-influenced trait and a phylogenetic tree | Ancestral reconstruction and calls for clusters and singletons | 
| How often does a genome-informed trait emerge and spread across a healthcare network? | Ancestral state reconstruction and tracing of estimates| asr() + asr_cluster_detection() + asr_cluster_analysis() | Dataframe with genome-influenced trait and a phylogenetic tree | Descriptive statistics on cluster calls and phylogenetic singletons | 
| Genetic features associated with trait/lineage emergence/spread? | Ancestral state reconstruction of phenotype and genotype | asr() + synchronous_detection() + synchronous_permutation_test() | Dataframe with genome-influenced trait + genotypes of interest and a phylogenetic tree | Two traits with synchronous episodes of gain or loss 
| What genotypes are subsequently gained/lost upon acquisition of a trait (i.e., potential compensatory or revertant mutations) | Ancestral state reconstruction of phenotype and genotypes | asr() + downstream_gain_loss()  + downstream_permutation_test() | Dataframe with genome-influenced trait + genotypes of interest and a phylogenetic tree | Genotypes classified as downstream mutations from a traits gain event |  
| What descriptive characteristics are associated with phenotypic emergence and spread? | Association testing for characteristics and phenotypic emergence (singletons) and clusters | phyloaware_regression() | Cluster calls from asr_cluster_detection() algorithm and a dataframe with characteristics of interest | Statistical association testing results for characteristics of interest | 
| Who is my isolate's closest related strain (i.e., nearest neighbor)? | Evaluation of phylogenetic distance | nearest_neighbor_algorithm() | Dataframe with genome-influenced trait and a phylogenetic tree | An isolate's nearest neighbor |  
