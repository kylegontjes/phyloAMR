# phyloAMR: an R package to perform PHYLOgenetic analysis of AntiMicrobial Resistance 

[![R CMD Build](https://github.com/kylegontjes/phyloAMR/actions/workflows/ci.yml/badge.svg)](https://github.com/kylegontjes/phyloAMR/actions/workflows/ci.yml) 
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![Visit Website](https://img.shields.io/badge/Website-Visit-blue)]([\[https://example.com](https://kylegontjes.github.io/phyloAMR/)])


**Description**

phyloAMR is an R package that characterizes the evolutionary history of genome-influenced traits and investigate phenotype-genotype associations. This package can be applied to genome-influenced traits, including antibiotic resistance, infection severity, source of infection, and genotypes.  

Leveraging corHMM's ancestral state reconstruction, this package characterizes a trait's gain, loss, and continuation across a phylogeny. Inferred ancestral states can address numerous research questions:  
1. How often does an antibiotic-resistant bacterium emerge and spread across a healthcare network?
2. Are certain genotypes associated with the gain or loss of antibiotic resistance?  
3. Do compensatory or revertant mutations follow the acquisition of a genome-informed trait? 
4. What patient characteristics are associated with the emergence and spread of antibiotic-resistant organisms?

Additional phylogenetically-informed algorithms, that do not use ancestral state reconstruction, have been developed:
1. **nearest_neighbor_algorithm()**: Identification of an isolate's nearest neighbor on the phylogenetic tree. Can be used to identify an isolate's nearest neighbor without a feature (i.e., susceptible neighbor of a resistant isolate). This method can be combined with nearest_neighbor_analysis() and nearest_neighbor_summary_statistic() to compare characteristics of the isolate and their nearest neighbor. This analysis framework is useful for identifying the influence of genetic mutations on bacterial phenotypes (i.e., acquisition of resistance mutation and a resistance measurement).
2. **phyloaware_regression()**: Phylogenetically-informed regression to test for associations between descriptive characteristics and the emergence and spread of a genome-influenced trait.

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
| How often is the trait gained and lost across the phylogeny? | Ancestral state reconstruction and transition statistics | asr() + asr_transition_analysis() | Ancestral reconstruction states from asr() | Descriptive statistics on trait gain, loss, and continuation |
| Do isolates belong to episodes of trait emergence or spread? | Ancestral state reconstruction and tree traversal of estimates | asr() + asr_cluster_detection() | Dataframe with genome-influenced trait and a phylogenetic tree | Ancestral reconstruction and calls for phylogenetic clusters and singletons | 
| How often does a genome-informed trait emerge and spread across a healthcare network? | Ancestral state reconstruction and tracing of estimates| asr() + asr_cluster_detection() + asr_cluster_analysis() | Dataframe with genome-influenced trait and a phylogenetic tree | Descriptive statistics on phylogenetic cluster calls and singletons | 
| Genetic features associated with trait/lineage emergence/spread? | Ancestral state reconstruction of phenotype and genotype | asr() + synchronous_transitions() + synchronous_permutation_test() | Dataframe with genome-influenced trait + genotypes of interest and a phylogenetic tree | Two traits with synchronous episodes of gain or loss |
| What genotypes are subsequently gained/lost upon acquisition of a trait (i.e., potential compensatory or revertant mutations) | Ancestral state reconstruction of phenotype and genotypes | asr() + downstream_transitions()  + downstream_permutation_test() | Dataframe with genome-influenced trait + genotypes of interest and a phylogenetic tree | Genotypes classified as downstream mutations from a traits gain event |  
| What descriptive characteristics are associated with phenotypic emergence and spread? | Association testing for characteristics and trait emergence (singletons) and spread (clusters) | phyloaware_regression() | Cluster calls from asr_cluster_detection() algorithm and a dataframe with characteristics of interest | Statistical association testing results for characteristics of interest | 
| Who is my isolate's closest related strain (i.e., nearest neighbor)? | Evaluation of phylogenetic distance | nearest_neighbor_algorithm() | Dataframe with genome-influenced trait and a phylogenetic tree | An isolate's nearest neighbor |  

**Citation**

Kyle J Gontjes, Aryan Singh, Sarah E Sansom, James D Boyko, Stephen A Smith, Ebbing Lautenbach, Evan Snitkin, Phylogenetic Context of Antibiotic Resistance Provides Insights into the Dynamics of Resistance Emergence and Spread, The Journal of Infectious Diseases, 2025;, jiaf478, https://doi.org/10.1093/infdis/jiaf478
