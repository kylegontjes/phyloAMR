**Real dataset for answering a biological question using the tool**

The package's data were abstracted from the following preprint [manuscript](https://doi.org/10.1101/2021.06.11.21258758), with minor changes to fit phylosuite's specifications. This project collected 413 carbapenem-resistant Klebsiella pneumoniae isolates from long-term acute care hospitals across the United States of America. Whole-genome sequencing and antimicrobial susceptibility testing were performed on each isolate, as described by [Han et al., 2019](https://doi.org/10.1128/ac.01622-19) and [Zhang et al., 2023](https://doi.org/10.1017/ice.2022.185), respectively. 

This dataset includes an antibiotic susceptibility profile for the last-resort antibiotic, colistin, and a collection of hand-curated resistance-associated genotypes. This will permit the characterization of the evolutionary history of these traits, alongside the evaluation of phenotype-genotype associations using our ancestral state reconstruction-informed algorithms. 

Specifically, we can ask the following question: 
1. How many lineages of colistin resistant organisms exist in this population? What is the antibiotic's dynamics of emergene and spread?
2. What genotypes are associated with resitance emergence and spread? Are there any genotypes associated with phylogenetic emergence or are there any downstream compensatory mutations? 

Knowing what we know about colistin resistance, we expect for there to be numerous emergence and spread events. Given the cost of colistin resistance, I anticipate several reversion events and the identification of genotype associated with both gain (i.e., synchronous gain events) and loss (i.e., synchronous trait loss and genotype gain events). I anticipate that this tool will help us better characterize our genotypes in this paper using a data-driven approach. 

**Attached dataset for using the tool**

**1. tr:** A Newick phylogenetic tree of carbapenem-resistant _Klebsiella pneumoniae_ belonging to the epidemic lineage sequence type 258, constructed using IQ-Tree and the GTR+G substitution model, with 413 tips. 

**2. genotype_mat:** Dataframe with 413 rows and 76 columns. Row names correspond to isolate names, as found in the tr$tip.name and df$tip_name_var. Columns correspond to colistin-related genotypes identified in the aforementioned study. These genotypes can be leveraged to test for synchronous and downstream phenotype-genotype transition events with the synchronous_detection() and downstream_gain_loss() functions, respectively. 

For more information about the genotypes, please review the preprint and the project's GitHub repository (https://github.com/Snitkin-Lab-Umich/ltach-crkp-colistin-ms).

**3. df:** Dataframe with 413 rows and 4 columns. Row names correspond to isolates. The following columns were included in the dataset. 
  * **isolate_no:** Isolate name, corresponds to tr$tip.name and rownames in genotype_mat
  * **Patient_ID:** Patient identifier
  * **clades:** If the isolate belongs to clades I, IIA, or IIB of the phylogenetic tree  
  * **colistin_ns:** Non-susceptibility to a last-line antibiotic, colistin, as inferred by broth microdilution test. Of 413 isolates, 136 were non-susceptible to colistin. This genome-influenced trait can be leveraged to study bacterial evolution and phenotype-genotype associations, especially given the corresponding genotypic data in the genotype_mat matrix. 
