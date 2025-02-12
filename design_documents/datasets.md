To demonstrate the utility of this package, we are providing a real-world dataset. Specifically, the package's test data was abstracted from the following preprint [manuscript](https://doi.org/10.1101/2021.06.11.21258758), with minor changes to fit phylosuite's specifications. This project collected 413 carbapenem-resistant Klebsiella pneumoniae isolates from long-term acute care hospitals in the United States of America. Whole-genome sequencing and antimicrobial susceptibility testing were performed on each isolate, as described by [Han et al., 2019](https://doi.org/10.1128/ac.01622-19) and [Zhang et al., 2023](https://doi.org/10.1017/ice.2022.185), respectively. 

**The following pieces of data are provided in this package:**

**1. tr:** A Newick phylogenetic tree, constructed using IQ-Tree and the GTR+G substitution model, with 413 tips.

**2. genotype_mat:** Dataframe with 413 rows and 76 columns. Row names correspond to isolate names, as found in the tr$tip.name and df$tip_name_var. Columns correspond to colistin-related genotypes identified in the prior study. See the preprint and the project's GitHub repository for more information about the genotypes.

**3. df:** Dataframe with 413 rows and 4 columns. Row names correspond to isolates. The following columns were included in the dataset. 
  * **isolate_no:** Isolate name, corresponds to tr$tip.name and rownames in genotype_mat
  * **Patient_ID:** Patient identifier
  * **clades:** If the isolate belongs to clades I, IIA, or IIB of the phylogenetic tree  
  * **colistin_ns:** Non-susceptibility to a last-line antibiotic, colistin, as inferred by broth microdilution test. A genome-influenced trait that can be leveraged to study bacterial evolution and phenotype-genotype associations, especially with the corresponding genotypic data in the genotype_mat matrix. Of 413 isolates, 136 were non-susceptible to colistin.
