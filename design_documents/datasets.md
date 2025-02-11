To demonstrate the utility of this package, we are providing a real-world dataset. Specifically, the package's test data was manipulated from the following pre-print [manuscript](https://doi.org/10.1101/2021.06.11.21258758). This project collected 413 carbapenem-resistant Klebsiella pneumoniae isolates from long-term acute care hospitals in the United States of America. Whole-genome sequencing and antimicrobial susceptibility testing was performed.  

Below, we will describe the avialable data:

1. tr: A Newick formated phylogenetic tree, constructed using IQ-Tree and the GTR+G substitution model, with 413 tips.
2. genotype_mat: Dataframe with 413 rows and 76 columns. Row names correspond to isolate names, as found in the tr$tip.name and df$tip_name_var. Columns correspond to colistin-related genotypes identified in the prior study. For more information about the genotypes, see the preprint and project's GitHub repository.
3. df: Dataframe with 413 rows and 4 columns. Row names correspond to isolates. The following columns were included in the dataset. 
* isolate_no: Isolate name, corresponds to tr$tip.name and rownames in genotype_mat
* Patient_ID: Patient identifier
* clades: If the isolate belongs to clades I, IIA, or IIB of the phylogenetic tree  
* colistin_ns: Non-susceptibility to last line antibiotic, colistin, as inferred by broth microdilution test. A genome-influenced trait that can be leveraged to study bacterial evolution and phenotype-genotype associations, especially with the corresponding genotypic data in the genotype_mat matrix. 