#' Dataframe with relevant data
#'
#' 413 carbapenem-resistant Klebsiella pneumoniae isolates from long-term acute care hospitals across the United States of America. Whole-genome sequencing and antimicrobial susceptibility testing were performed on each isolate
#'
#' @format ## `df`
#' Dataframe with 413 rows and 4 columns
#' @source <https://github.com/Snitkin-Lab-Umich/ltach-crkp-colistin-ms>
"df"

#' Phylogenetic tree
#'
#' A Newick phylogenetic tree of carbapenem-resistant Klebsiella pneumoniae belonging to the epidemic lineage sequence type 258, constructed using IQ-Tree and the GTR+G substitution model, with 413 tips.
#'
#' @format ## `tr`
#' A tree with 413 trips:
#' @source <https://github.com/Snitkin-Lab-Umich/ltach-crkp-colistin-ms>
"tr"

#' Genotype matrix
#'
#' Row names correspond to isolate names, as found in the tr$tip.name and df$tip_name_var. Columns correspond to colistin-related genotypes identified in the aforementioned study. These genotypes can be leveraged to test for synchronous and downstream phenotype-genotype transition events with the synchronous_detection() and downstream_gain_loss() functions, respectively
#'
#' @format ## `genotype_mat`
#' Dataframe with 413 rows and 76 columns
#' @source <https://github.com/Snitkin-Lab-Umich/ltach-crkp-colistin-ms>
"genotype_mat"
