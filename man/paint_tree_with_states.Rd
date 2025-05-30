% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paint_tree_with_states.R
\name{paint_tree_with_states}
\alias{paint_tree_with_states}
\title{Paint a phylogenetic tree with ancestral and tip states}
\usage{
paint_tree_with_states(
  parent_child_df,
  tr,
  scale_breaks = c(1, 0),
  scale_labels = c("Present", "Absent"),
  scale_colors = c("red", "black"),
  legend_name = "Trait"
)
}
\arguments{
\item{parent_child_df}{Parent child dataframe from asr()}

\item{tr}{Phylogenetic tree}

\item{scale_breaks}{Values for breaks (marginal = c(1,0.5,0); joint = c(1,0))}

\item{scale_labels}{Labels for present and absent, marginal should be c('Present','Unsure','Absent')}

\item{scale_colors}{Colors for states}

\item{legend_name}{Name for legend}
}
\value{
Tree object with painted states
}
\description{
Function to paint a tree with ancestral and tip statess
}
