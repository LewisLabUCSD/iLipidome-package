% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/species_sub_extract.R
\name{species_sub_extract}
\alias{species_sub_extract}
\title{maps species substructures in each pathway with fold changes from the
“unprocessed_data_test” result and extracts them}
\usage{
species_sub_extract(
  lipid_substructure,
  unprocessed_data_result,
  type = "species",
  pct_limit = 0.3,
  exo_lipid = NULL
)
}
\arguments{
\item{lipid_substructure}{Lipid species substructure table. Output of
"species_sub_transform".}

\item{unprocessed_data_result}{Differential expression for unprocessed
lipidomics data. Output of "unprocessed_data_test".}

\item{type}{"species" or "class" to indicate which type of substructure.}

\item{pct_limit}{A threshold of non-NA percent (0~1) in one biosynthetic
route can be set to control decompostion process and reduce artifacts. If
the proportion of non-missing substructures exceeds the threshold in one
biosynthetic route, the target lipid species will not be decomposed to
substructures.}

\item{exo_lipid}{An character vector specifying exogenous lipid addition.
The exogenous lipids and their adjacent nodes in FA network will not be
decomposed.}
}
\value{
lipid species
}
\description{
maps species substructures in each pathway with fold changes from the
“unprocessed_data_test” result and extracts them
}
