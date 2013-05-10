\name{list2score}
\alias{list2score}

\title{
Create a score vector from a List of IDs
}
\description{
Given a list of probe or gene IDs and a platform (chip type or
organism), create a binary vector for input into \code{allez}
}
\usage{
list2score(idlist, lib)
}

\arguments{
  \item{idlist}{character vector, probe_id or gene_id}
  \item{lib}{character, name of data package correpsonding to microarry
    platform or organism (e.g. "hgu133plus2" or "org.Hs.eg")}
}
\value{
Named numeric vector of {0,1}, with names corresponding to all probe_id on the
platform or gene_id in the organism database.
}
\author{
Aimee Teo Broman
}

\examples{
library(hgu133plus2.db)
ids <- sample(keys(hgu133plus2SYMBOL),500)

scores <- list2score(ids,lib="hgu133plus2")
}
