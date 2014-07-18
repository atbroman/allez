\name{list2score}
\alias{list2score}

\title{
Create a score vector from a list of genes
}
\description{
Given a list of genes, a platform (chip type or organism), and GO or
KEGG category, create a binary vector for input into \code{allez}
}
\usage{
list2score(idlist, lib, idtype = c("SYMBOL", "ENTREZID", "PROBEID", "REFSEQ", "ENSEMBL", "ACCNUM", "UNIPROT", "PMID"), sets = c("GO", "KEGG"))
}

\arguments{
  \item{idlist}{character vector of gene symbols, entrez ID,
  probe ID (for chip platform only), Refseq ID, Ensembl ID, Accession
  number, Uniprot ID, or PMID.
}
  \item{lib}{character string, name of data package corresponding to microarray
  platform or organism (e.g. "hgu133plus2" or "org.Hs.eg", do not
  include ".db")
}
  \item{idtype}{character string, one of "SYMBOL", "ENTREZID", "PROBEID",
  "REFSEQ", "ENSEMBL", "ACCNUM", "UNIPROT", or "PMID"
}
  \item{sets}{character string, "GO" or "KEGG", for type of gene set to
  be used in \code{allez}
}
}
\value{
  Named numeric vector of {0,1}, with names corresponding to all probe_id on the
  platform or gene_id in the organism database that are connected to a
  gene set in GO or KEGG.
}
\author{
Aimee Teo Broman
}
\seealso{
\code{allez}
}
\examples{
library(hgu133plus2.db)
idlist <- sample(mappedRkeys(hgu133plus2SYMBOL),500)

scores <- list2score(idlist,lib="hgu133plus2",idtype="SYMBOL", sets="KEGG")
}
