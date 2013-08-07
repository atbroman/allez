\name{allezMat}
\alias{allezMat}

\title{
Creates an indicator matrix of functional sets vs. genes
}
\description{
Given output list from \code{allez}, this function will create an
indicator matrix of the highest ranking sets; sets IDs are listed as
columns, genes as rows.
}
\usage{
allezMat(allez.out, n.low = 5, n.upp = 500, n.cell = 0, zthr = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allez.out}{list output from \code{allez}}
  \item{n.low}{numeric, smallest number of genes in a set to consider
    for ranking}
  \item{n.upp}{numeric, largest number of genes in a set to consider for
    ranking}
  \item{n.cell}{numeric, minimum number of genes in a set with score>0,
  mainly used in conjunction with binary scores}
  \item{zthr}{numeric, keep sets with z-scores at this level or higher}
}
\value{
A matrix of {0,1}, with genes listed as rows and sets as columns; matrix
element is a 0 if the gene is not in the set, and 1 if the gene is part
of the set.
}
\author{
Aimee Teo Broman and Michael Newton
}
\examples{
data(npc)
scores  <- (1/2)*sqrt(28)*log((1-npc)/(1+npc))
npc.kegg <- allez(scores, lib="hgu133plus2", sets="KEGG")

amat <- allezMat(npc.kegg,zthr=5)
head(amat)
}
