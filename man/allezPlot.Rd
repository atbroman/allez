\name{allezPlot}
\alias{allezPlot}

\title{
Display an image of gene scores by functional sets
}
\description{
Plots the highest ranking functional sets against genes, with highest
gene scores plotted first.  Gene scores are indicated by darker color.
}
\usage{
allezPlot(allez.out, n.low = 5, n.upp = 500, zthr = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allez.out}{list output from \code{allez} function, or
    concatenated output from \code{allezC}.}
  \item{n.low}{numeric, smallest number of genes in a set to consider
    for ranking}
  \item{n.upp}{numeric, largest number of genes in a set to consider for
    ranking} 
  \item{zthr}{numeric, keep sets with z-scores at this level or higher}
  \item{\dots}{any extra arguments for \code{image}, \code{axis}, or
    \code{text}}
}
\details{
\code{allezPlot} creates an indicator matrix of the highest ranking sets
using \code{allezMat}, adds gene score information, and orders the genes
by highest ranking set and gene score.
}
\author{
Aimee Teo Broman and Michael Newton
}
\seealso{
  \code{allezMat}
  \code{allez}
}
\examples{
data(npc)
scores  <- sqrt(28) * log((1-npc)/(1+npc))

npc.kegg <- allez( scores=scores, lib="hgu133plus2", sets="KEGG")

allezPlot(npc.kegg,zthr=5)

\dontrun{
npc.go <- allez( scores=scores, lib="hgu133plus2", sets="GO")
npc.all <- allezC(npc.go, npc.kegg)

allezPlot(npc.all,zthr=6)
}
}