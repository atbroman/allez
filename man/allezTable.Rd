\name{allezTable}
\alias{allezTable}

\title{
Extract a table of top-ranked functional sets from allez output
}
\description{
Thresholds allez \code{setscores} element by number of genes in the set
and z-score; orders sets from highest to lowest set mean.
}
\usage{
allezTable(allez.out, n.low = 5, n.upp = 500, zthr = 6, type = c("gene_id", "symbol"), in.set = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allez.out}{list object, output from allez function}
  \item{n.low}{numeric, smallest number of genes in a set to
    consider for ranking}
  \item{n.upp}{numeric, largest number of genes in a set to
    consider for ranking}
  \item{zthr}{numeric, threshold z-scores at this level}
  \item{type}{string, either \code{gene_id} or \code{symbol} will list
    entrez ID or gene symbol, respectively for genes in the set.} 
  \item{in.set}{TRUE or FALSE; if TRUE, will list genes in set that had
    scores > 0}
}
\value{
  A subset of \code{allez.out}$setscores, thresholded by z-score and
  number of genes in a set, ordered by set.means.  \code{genes} column
  is ordered by gene scores.
}
\author{
Aimee Teo Broman and Michael Newton
}
\seealso{
\code{allez}
}
\examples{
data(npc)
scores  <- sqrt(28) * log((1-npc)/(1+npc))
npc.kegg <- allez(scores, lib="hgu133plus2", sets="KEGG")

allezTable(npc.kegg,type="symbol")

}
