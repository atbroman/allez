\name{allezTable}
\alias{allezTable}

\title{
Extract a table of top-ranked functional sets from allez output
}
\description{
Thresholds allez \code{setscores} element by number of genes in the set
and z-score; orders sets from highest to lowest set mean.  A functional set is
removed if it contains the entire set of input gene scores.
}
\usage{
allezTable(allez.out, n.low = 5, n.upp = 500, n.cell=0, zthr = 5, symbol = FALSE, in.set = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allez.out}{list object, output from allez function}
  \item{n.low}{numeric, smallest number of genes in a set to
    consider for ranking}
  \item{n.upp}{numeric, largest number of genes in a set to
    consider for ranking}
  \item{n.cell}{numeric, minimum number of genes in a set with score>0,
    mainly used for binary scores}
  \item{zthr}{numeric, threshold z-scores at this level}
  \item{symbol}{logical; if TRUE, list gene symbols in set, otherwise
    list gene or probe ID} 
  \item{in.set}{logical; if TRUE, will list genes in set that had
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
scores  <- (1/2)*sqrt(28)*log((1-npc)/(1+npc))
npc.kegg <- allez(scores, lib="hgu133plus2", sets="KEGG")

allezTable(npc.kegg,symbol=TRUE)

}
