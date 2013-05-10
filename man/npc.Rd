\name{npc}
\alias{npc}
\docType{data}
\title{
Spearman correlations between viral gene EBNA1 and human genes
}
\description{
Spearman correlations between the viral gene EBNA1 and all 54675
host genes (from Sengupta et al 2006) in a study of NasoPharyngeal
Carcinoma

Preliminary data analysis in Sengupta et al indicated something interesting
about host genes that are negatively correlated with the virus.
}
\usage{data(npc)}
\format{
  The format is:
 Named num [1:54675] -0.00806 0.37339 -0.41371 -0.27903 -0.22782 ...
 - attr(*, "names")= chr [1:54675] "1007_s_at" "1053_at" "117_at" "121_at" ...
}
\examples{
data(npc)

\dontrun{
scores  <- sqrt(28)*log((1-npc)/(1+npc))
npc.go <- allez(scores, lib="hgu133plus2")
}
}
\keyword{datasets}
