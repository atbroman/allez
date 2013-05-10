\name{allezC}
\alias{allezC}

\title{
Concatenate allez output
}
\description{
  Concatenate output of \code{allez} into a list format; must be done with
  the same probe scores, used to combine GO and KEGG outputs together,
  before inputting into \code{allezPlot}.
}
\usage{
allezC(...)
}
\arguments{
  \item{\dots}{2 or more list objects outputted from the \code{allez} function
}
}
\value{A list in the same format as \code{allez} output}

\author{
Aimee Teo Broman
}
\seealso{
  \code{allez}
  \code{allezPlot}
}
\examples{
\dontrun{
data(npc)
scores  <- sqrt(28) * log((1-npc)/(1+npc))

## GO takes a long time ##
npc.go <- allez( scores=scores, lib="hgu133plus2", sets="GO")
npc.kegg <- allez( scores=scores, lib="hgu133plus2", sets="KEGG")

npc.all <- allezC(npc.go,npc.kegg)
}
}
