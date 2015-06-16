\name{allez}
\alias{allez}

\title{ Random-set calibration of gene-set statistics}
\description{
  For each set (i.e. category)  in a collection defined by Gene Ontology
  (GO) or the Kyoto encyclopedia of genes and genomes (KEGG), 
 \code{allez} computes a standardized score that measures
 how unusual microarray measurements in that set are compared
 to measurements in same-sized random subsets of the microarray-level data.
 \code{allez} may be used to assess the enrichment of a category
 for genes that are interesting, for example, owing to differential
 expression across groups or owing to substantial correlation with
 some other phenotype.  Additionally \code{allez} may be used to
 identify sets that show unusual variance characteristics.
}
\usage{
allez(scores, lib, idtype = c("ENTREZID", "SYMBOL"),library.loc = NULL, sets = c("GO", "KEGG"), 
locallist = NULL, collapse = c("full","partial","none"), reduce = NULL, setstat = c("mean", "var"),
universe = c("global", "local"), transform = c("none", "binary", "rank",
"nscore"), cutoff = NULL, annotate = TRUE, ...)
}

\arguments{
  \item{scores}{numeric vector of microarray-level or
  organism-level scores,  usually measuring differential expression
  among groups or the relationship of expression values to some other
  variable: possibly log fold change, t-statistic, indicator of
  significant differential expression (i.e. gene list), 
  posterior probability of differential expression, or correlation
  with some phenotype.
  For microarray-level scores, the vector must be named with
  manufacturer probe names; for example, \code{scores} must be 
  named by the Affymetrix probe set IDs.
  For organism-level scores, the vector must be named with Entrez Gene
  Identifiers. 
  }
  \item{lib}{character string, name of data package corresponding to
  microarray (e.g. \code{"hgu133plus2"}), or organism
  (\code{"org.Hs.eg"}). A vector of package names is allowed in case
  multiple chips are used (e.g \code{"moe430a"} and \code{"moe430b"} 
  cover the mouse genome.) }
\item{idtype}{idtype could be either \code{"ENTREZID"} (default) or 
	\code{"SYMBOL"}. It should match
names of the scores vector.}
\item{library.loc}{
%%     ~~Describe \code{library.loc} here~~
}
\item{sets}{character string, describing the collection of
  sets. \code{"GO"} (default) or \code{"KEGG"}}
\item{locallist}{list contains in-house gene sets of interest. Default is NULL.
Each element in the list represents a gene set. Each element should be a character
string containing genes' entrez IDs (if \code{idtype="ENTREZID")} or
gene symbols (if \code{idtype="STMBOL")}. Element names will be used as 
set names. If \code{locallist} is not NULL, \code{universe="local"} will be disabled.}
\item{collapse}{character string, describing the method for reducing
  probe-level data to gene-level data. \code{"full"} (default) uses the
  function \code{reduce} to reduce from 
  probe level to gene level. This uses the
  \code{ENTREZID} class associated with the data
  package/s.  \code{"none"} probe level results
  \code{"partial"} adjustment to the z-scores computed in \code{"full"}
  by the factor sqrt{ n.genes/n.probes }. }

\item{reduce}{ function used in \code{collapse}. }

\item{setstat}{ \code{"mean"} or \code{"var"} indicating the function
   used to compute unstandardized set scores from microarray level
   scores. }

 \item{universe}{ If \code{"global"}, each set's score is compared to
   the score of a random set taken from all annotated genes (if
   \code{collapse="full"} or scores are organism-level) or
   probes (otherwise). If \code{"local"}, comparison is made not to the
   entire set of annotated probes/genes but to those defining different
   parents in the GO DAG (so it's not applicable with KEGG).  With
   multiple parents per set, the one with the largest Z score is
   reported.}

 \item{transform}{ optional transformation of microarray-level data.
    If \code{"binary"}, a cutoff must
    also be supplied. Transformation is done after collapsing to
    gene level if \code{collapse="full"}.   \code{"nscore"} is the
    normal scores transformation \code{ qnorm( rank( ./(G+1) ) ) },
    which makes the gene set means very close to normal, and thus
    improves the z-score quality  }

  \item{cutoff}{ numeric cutoff when \code{transform = "binary"};
      selected genes are those with probe/gene score larger than (or equal to)
       \code{cutoff} }
    
  \item{annotate}{ logical, whether to include set names in the output }

  \item{\dots}{ Other arguments, passed on \code{reduce}. }
}

\details{
\code{allez} uses formulas for both the expected value and
 variance of a sample mean or sample variance computed on
 a random subset of fixed microarray-level data. These formulas
 enable it to standardize observed scores computed on 
 categories from GO or KEGG, and thus to make the different
 categories comparable in terms of how unusual they are compared
 to random sets.  Facilities allow various microarray-level
 scores, various set-level scoring methods (mean, variance), 
 various reductions to gene-level, and various calibrations 
(i.e. for GO, should we use the whole annotated collection of 
  genes/probes or just those in a GO parent).
}
\value{returns a list, containing components:
  \item{setscores}{data frame with rows for
  gene sets and columns for summary information scoring these
  sets. Information includes set name (if \code{annotate = TRUE}),
  the set sample mean ( or set sample variance if \code{setstat="var"} ),
 the set size (number of probes if \code{collapse="none"} and number
 of genes otherwise), and the appropriate Z-score (nominally distributed
 as a standard normal variate). If a functional set contains the entire
 set of inputted gene scores, the z-score will output NA.}
\item{aux}{a list, with auxiliary information from the calculation, 
 including a data frame \code{set.data}, which contains, in table
 format, the map of each gene to a gene
 set and the gene score. The \code{globe} vector is the full complement
 of annotated probes (or genes) (this is the `universe' if
 \code{universe="global"}).}

 If  \code{universe="local"} and \code{sets="GO"}, \code{aux} also
 contains a matrix recording set-level results for all parents of
 every child set (the \code{setscores} in this case reports only the
 largest Z-score among all the parents).

}
\references{Newton, M.A., Quintana, F.A., den Boon, J.A., Sengupta, S.,
  and Ahlquist, P. (2007). Random-set methods identify
  distinct aspects of the enrichment signal in gene-set analysis.
  Annals of Applied Statistics, 1, 85-106.

  Sengupta, S., den Boon, J.A., Chen, I.H., Newton, M.A. et al. (2006).
  Genome-wide expression profiling reveals EBV-associated inhibition 
  of MHC class I expression in nasopharyngeal carcinoma. Cancer Research, 66,
  7999-8006.}


\author{ Michael Newton, Deepayan Sarkar, Aimee Teo Broman, Subhrangshu Nandi}

\examples{
data(npc)
scores  <- (1/2)*sqrt(28)*log((1-npc)/(1+npc))

npc.kegg <- allez( scores=scores, lib="hgu133plus2", sets="KEGG")
allezTable(npc.kegg)[,1:4]
}
