
## Make an annotation design matrix
## rows=names(allez.out$aux$globe)
## cols=subset of GO/KEGG categories
allezMat <- function(allez.out,
                     n.low=5,
                     n.upp=500,
                     zthr=3){
## Subset of GO/KEGG terms ##
  ok <- (allez.out$setscores$n.genes >= n.low) &
         (allez.out$setscores$n.genes <= n.upp) &
         (allez.out$setscores$z.score >= zthr)
## allez.out$aux$set.data: 1st col = set id; 2nd col = gene id ##
## mat: genes by GO/KEGG category, 0 if not in cat, 1 if in category ##
  mat <- sapply(rownames(allez.out$setscores)[ok],function(x)
       as.numeric(names(allez.out$aux$globe) %in%
       allez.out$aux$set.data[allez.out$aux$set.data[,1]==x,2]))
  rownames(mat) <- names(allez.out$aux$globe)
  mat
}
