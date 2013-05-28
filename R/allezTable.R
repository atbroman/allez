## Outputs top GO/KEGG categories ##

allezTable <- function(allez.out,
                       n.low=5,
                       n.upp=500,
                       n.cell=0,
                       zthr=5,
                       type=c("gene_id","symbol"),
                       in.set=FALSE){
   type <- match.arg(type)
   ## Number of genes in list and functional set ##
   nc <- sapply(rownames(allez.out$setscores), function(x){
            s <- allez.out$aux$set.data[
                 allez.out$aux$set.data[,1] %in% x,"gscores"]
            sum(s>0 & !is.na(x))})
   ok <- (allez.out$setscores$set.size >= n.low) &
         (allez.out$setscores$set.size <= n.upp) &
         (allez.out$setscores$z.score >= zthr) &
         (nc >= n.cell)
   allez.table <- allez.out$setscores[ok,
                  -grep("sd",names(allez.out$setscores))]
   ## In allez.out$aux$set.data ##
   allez.table$genes <- sapply(rownames(allez.table), function(x){
     y <- allez.out$aux$set.data[allez.out$aux$set.data[,1]==x,]
     paste(y[order(y[,"gscores"],decreasing=
           allez.out$setscores[x,"z.score"]>0),type],collapse=";")})
   if(in.set==TRUE)
     allez.table <- cbind(allez.table, in.set=nc[rownames(allez.table)],
          in.genes=sapply(rownames(allez.table), function(x){
            y <- allez.out$aux$set.data[allez.out$aux$set.data[,1] %in% x,]
            paste(y[y[,4]>0,type][order(y[y[,4]>0,"gscores"],
            decreasing=allez.out$setscores[x,"z.score"]>0)],collapse=";")}))
    ##allez.table$in.set <- allez.table$set.mean*allez.table$n.genes
   ord <- order(allez.table$set.mean,decreasing=TRUE)
   allez.table[ord,]
 }
