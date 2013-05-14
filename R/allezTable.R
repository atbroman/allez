## Outputs top GO/KEGG categories ##

allezTable <- function(allez.out,
                       n.low=5,
                       n.upp=500,
                       zthr=5,
                       type=c("gene_id","symbol"),
                       in.set=FALSE){
   type <- match.arg(type)
   ok <- (allez.out$setscores$n.genes >= n.low) &
         (allez.out$setscores$n.genes <= n.upp) &
         (allez.out$setscores$z.score >= zthr)
   allez.table <- allez.out$setscores[ok,
                  -grep("sd",names(allez.out$setscores))]
   ## In allez.out$aux$set.data ##
   allez.table$genes <- sapply(rownames(allez.table), function(x){
     y <- allez.out$aux$set.data[allez.out$aux$set.data[,1]==x,]
     paste(y[order(y[,"gscores"],decreasing=
           allez.out$setscores[x,"z.score"]>0),type],collapse=";")})
   if(in.set==TRUE)
     allez.table <- cbind(allez.table,
       do.call(rbind,lapply(rownames(allez.table),function(x){
          y <- allez.out$aux$set.data[allez.out$aux$set.data[,1]==x,]
          data.frame(in.set=sum(y[,4]>0),
          in.genes=paste(y[y[,4]>0,type][order(y[y[,4]>0,"gscores"],
          decreasing=allez.out$setscores[x,"z.score"]>0)],collapse=";"))})))
    ##allez.table$in.set <- allez.table$set.means*allez.table$n.genes
   ord <- order(allez.table$set.means,decreasing=TRUE)
   allez.table[ord,]
 }
