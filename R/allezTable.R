## Outputs top GO/KEGG categories ##

allezTable <- function(allez.out,
                       n.low=5,
                       n.upp=500,
                       n.cell=0,
                       zthr=5,
                       symbol=FALSE,
                       in.set=FALSE){
  ## gene list of gene_id, probe_id, or symbol, from set.data ##
  idcol <- ifelse(symbol,3,2)
  ## z.score column ##
  zcol <- grep("z.score",colnames(allez.out$setscores))[1]

   ## Number of genes in list and functional set ##
  nc <- tapply(allez.out$aux$set.data$gscores,
               allez.out$aux$set.data[,1],
               function(x) sum(x>0 & !is.na(x)))
  G <- length(allez.out$aux$globe)

  ## If set.size==G then z.score=NA ##
  ok <- (allez.out$setscores$set.size >= n.low) &
    (allez.out$setscores$set.size <= n.upp) &
      (allez.out$setscores$set.size < G) &
      (allez.out$setscores[,zcol] >= zthr) &
      (nc[rownames(allez.out$setscores)] >= n.cell)
  allez.table <- allez.out$setscores[ok,
                 -grep("sd",names(allez.out$setscores))]
   
   ## Subset set.data ##
  set.data <- allez.out$aux$set.data[
              allez.out$aux$set.data[,1] %in% rownames(allez.table),]
  set.data <- set.data[order(set.data$gscores,decreasing=TRUE),]

   ## rownames(genes) == rownames(allez.table) ##
  genes <- data.frame(
             pos=tapply(set.data[,idcol],set.data[,1],paste,collapse=";"),
             neg=tapply(set.data[,idcol],set.data[,1],function(x)
               paste(rev(x),collapse=";")))
  allez.table$genes <- if(nrow(allez.table)>0) genes[cbind(1:nrow(allez.table),
      ifelse(allez.table[,grep("z.score",names(allez.table))[1]]>0,1,2))] else
      character(0)

  if(in.set==TRUE){
    set.data <- set.data[set.data$gscores>0,]
    genes <- data.frame(
             pos=tapply(set.data[,idcol],set.data[,1],paste,collapse=";"),
             neg=tapply(set.data[,idcol],set.data[,1],function(x)
               paste(rev(x),collapse=";")))
    allez.table <- cbind(allez.table,
                    in.set=nc[rownames(allez.table)],
                    in.genes=if(nrow(allez.table)>0)
                       genes[cbind(1:nrow(allez.table),
                       ifelse(allez.table[,grep("z.score",
                         names(allez.table))[1]]>0,1,2))] else
                          character(0))
   }
    ##allez.table$in.set <- allez.table$set.mean*allez.table$n.genes
  ord <- order(allez.table$set.mean,decreasing=TRUE)
  allez.table[ord,]
 }
