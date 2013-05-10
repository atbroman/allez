## If just a list of Probe/Entrez ID's ##

# idlist <- sample(keys(hgu133plus2ENTREZID),100)

list2score <- function(idlist,lib){
    k <- if(substr(lib,1,3)=="org") "SYMBOL" else "ENTREZID"
    gnames <- keys(get(paste(lib,k,sep="")))   
    g <- as.numeric(gnames %in% idlist)
    names(g) <- gnames
    g
  }
