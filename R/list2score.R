
list2score <- function(idlist,lib,
              idtype=c("SYMBOL","ENTREZID","PROBEID","REFSEQ",
                  "ENSEMBL","ACCNUM","UNIPROT","PMID"),
              sets=c("ALL","GO","KEGG")){

    idtype <- match.arg(idtype)
    sets <- match.arg(sets)

    gstype <- if(substr(lib,1,3)=="org") "EG" else "PROBE"
    gs <- c(GO=paste("GO2ALL",gstype,"S",sep=""),
            KEGG=paste("PATH2",gstype,sep=""))

    allgenes <- mappedRkeys(get(paste(lib,idtype,sep="")))
    if(any(!(i <- idlist %in% allgenes))){
        if(all(i==FALSE)) stop(paste("idlist containing ", idtype,
                   "s not found in ", lib,".db",sep="")) else
        warning(paste(paste(idlist[!i], collapse=", "),
                      " not found in ", lib,".db",sep=""))
        idlist <- idlist[i]
    }
    ## PROBE level if chip.db, EG level if org.db ##
    ids <- if(substr(lib,1,3)=="org" & idtype=="ENTREZID") idlist else
           if(substr(lib,1,3)!="org" & idtype=="PROBEID") idlist else
           toTable(revmap(get(paste(lib,idtype,sep="")))[idlist])[,1]

    ## ALL=all EG or PROBE keys in platform ##
    ## GO, KEGG = all EG or PROBE in GO or PATH sets ##
    allids <- if(sets=="ALL") keys(get(paste(lib,idtype,sep=""))) else
        mappedLkeys(get(paste(lib,gs[sets],sep="")))   
    
    g <- as.numeric(allids %in% ids)
    names(g) <- allids
    g
  }
