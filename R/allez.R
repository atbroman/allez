
## Need to implement GO / local ##

allez <- function (scores,
                   lib,
                   library.loc=NULL,
                   sets = c("GO","KEGG"),
                   collapse = c("full", "partial", "none"),
                   reduce = NULL,
                   setstat = c("mean", "var"),
                   universe = c("global", "local"),
                   transform = c("none", "binary", "rank", "nscore"),
                   cutoff = NULL,
                   annotate = TRUE,
                   ...)
{
  
  sets <- match.arg(sets)
  universe <- match.arg(universe)
  transform <- match.arg(transform)
  collapse <- match.arg(collapse)
  setstat <- match.arg(setstat)

  stopifnot(any(!is.na(scores)))
  scorenames <- names(scores)

  if (any(duplicated(scorenames))) 
    stop("input IDs must be unique")
  if (any(is.na(scores)))
    warning("scores containing NA's will be ignored")
  if( !is.numeric(scores) ) 
    stop("scores must be numeric")
  vv <- apply( X=as.matrix(scores), MARGIN=2, FUN=var, na.rm=TRUE )
  if( min(vv) == 0 )
    stop("no variance in at least one profile" ) 
  if(universe=="local"){
    if(sets=="KEGG")
      stop("local universe only available for GO")
    if(collapse=="partial")
      stop("partial correction not implemented for 'local'")
  }
  if(transform=="binary" & !is.numeric(cutoff))
    stop("cutoff must be numeric when transform = 'binary'")

  ## Default reduce ##
  uscores <- unique(scores)
  if(is.null(reduce))
   reduce <- if(length(uscores)==2 & identical(uscores,c(0,1))) max else median
  if(length(uscores)==2 & !identical(uscores,c(0,1)))
    warning("if scores are binary, please convert to {0,1}")

  message("Loading necessary libraries...")
  fn_loadSetLibraries( sets=sets )
  fn_loadPlatformLibraries( Libraries=lib, library.loc=library.loc )

  set_id <- switch(sets, GO="go_id", KEGG="path_id")
  is.org <- substr(lib,1,3)=="org"
  orgpkg <- ifelse(is.org,lib[1],get(paste(lib[1],"ORGPKG",sep="")))

  ## ANNOTATION ##
  message("Converting annotations to data.frames ...")
  if(!is.org & collapse %in% c("none","partial")){
     set2probe <- toTable(getDataEnv(name=ifelse(sets=="GO",
                        "GO2ALLPROBES","PATH2PROBE"),lib=lib[1]))
     probe.symbol <- toTable(getDataEnv(name="SYMBOL",lib=lib[1]))
     set2probe <- cbind(set2probe, probe.symbol[
        match(set2probe$probe_id, probe.symbol$probe_id),"symbol",drop=FALSE])
     ## remove probes not in scores vector ##
     set2probe <- set2probe[set2probe$probe_id %in% names(scores),]
   }
  if(collapse != "none"){
    ## Use org info for ENTREZ TO GO/KEGG ID ##
    set2eg <- toTable(getDataEnv(name=ifelse(sets=="GO",
              "GO2ALLEGS","PATH2EG"), lib=orgpkg))
    org.symbol <- toTable(getDataEnv(name="SYMBOL",lib=orgpkg))
    set2eg <- cbind(set2eg, org.symbol[
              match(set2eg$gene_id,org.symbol$gene_id),"symbol",drop=FALSE])
    ## Remove gene_ids not on microarray or scores##
     if(!is.org){
       probe2eg <- toTable(getDataEnv(name="ENTREZID",lib=lib[1]))
       egs <- unique(probe2eg[probe2eg$probe_id %in% names(scores),"gene_id"])
       set2eg <- set2eg[set2eg$gene_id %in% egs,]
     } else set2eg <- set2eg[set2eg$gene_id %in% names(scores),]
  }

  ## SCORES ##
  if(collapse != "none" & !is.org){
      message("Reducing probe data to gene data...")
      pscores <- data.frame(probe2eg,scores=scores[probe2eg$probe_id])
      scores <- unlist(tapply(pscores$scores,as.character(pscores$gene_id),
                       FUN=reduce,simplify=FALSE))
    }
  gscores <- switch(transform,
             none = scores,
             rank = rank(scores),
             nscore = qnorm( (rank(scores))/(length(scores)+1) ),
             binary = {
               warning("cutoff used at collapsed gene level not probe level")
               1 * (scores >= cutoff) })
  
## This needs to be changed for GO ?? check old allez code ##
  globe <- switch(sets,GO=gscores,
                  KEGG=gscores)

  mu.globe <- mean(globe)
  sigma.globe <- sd(globe)
  G <- length(globe)
  
  set.data <- if(!is.org & collapse=="none"){
    set2probe <- unique(set2probe[,c(set_id,"probe_id","symbol")])
    data.frame(set2probe,gscores[set2probe$probe_id])
  } else {
    set2eg <- unique(set2eg[,c(set_id, 'gene_id', 'symbol')])
    data.frame(set2eg,gscores=gscores[set2eg$gene_id])
  }

  set.mean <- unlist(tapply(set.data$gscores,set.data[[set_id]],
                      mean,na.rm=TRUE,simplify=FALSE))
  set.sd <- unlist(tapply(set.data$gscores,set.data[[set_id]],
                    sd,na.rm=TRUE,simplify=FALSE))
  ## ALL set.size < G, see Annotation section ##
  set.size <- table(set.data[[set_id]])
  class(set.size) <- "array"

 if (universe == "global"){
    if (setstat == "mean") {
      dd <- sigma.globe * fact(G=G, m=set.size)
      z.score <- (set.mean - mu.globe)/dd
    }
  
    if (setstat == "var") {
      ok <- set.size > 3
      E.globe <- fn_getE.Globe(globe=globe)
      sigma1 <- sapply(set.size[ok], sigma.fun, 
                    E = E.globe, esig2 = sigma.globe^2)
      z.score <- (set.sd[ok]^2 - (sigma.globe^2))/sigma1
      sigma1.Nandi <- apply(X=as.matrix(set.size[ok]), MARGIN=1,
                          FUN=sigma.fun.Nandi, 
                          E = E.globe, esig2 = sigma.globe^2)
      z.score.Nandi <- (set.sd[ok]^2 - (sigma.globe^2))/sigma1.Nandi
    }
  
    if (!is.org & collapse == "partial") {
      set.np <-  table(unique(set2probe[,c(set_id,"probe_id")])[,1])
      class(set.np) <- "array"
      z.score <- z.score * sqrt(set.size[names(z.score)]/
                                set.np[names(z.score)])
    }
  
    res <- cbind(data.frame(set.mean=set.mean, set.sd=set.sd,
                      set.size=set.size)[names(z.score),],
                      z.score=z.score)
   #? if(!is.org & collapse %in% c("none","partial"))
   #?   names(res)[3] <- "n.probeset"
    if(!is.org & collapse=="partial") names(res)[4] <- "adjusted z.score"
  }

  if(universe=="local"){
    set.names <- if(collapse=="none" & !is.org)
      tapply(set2probe[,1], set2probe[,2],c) else
      tapply(set2eg[,1],set2eg[,2],c)

    go.id <- unique(set.data[,1])

    message("Checking parent/child relationships in GO...")
    ## parent info
    p0 <- rbind(toTable(GOBPPARENTS),
              toTable(GOMFPARENTS),
              toTable(GOCCPARENTS))
    ## only parents annotated on chip ##
    parents <-  p0[p0[,1] %in% go.id & p0[,2] %in% go.id,]
    names(parents)[2] <- "go_parent"
    ## proper subset ##
    in.subset <- apply(parents,1,function(x)
                     all(set.names[[x[1]]] %in% set.names[[x[2]]]))

    parents <- cbind(parents,
        data.frame(set.mean=set.mean[parents[,1]], #v.stat if setstat="mean"
            set.sd=set.sd[parents[,1]], #v.stat^2 if setstat="var"
            set.size=set.size[parents[,1]], #v.nset
            parent.means=set.mean[parents[,2]], #parent.means | v.means
            parent.sds=set.sd[parents[,2]], #parent.sds | v.sds
            parent.size=set.size[parents[,2]]))[in.subset,]
                 #parent.np | n.par | p.np | v.np

    message("Computing enrichment ..." )
    if( setstat == "mean" ){
      den <- parents$parent.sds*fact(parents$parent.size,
                                   parents$set.size)
      den.globe <- sigma.globe*fact(G, parents$set.size) 
      z1 <- (parents$set.mean - parents$parent.means)/den   # local z score
      z0 <- (parents$set.mean - mu.globe )/den.globe   # global z score
      res.full <- data.frame(parents, local.zscore=z1, global.z.score=z0 )
    }

    if( setstat=="var" ){
      E.set <- do.call(rbind,tapply(set.data[,"gscores"],
                                  set.data[,1], fn_getE.Globe))
      colnames(E.set) <- paste("M",1:5,sep="")
      E.par <- E.set[parents[,2],] ## EE
    ## remove some iffy cases
      ok <- (parents$set.size > 3) & (parents$parent.sds > 0) &
        (parents$set.size < parents$parent.size-2 )  ## a bit of room
    ## denominator, local scoring
      den <- apply(cbind(parents$set.size, parents$parent.sds^2, E.par)[ok,],
               1, function(x) sigma.fun(x[1],x[3:7],x[2]))
      z1 <- (parents$set.sd^2 - parents$parent.sds^2)[ok]/den  # local z score

    ## denominator, global scoring
      den.globe <- sapply(parents$set.size[ok],sigma.fun,
                          E=E.globe, esig2=sigma.globe^2)
      z0 <- (parents$set.sd^2[ok] - sigma.globe^2 )/den.globe # global z score
    ## Return full results, use local.max() to pull out "best" local.zscore ##
      res.full <- data.frame(parents[ok,], local.zscore=z1,
                                 global.zscore=z0 )
    }
    res <- local.max(res.full)
  }

  aux <- list(set.data = set.data, globe = globe)
  if(universe=="local") aux$res.full <- res.full
  
  if (!annotate) {
    main <- res
  }
  if (annotate) {
    message("Labeling output ...")
    if (sets == "GO") {
      gterms <- toTable(GOTERM)
      main <- data.frame(gterms[match(rownames(res),gterms$go_id),
                                c("Term","Ontology")],res)
      rownames(main) <- rownames(res)
    }
    if (sets == "KEGG") {
      kterms <- toTable(KEGGPATHID2NAME)
      main <- data.frame(kterms[match(rownames(res),kterms$path_id),
              "path_name",drop=FALSE], res)
      rownames(main) <- rownames(res)
    }
  }
  out <- list(setscores = main, aux = aux, call = match.call())
  return(out)
}
