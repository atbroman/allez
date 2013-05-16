
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
                   max.n=NULL, ...)
{
  stopifnot(any(!is.na(scores)))
  scorenames <- names(scores)
  
  if (any(duplicated(scorenames))) 
    stop("input IDs must be unique")
  if (any(is.na(scores)))
    warning("scores containing NA's will be ignored")
  if( !is.numeric(scores) ) 
    stop("scores must be numeric")
  vv <- apply( X=as.matrix(scores), MARGIN=2, FUN=var, na.rm=TRUE )
  if( min(vv) == 0 ) stop("no variance in at least one profile" ) 
  
  sets <- match.arg(sets)
  universe <- match.arg(universe)
  transform <- match.arg(transform)
  collapse <- match.arg(collapse)
  setstat <- match.arg(setstat)

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
     n.probes <-  table(unique(set2probe[,c(set_id,"probe_id")])[,1])
   }
  if(collapse != "none"){
    ## Use org info for ENTREZ TO GO/KEGG ID ##
    set2eg <- toTable(getDataEnv(name=ifelse(sets=="GO",
              "GO2ALLEGS","PATH2EG"), lib=orgpkg))
    org.symbol <- toTable(getDataEnv(name="SYMBOL",lib=orgpkg))
    set2eg <- cbind(set2eg, org.symbol[
              match(set2eg$gene_id,org.symbol$gene_id),"symbol",drop=FALSE])
    n.genes <- table(unique(set2eg[,c(set_id, 'gene_id')])[,1])
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
               if (!is.numeric(cutoff))
                 stop("cutoff must be numeric when transform = 'binary'")
                 warning("cutoff used at collapsed gene level not probe level")
               1 * (scores >= cutoff) })
  
## This needs to be changed for GO ?? check old allez code ##
  globe <- switch(sets,GO=gscores,
                  KEGG=gscores)

  mu.globe <- mean(globe)
  sigma.globe <- sd(globe)
  G <- ifelse(is.null(max.n),length(globe),max.n)
  
  set.data <- if(!is.org & collapse=="none"){
         set2probe <- unique(set2probe[,c(set_id,"probe_id","symbol")])
         data.frame(set2probe,gscores[set2probe$probe_id])} else {
           set2eg <- unique(set2eg[,c(set_id, 'gene_id', 'symbol')])
           data.frame(set2eg,gscores=gscores[set2eg$gene_id])}

  set.means <- unlist(tapply(set.data$gscores,set.data[[set_id]],
                      mean,na.rm=TRUE,simplify=FALSE))
  set.sds <- unlist(tapply(set.data$gscores,set.data[[set_id]],
                    sd,na.rm=TRUE,simplify=FALSE))
  set.sizes <- table(set.data[[set_id]])
  
  if (setstat == "mean") {
    ok <- (set.sizes < G)
    dd <- sigma.globe * fact(G=G, m=set.sizes[ok])
    z.score <- (set.means[ok] - mu.globe)/dd
  }
  
  if (setstat == "var") {
    ok <- (set.sizes < G) & (set.sizes > 3)
    E.globe <- fn_getE.Globe(globe=globe)
    sigma1 <- apply(X=as.matrix(set.sizes[ok]), MARGIN=1, FUN=sigma.fun, 
                    E = E.globe, esig2 = sigma.globe^2)
    z.score <- (set.sds[ok]^2 - (sigma.globe^2))/sigma1
    sigma1.Nandi <- apply(X=as.matrix(set.sizes[ok]), MARGIN=1, FUN=sigma.fun.Nandi, 
                          E = E.globe, esig2 = sigma.globe^2)
    z.score.Nandi <- (set.sds[ok]^2 - (sigma.globe^2))/sigma1.Nandi
  }
  
  if (!is.org & collapse == "partial") {
    if (universe == "local")
      warning("partial correction not implemented for `local'", call. = FALSE)
    z.score <- z.score * sqrt(n.genes[ok]/n.probes[ok])
  }
  
  res <- data.frame(set.means=set.means[ok],set.sds=set.sds[ok],
         n.genes=set.sizes[ok], z.score)
  if(!is.org & collapse %in% c("none","partial")) names(res)[3] <- "n.probeset"
  if(!is.org & collapse=="partial") names(res)[4] <- "adjusted z.score"

  aux <- list(set.data = set.data, globe = globe)
  if ((universe == "local") & (sets == "KEGG")) {
    warning("local universe only applicable with GO", call. = FALSE)
  }
  
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
