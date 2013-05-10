allez <- function (scores,
                   lib,
                   library.loc = NULL,
                   sets = c("GO","KEGG"),
                   collapse = c("full"),
                   reduce = median,
                   setstat = c("mean", "var"),
                   universe = c("global", "local"),
                   transform = c("none", "binary", "rank", "nscore"),
                   cutoff = NULL,
                   annotate = TRUE, ...)
{
  stopifnot(any(!is.na(scores)))
  
  if (any(duplicated(names(scores)))) 
    stop("input IDs must be unique")
  if (any(is.na(scores)))
    warning("scores contain NA's (will be ignored)")
  if( !is.numeric(scores) ) 
    stop("scores must be numeric")
  vv <- apply( X=as.matrix(scores), MARGIN=2, FUN=var, na.rm=TRUE )
  if( min(vv) == 0 ) stop("no variance in at least one profile" ) 
  
  sets <- match.arg(sets)
  universe <- match.arg(universe)
  transform <- match.arg(transform)
  collapse <- match.arg(collapse)
  setstat <- match.arg(setstat)
  
  message("Loading necessary libraries...")
  fn_loadSetLibraries( sets=sets )
  fn_loadPlatformLibraries( Libraries=lib )
  
  message("Converting annotations to lists ...")

  is.org <- substr(lib,1,3)=="org"
  orgpkg <- ifelse(is.org,lib[1],get(paste(lib[1],"ORGPKG",sep="")))

  ## Use org info for ENTREZ TO GO/KEGG ID ##
  set2eg <- switch(sets,
            GO = toTable(getDataEnv(name="GO2ALLEGS",lib=orgpkg)),
            KEGG = toTable(getDataEnv(name="PATH2EG",lib=orgpkg)))
  org.symbol <- toTable(getDataEnv(name="SYMBOL",lib=orgpkg))
  set2eg <- cbind(set2eg,
                  org.symbol[match(set2eg$gene_id,org.symbol$gene_id),
                  "symbol",drop=FALSE])
  ## Remove gene_ids not on microarray ##
  if(!is.org){
    probe2eg <- toTable(getDataEnv(name="ENTREZID",lib=lib[1]))
    egs <- unique(probe2eg$gene_id)
    set2eg <- set2eg[set2eg$gene_id %in% egs,]
  }

  set_id <- switch(sets, GO="go_id", KEGG="path_id")
  n.genes <- table(unique(set2eg[,c(set_id, 'gene_id')])[,1])
  n.probes <- if(is.org) NA else 
              table(unique(toTable(getDataEnv(
                   name=ifelse(sets=="GO","GO2ALLPROBES","PATH2PROBE"),
                   lib=lib[1]))[,c(set_id,"probe_id")])[,1])

  ## I think do a special collapse for "none" 
  
  if( collapse == "full" | is.org)
  {
    ## Reduce to gene level from probe level
    message("Reducing probe data to gene data...")

    scores.df <- data.frame(scores)
    scores.df[[paste(ifelse(is.org,"gene","probe"),"id",sep="_")]] <-
      rownames(scores.df)

    if(!is.org){
      probe2eg <- cbind(probe2eg, scores.df[
            match(probe2eg$probe_id,scores.df$probe_id),"scores",drop=FALSE])
      gscores.matrix <- aggregate(x=probe2eg$scores,
            by=list(as.character(probe2eg$gene_id)), FUN=reduce)
    } else gscores.matrix <- scores.df[,2:1]
    gscores <- gscores.matrix[,2]
    names(gscores.matrix) <- c('gene_id', 'gscores')
    names(gscores) <- gscores.matrix$gene_id
    
    gscores <- switch(transform,
               none = gscores,
               rank = rank(gscores),
               nscore = qnorm( (rank(gscores))/(length(gscores)+1) ),
               binary = {  
                 if (!is.numeric(cutoff))
                  stop("cutoff must be numeric when transform = 'binary'")
                  warning("cutoff used at collapsed gene level not probe level")
                  1 * (gscores >= cutoff) })
  }

## This needs to be changed for GO ?? check old allez code ##
  globe <- switch(sets,GO=gscores,
                  KEGG=gscores)

  mu.globe <- mean(globe)
  sigma.globe <- sd(globe)
  G <- length(globe)      ## This line corresponds to line # 178 of WorkingScript_Allez.R
  ## WHy isn't this already unique?  I don't know ##
  set2eg <- unique(set2eg[,c(set_id, 'gene_id', 'symbol')])
  setdata <- cbind(set2eg,
          gscores.matrix[match(set2eg$gene_id,gscores.matrix$gene_id),
                "gscores",drop=FALSE])

  #setdata <- setdata[order(setdata[,set_id], setdata$symbol),]
  
  set.means <- fn_getSetStats(SetData=setdata, Stat=mean, set_id=set_id)
  set.sds <- fn_getSetStats(SetData=setdata, Stat=sd, set_id=set_id)
  set.sizes <- fn_getSetStats(SetData=setdata, Stat=length, set_id=set_id)
  
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
  
## Where does this go?  n.probes is not defined ##  
  if (collapse == "partial") {
    z.score <- z.score * sqrt(n.genes[ok]/n.probes[ok])
    if (universe == "local") {
      warning("partial correction not implemented for `local'", 
              call. = FALSE)
    }
  }
  clabs <- switch(collapse,
        none = c("set.means", "set.sd", "n.probesets", "z.score"),
        full = c("set.means", "set.sd", "n.genes", "z.score"),
        partial = c("set.means", "set.sd", "n.probesets", "adjusted z.score"))
  
  res <- cbind(set.means[ok], set.sds[ok], set.sizes[ok], z.score)
  dimnames(res)[[2]] <- clabs
  aux <- list(set.data = setdata, globe = globe)
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
