## Function to combine allez output, for same probescores / probe package ##
## Usually GO and KEGG outputs ##
allezC <- function(...){
  dots <- list(...)
  if(length(dots)==1)
    return(dots[[1]])

  ## Check that entrez-level reductions are the same ##
  globe <- dots[[1]]$aux$globe
  for(i in 2:length(dots)){
    ind <- match(names(globe),names(dots[[i]]$aux$globe))
    if(!identical(globe[!is.na(ind)],dots[[i]]$aux$globe[ind[!is.na(ind)]]))
      stop("The input must be from the same probescores data\n")
  }

  setscores <- set.data <- globe <- NULL
  for(i in 1:length(dots)){
    ss <- switch(names(dots[[i]]$setscores)[1],
          "Term"= dots[[i]]$setscores,
          "path_name"={x <- data.frame(
                       Term = dots[[i]]$setscores[,1],
                       Ontology = "KEGG",
                       dots[[i]]$setscores[,2:5])
                       x})
    setscores <- rbind(setscores,
          ss[is.na(match(rownames(ss),rownames(setscores))),])
    names(dots[[i]]$aux$set.data)[1] <- "set_id"
    set.data <- unique(rbind(set.data, dots[[i]]$aux$set.data))

    globe <- c(globe,
          dots[[i]]$aux$globe[is.na(
              match(names(dots[[i]]$aux$globe),names(globe)))])
  }
  list(setscores=setscores,
       aux=list(set.data=set.data,globe=globe),
       call=match.call())
}

