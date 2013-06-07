## This is the function library for allez
################################################################################

################################################################################
#### This function loads required "set" libraries for allez                 ####
################################################################################
fn_loadSetLibraries <- function(sets){
  require(annotate)
  options(warn=-1)
  switch(sets, GO = {
    if( !try(require(GO.db), silent=TRUE) ){
      stop( "Please install 'GO.db' " )
    } else{
      require(GO.db)
      message( 'Loaded Package GO.db' )
    }
  }, KEGG = {
    if(!try( require(KEGG.db), silent=TRUE )){
      stop( "Please install 'KEGG.db' " )
    } else{
      require(KEGG.db)
      message( 'Loaded Package KEGG.db' )
    }
  }, REACTOME = {
    if(!try(require(reactome.db), silent=TRUE)){
      stop( "Please install 'reactome.db' " )
    } else{
      require(reactome.db)
      message( 'Loaded Package reactome.db' )
    }
  },
         {
           stop( "Entry should be GO, KEGG or REACTOME, all uppercase")
         })
  options(warn=1)
}
################################################################################

################################################################################
#### This function loads required "platform" libraries for allez            ####
################################################################################
fn_loadPlatformLibraries <- function(Libraries, library.loc=library.loc){
  for(lib in Libraries){
    if(substr(x=lib, start=nchar(lib)-2, stop=nchar(lib)) == '.db'){
      lib2 <- lib
    } else{
      lib2 <- paste(lib, 'db', sep='.')
    }
    if(!try(require(package=lib2, character.only = TRUE,
                    lib.loc=library.loc), silent=TRUE)){ 
      stop( paste('Please install', lib2) )
    } else{
      require(package=lib2, character.only = TRUE, lib.loc=library.loc)
      message( paste( 'Loaded Package', lib2 ))
    }
    
  }
}
  
################################################################################
#### This function pastes annotation library and map name                   ####
################################################################################
getDataEnv <- function(lib, name) {
  get( x=paste(lib, name, sep = ""), mode = "any")
}

################################################################################
#### Some functions for modified standard deviation calculation             ####
################################################################################
fn_getE.Globe <- function(globe){
  G <- length(globe)
  T1 <- sum(globe^4)
  T2 <- sum(globe^2)^2
  T3 <- sum(globe) * sum(globe^3)
  T4 <- (sum(globe)^2) * sum(globe^2)
  T5 <- sum(globe)^4
  M1 <- T1/G
  M2 <- (T2 - T1)/(G^2 - G)
  M3 <- (T3 - T1)/(G^2 - G)
  M4 <- (T4 - 2 * T3 - T2 + 2 * T1)/(G * (G - 1) * (G - 2))
  M5 <- (T5 - 6 * T4 + 8 * T3 + 3 * T2 - 6 * T1)/(G * (G - 1) * (G - 2) * (G - 3))
  E.globe <- c(M1, M2, M3, M4, M5)
  return(E.globe)
}

fact <- function(G, m) sqrt(((G - m)/(G - 1))/m)

## sigma.fun <- function(m, E, esig2) {
##   vector <- c(m + 1/m - 2, (m - 1) * (m + 3/m - 2), -4 * (m - 1) * (1 - 1/m), 
##               -2 * (m - 1) * (m - 2) * (1 - 3/m), (m - 1) * (m - 2) * (m - 3)/m)
##   sd.x <- sqrt((sum(vector * E)/((m - 1)^2) - esig2^2))
##   sd.x
## }

# This function replaces the one above. Notice this has one extra argument
sigma.fun <- function(m, E, G){
  Beta1 <- c(1, -3, -4, 12, -6)
  Beta2 <- c(0, 1, 0, -2, 1)
  var.x <- (1/m - 1/G)*Beta1 %*% E + (2/(m - 1) - 2/(G - 1))*Beta2 %*% E
  sd.x <- sqrt(var.x)
  sd.x
}

## Local Universe ##
local.max <- function(setscores){
#  if(!as.logical(grep("local.zscore",colnames(setscores))))
#    stop("setscores must be from universe='local'")
  do.call(rbind,by(setscores,setscores[,1],
     function(x) x[which.max(x$local.zscore),]))
}
