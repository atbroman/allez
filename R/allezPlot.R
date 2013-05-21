
## Order design matrix ##
## 1. Column order by GO/KEGG sum
## 2. Row order by GO/KEGG category, gene reduction value
ordMat <- function(aMat,allez.out){
  rind <- cind <- character(0)
  for(i in 1:ncol(aMat)){
    mat <- if(i==1) aMat*allez.out$aux$globe else
           aMat[-match(rind,rownames(aMat)),
                -match(cind,colnames(aMat)),drop=FALSE]*
           allez.out$aux$globe[-match(rind,names(allez.out$aux$globe))]
    smax <- which.max(apply(mat,2,sum)) ## set with highest sum
    cind <- c(cind,colnames(mat)[smax])
    rord <- order(mat[,smax],decreasing=TRUE)
    rind <- c(rind,rownames(mat)[rord][mat[rord,smax]>0])
  }
  apply(aMat[rind,cind,drop=FALSE],2,"*",allez.out$aux$globe[rind])
}

allezplot <- function(aOrd, allez.out,
              glab=c("none","gene_id","symbol"), ...){
  require(GO.db)
  require(KEGG.db)

  goterm <- toTable(GOTERM)[,c("go_id","Term")]
  names(goterm) <- c("id","term")
  keggterm <- toTable(KEGGPATHID2NAME)
  names(keggterm) <- c("id","term")
  allterm <- rbind(goterm,keggterm)
  term <- allterm[match(colnames(aOrd),allterm$id),]

  glab <- match.arg(glab)

  p <- par(no.readonly=TRUE)

  xlabs <- switch(glab,
          "none" = 1:nrow(aOrd),
          "gene_id" = rownames(aOrd),
          "symbol" = allez.out$aux$set.data[
            match(rownames(aOrd),allez.out$aux$set.data[,2]),"symbol"])

  ## Space for text in xlim ##
  xpos <- apply(aOrd,2,function(x)
       (1:length(x))[x>0][which.max((1:length(x))[x>0])])
  mwidth1 <- strwidth(xlabs,units="inches")
  mwidth2 <- strwidth(term$id,units="inches")
  twidth <- strwidth(term$term,units="inches")+0.5*par("cin")[1]

  if(glab=="none"){
    par(mai=c(par("mai")[1],max(mwidth2)+par("mgp")[2]*par("cin")[2],
          par("mai")[3:4]))
    image(1:nrow(aOrd),1:ncol(aOrd),aOrd[,ncol(aOrd):1],
          xlab="",ylab="", yaxt="n",col=gray(seq(1,0,length=64)),
          xlim=c(0.5,max((xpos+0.5)/(1-twidth/par("pin")[1]))))
  } else {
    pwidth <- (xpos+1)*par("cin")[2]+twidth ## plot width, inches ##
    par(pin=c(max(pwidth),par("pin")[2]),
        mai=c(max(mwidth1)+par("mgp")[2]*par("cin")[2],
          max(mwidth2)+par("mgp")[2]*par("cin")[2],par("mai")[3:4]))
    if(max(pwidth)/diff(par("plt")[1:2])>par("fin")[1])
      warning(paste("Increase width of figure to at least",
              round(max(pwidth)/diff(par("plt")[1:2]),2),
              "inches and re-run allezPlot"))
    image(1:nrow(aOrd),1:ncol(aOrd),aOrd[,ncol(aOrd):1],
          xlab="",ylab="", yaxt="n", xaxt="n",col=gray(seq(1,0,length=64)),
          xlim=c(0.5,max((xpos+0.5)/(1-twidth/par("pin")[1]))))
    axis(side=1,at=1:nrow(aOrd),labels=xlabs,las=2, ...)
  }
  axis(side=2,at=1:ncol(aOrd),labels=colnames(aOrd)[ncol(aOrd):1],las=1, ...)
  text(xpos+0.5,ncol(aOrd):1,term$term,pos=4, ...)
  par(p)
}

allezPlot <- function(allez.out,
                     n.low=5,
                     n.upp=500,
                     zthr=3,
                     gmax=20,
                     glab=c("none","gene_id","symbol"),
                     ...){
aMat <- allezMat(allez.out,n.low,n.upp,zthr)
aOrd <- ordMat(aMat,allez.out)
allezplot(aOrd,allez.out,
          glab=ifelse(nrow(aOrd)<=gmax,glab,"none"),
          ...)
}
