
## Order design matrix ##
## 1. Column order by GO/KEGG sum
## 2. Row order by GO/KEGG category, gene reduction value
ordMat <- function(aMat,allez.out){
  rind <- cind <- character(0)
  for(i in 1:ncol(aMat)){
    mat <- if(i==1) apply(aMat,2,"*",allez.out$aux$globe) else
           apply(aMat[-match(rind,rownames(aMat)),
                -match(cind,colnames(aMat)),drop=FALSE],2,"*",
           allez.out$aux$globe[-match(rind,names(allez.out$aux$globe))])
    smax <- which.max(apply(mat,2,sum)) ## set with highest sum
    cind <- c(cind,colnames(mat)[smax])
    rord <- order(mat[,smax],decreasing=TRUE)
    rind <- c(rind,rownames(mat)[rord][mat[rord,smax]>0])
  }
  apply(aMat[rind,cind,drop=FALSE],2,"*",allez.out$aux$globe[rind])
}

ordMat2 <- function(aMat,allez.out){
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

allezplot <- function(aOrd,...){
  require(GO.db)
  require(KEGG.db)

  goterm <- toTable(GOTERM)[,c("go_id","Term")]
  names(goterm) <- c("id","term")
  keggterm <- toTable(KEGGPATHID2NAME)
  names(keggterm) <- c("id","term")
  allterm <- rbind(goterm,keggterm)
  term <- allterm[match(colnames(aOrd),allterm$id),]

  p <- par(no.readonly=TRUE)
  ## Space for text in xlim ##
  xpos <- apply(aOrd,2,function(x)
       (1:length(x))[x>0][which.max((1:length(x))[x>0])])

  mwidth <- strwidth(term$id,units="inches")
  mai <- par("mai")
  par(mai=c(mai[1],mai[4]+max(mwidth),mai[3:4]))

  twidth <- (strwidth(term$term,units="inches")+
             0.5*par("cin")[1])/par("pin")[1]            
  image(1:nrow(aOrd),1:ncol(aOrd),aOrd[,ncol(aOrd):1],
   xlab="",ylab="", yaxt="n",col=gray(seq(1,0,length=64)),
   xlim=c(0.5,max((xpos+0.5)/(1-twidth))),...)
  axis(side=2,at=1:ncol(aOrd),labels=colnames(aOrd)[ncol(aOrd):1],las=1,...)
  text(xpos+0.5,ncol(aOrd):1,term$term,pos=4,...)
  par(p)
}

allezPlot <- function(allez.out,
                     n.low=5,
                     n.upp=500,
                     zthr=3,
                     ...){
system.time(  aMat <- allezMat(allez.out,n.low,n.upp,zthr))
system.time(  aOrd <- ordMat(aMat,allez.out))
  allezplot(aOrd,...)
}
