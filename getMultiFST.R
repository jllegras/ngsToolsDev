getMultiFST<-function(filein, fileout, from_known=FALSE) {

  # this script converts the output of ngsFST and compute multiple-site FST and rewrite the file with this new values of FST (to be used as -firstfile) 

  val=c()
  if (from_known) { # this is for compatibility with gimmeFST
    rt=read.table(filein,head=T)
    val=sum(rt[,8])/sum(rt[,9])
  } else {
    rt=read.table(filein,head=F)
    val=sum(rt[,1])/sum(rt[,2])
  }

  write.table(cbind(rep(0,nrow(rt)),rep(0,nrow(rt)),rep(0,nrow(rt)),rep(val,nrow(rt)),rep(0,nrow(rt))), sep="\t", quote=F, row.names=F, col.names=F, file=fileout)

}
