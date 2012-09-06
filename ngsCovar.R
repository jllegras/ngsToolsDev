library(optparse)
library(ggplot2)
library(reshape2)
library(scatterplot3d)
#library(RMTstat)

option_list <- list(make_option(c('-i', '--in_file'), action='store', type='character', default=NULL, help='Input distrib file'),
                    make_option(c('-a', '--in_annot'), action='store', type='character', default=NULL, help='File with IDs (one per line)'),
                    make_option(c('-m', '--mark_annot'), action='store', type='character', default=NULL, help='File with IDs to mark (requires --in_annot)'),
                    make_option(c('-r', '--rem_annot'), action='store', type='character', default=NULL, help='File with IDs to remove (requires --in_annot).'),
                    make_option(c('-n', '--no_names'), action='store_true', type='logical', default=FALSE, help='Do not plot names (requires --in_annot).'),
                    make_option(c('-s', '--n_sites'), action='store', type='integer', default=0, help='Num sites used (for tw test).')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

############################################################################################

#opt <- list(in_file = '../allsp-LC_chrALL.geno_covar', in_annot='/home/fgvieira/run/UCB/Rice/real_data/allsp-LC.annot', no_names=FALSE, n_sites=125556657)
#opt <- list(in_file = '/home/fgvieira/subset.50Msites.MINMAF_05.geno_covar', in_annot='/home/fgvieira/subset.50Msites.annot', no_names=TRUE, mark_annot='/home/fgvieira/subset.50Msites.annot_MARK', n_sites=0)

twtest <- function(eigen, nsites) {
  n <- nsites;
  values <- eigen$values;
  values[which(values<0)] <- 0;
  m <- length(values);
  
  mu <- ( sqrt(n-1) + sqrt(m) )^2 / n;
  sig <- ( ( sqrt(n-1) + sqrt(m) ) / n ) * ( (1/sqrt(n-1)) + (1/sqrt(m)) )^(1/3);
  x <- (values - mu) / sig;
  
  pv = rep(NA, length(values));

  return(1-ptw(x));
}

############################################################################################

colors <- c("blue","red","green","cyan","black","pink","yellow")
group_col <- NULL

if(is.null(opt$in_annot) && (! is.null(opt$mark_annot) || ! is.null(opt$mark_annot)) ) {
  cat("ERROR: \"mark_annot\" and \"rem_annot\" require \"in_annot\"!", fill=TRUE)
  exit()
}



### Read data
cat('# Reading input file...', fill=TRUE)
covmat <- read.table(opt$in_file, stringsAsFac=FALSE)

annot <- data.frame(sample("", ncol(covmat), replace=TRUE),
                    sample("", ncol(covmat), replace=TRUE),
                    sample("black", ncol(covmat), replace=TRUE))

### Read annots
cat('# Reading annotations file...', fill=TRUE)
if( !is.null(opt$in_annot) ) {
  annot <- read.table(opt$in_annot, stringsAsFac=FALSE, sep="\t")
  group <- unique(annot[,1])
  for (i in 1:length(group)) {
    annot[annot[,1] == group[i],3] = colors[i]
  }

  rownames(covmat) = colnames(covmat) = annot[,"V2"]
}


### Read REM annots
if( ! is.null(opt$rem_annot) ) {
  cat('# Removing annotations...', fill=TRUE)
  rem_annot <- read.table(opt$rem_annot, stringsAsFac=FALSE)[,1]
  good_annot <- colnames(covmat)[!colnames(covmat) %in% rem_annot]

  annot <- annot[annot[,"V2"] %in% good_annot,]
  covmat <- covmat[good_annot, good_annot]
}



### Read MARK annots
if( ! is.null(opt$mark_annot) ) {
  cat('# Marking annotations...', fill=TRUE)
  mark_annot <- read.table(opt$mark_annot, stringsAsFac=FALSE)[,1]
}



### PCA
cat('# Calculating PCA...', fill=TRUE)
eigen = eigen(covmat, symm=TRUE)
if(opt$n_sites != 0){
  twtest(eigen, opt$n_sites)
}
eigenV <- eigen$vectors

cat('# Plotting data...', fill=TRUE)
pdf(paste(opt$in_file,".pdf",sep=""))
### Density plot
plot(density(as.matrix(covmat)))

### R plots
#plot(eigenV[,1], eigenV[,2], col=annot[,3], xlab="PC1", ylab="PC2", pch=20)
#text(eigenV[,1], eigenV[,2], labels=annot[,2], pos=3)
#legend("topright", group, col=colors, pch=20)
#plot(eigenV[,2], eigenV[,3], col=annot[,3], xlab="PC2", ylab="PC3", pch=20)
#text(eigenV[,2], eigenV[,3], labels=annot[,2], pos=3)
#legend("topright", group, col=colors, pch=20)
#plot(eigenV[,3], eigenV[,4], col=annot[,3], xlab="PC2", ylab="PC3", pch=20)
#text(eigenV[,3], eigenV[,4], labels=annot[,2], pos=3)
#legend("topright", group, col=colors, pch=20)

### GGPLOT2
df <- as.data.frame(cbind(eigenV[,c(1:6)], annot))
colnames(df) <- c("PC1","PC2","PC3","PC4","PC5","PC6","Species","ID","color")

x <- matrix(c("PC1","PC2",
              "PC2","PC3",
              "PC3","PC4"), ncol=2, byrow=TRUE)

for (i in 1:3) {
  plot <- ggplot(df, aes_string(x=x[i,1], y=x[i,2], label="ID")) + geom_point(aes(colour=Species))

  if( ! is.null(opt$mark_annot) ) {
    plot <- plot + geom_point(data=df[df$ID %in% mark_annot, ], size=5, colour="green", alpha=0.3)
  }

  if(! opt$no_names) {
    plot <- plot + geom_text(hjust=0.5, vjust=-0.5)
  }
  print(plot)
}


### Scatter 3D
scatterplot3d(eigenV[,1],eigenV[,2],eigenV[,3], type="h", color=annot[,3], xlab="PC1", ylab="PC2", zlab="PC3", pch=20)
scatterplot3d(eigenV[,3],eigenV[,2],eigenV[,1], type="h", color=annot[,3], xlab="PC3", ylab="PC2", zlab="PC1", pch=20)

dev.off()
quit()

library(rgl)
plot3d(eigenV[,1],eigenV[,2],eigenV[,3], main="3D Scatterplot", col=group_col, xlab="PC1", ylab="PC2", zlab="PC3", pch=19, size=3)

library(Rcmdr)
scatter3d(eigenV[,1],eigenV[,2],eigenV[,3], main="3D Scatterplot", xlab="PC1", ylab="PC2", zlab="PC3")
