
# Usage: Rscript infile.covar component1 component2 outfile.eps npop nind1 nind2 ... nindNPOP

# Read parameters
args <- commandArgs(trailingOnly = TRUE);
infile <- args[1];
comp1 <- as.integer(args[2]);
comp2 <- as.integer(args[3]);
outfile <- args[4];
npop <- as.integer(args[5]);
nind <- as.array(args[6:(npop+5)]);
rm(args);
cols <- c(); for (i in 1:npop) cols <- c(cols, rep(rainbow(npop)[i], nind[i]));

# Read input file
covar <- read.table(infile, stringsAsFact=F);
eig <- eigen(covar, symm=TRUE);

# Eigenvalues
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");

# Plot
x11();
plot(eig$vec[,comp1], eig$vec[,comp2], col=cols, xlab=paste("PC",comp1," (",signif(eig$val[comp1], digits=3)*100,"%)",sep="",collapse=""), ylab=paste("PC",comp2," (",signif(eig$val[comp2], digits=3)*100,"%)",sep="",collapse=""), main="PCA", sub=""); 
dev.copy2eps(file=outfile);
dev.off();


