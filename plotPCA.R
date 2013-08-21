
# Usage: Rscript infile.covar component1 component2 outfile.eps npop nind1 nind2 ... nindNPOP

library(ggplot2)

# Read parameters
args <- commandArgs(trailingOnly = TRUE);
infile <- args[1];
comp1 <- as.integer(args[2]);
comp2 <- as.integer(args[3]);
outfile <- args[4];
npop <- as.integer(args[5]);
nind <- as.array(args[6:(npop+5)]);
rm(args);

# Read input file
covar <- read.table(infile, stringsAsFact=F);
eig <- eigen(covar, symm=TRUE);

# Eigenvalues
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
id <- c(); for (i in 1:npop) id <- c(id, rep(paste("Pop",i), nind[i]));
PC$Pop <- id

title <- paste("PC",comp1," (",signif(eig$val[comp1], digits=3)*100,"%)"," / PC",comp2," (",signif(eig$val[comp2], digits=3)*100,"%)",sep="",collapse="")

ggplot() + ggtitle(title) + geom_point(data=PC, aes_string(x=paste("PC",comp1,sep=""), y=paste("PC",comp2,sep=""), color="Pop"))
ggsave(outfile)
unlink("Rplots.pdf", force=TRUE)
