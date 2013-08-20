
# Usage: Rscript infile.stat outfile.eps npop name_pop1 name_pop2

# Read parameters
args <- commandArgs(trailingOnly = TRUE);
infile <- args[1];
outfile <- args[2];
npop <- as.integer(args[3]);
nampop <- as.array(args[4:(npop+3)]);
rm(args);

# Read input file
values <- read.table(infile, stringsAsFact=F);
pos=values[,1]+(values[,2]-values[,1])/2;

# Plot
if (npop==1)
{
	x11();
	par(mfrow=c(2,1));
	plot(x=pos, y=values[,3], col="black", xlab="Position", ylab="Segregating sites", ty="b", main=nampop[1], sub="");
	plot(x=pos, y=values[,4], col="black", xlab="Position", ylab="Expected heterozygosity", ty="b", main="", sub="");	

} else {

        x11();
        par(mfrow=c(3,1));
        plot(x=pos, y=values[,3], col="red", xlab="Position", ylab="Segregating sites", ty="b", main="", sub="");
	lines(x=pos, y=values[,5], col="blue", ty="b");
	legend(legend=nampop, col=c("red", "blue"), x="topright", pch=1);
        plot(x=pos, y=values[,4], col="red", xlab="Position", ylab="Expected heterozygosity", ty="b", main="", sub="");
        lines(x=pos, y=values[,6], col="blue", ty="b");
        legend(legend=nampop, col=c("red", "blue"), x="topright", pch=1);
        plot(x=pos, y=values[,7], col="black", xlab="Position", ylab="Fixed differences", ty="b", main="", sub="");       

}

dev.copy2eps(file=outfile);
dev.off();

