
# Usage: Rscript infile.fst outfile.eps win step

# Read parameters
args <- commandArgs(trailingOnly = TRUE);
infile <- args[1];
outfile <- args[2];
win <- as.integer(args[3]);
step <- as.integer(args[4]);
rm(args);

# Read input file
values <- read.table(infile, stringsAsFact=F);
cat("Overall FST:",sum(values[,1])/sum(values[,2]),"\n");

# Windows
len=nrow(values);
start=seq(1, len-win, step);
end=start+win; # it requires that all windows have the same nr of sites, if less the window (usually at the end) is discarded
pos=start+(win/2);
fst=c(); for (i in 1:length(start)) fst[i]=sum(values[start[i]:end[i],1])/sum(values[start[i]:end[i],2]);

# Plot
x11();
plot(x=pos, y=fst, col="black", xlab="Position", ylab=expression(F[ST]), ty="b", main="", sub="", frame.plot=F);
dev.copy2eps(file=outfile);
dev.off();

