
Set of programs/scripts to analyze NGS data for multi-population genetics analyses.

### INSTALL

After you linked the git repository enter the folder and type 'make'. You need GSL library.

### INPUT FILES

These programs receive in input files produces by the software ANGSD.

A typical pipeline can be the following, assuming we have genotype likelihoods data for one pop in 'sim' format(e.g. from nsgSim). We assume 40 individuals.

This computes genotype posterior probabilities (.geno) as well as estimates of minor allele frequencies (.maf):
```angsd0.204/angsd.g++ -sim1 pop.glf.gz -nInd 40 -doGeno 64 -doMaf 2 -outfiles all -doMajorMinor 1 -doGlf 3```

This computes sample allele frequency posterior probabilities assuming no prior (.sfs):
angsd0.204/angsd.g++ -outfiles pop -realSFS 1 -sim1 pop.glf.gz -nInd 40

This estimates the overall SFS (.sfs.ml) using an ML approach:
angsd0.204/misc/optimSFS.gcc -binput pop.sfs -nChr 80 -outnames pop.sfs

This compute sample allele frequency posterior probabilities (.sfs.ml.norm):
angsd0.204/misc/sfstools.g++ -nChr 80 -priorFile pop.sfs.ml -sfsFile pop.sfs -dumpBinary 1 > pop.sfs.ml.norm


### ngsFST

Program to estimate FST from NGS data. It computes expected genetic variance components and estimate FST from those.
In input it receives posterior probabilities of sample allele frequencies (from angsd0.204 and sfstools) for each population. It may receive also a 2D-SFS as a prior and in this case in gets in input posterior probabilities with uniform prior (angsd with -realSFS 1 only). Additionally it can use a corrected product of marginal spectra as prior. In this case it receives in input posterior probabilities of sample allele frequencies (from angsd0.204 and sfstools) and also marginal SFS.

Run with no arguments for help.

Examples:
ngsTools/bin/ngsFST -postfiles pop1.sfs pop2.sfs -priorfile spectrum.txt -nind 20 20 -nsites 100000 -block_size 20000 -outfile pops.fst # using a 2D-SFS as prior, estimated using ngs2dSFS
ngsTools/bin/ngsFST -postfiles pop1.sfs.ml.norm pop2.sfs.ml.norm -nind 20 20 -nsites 100000 -block_size 20000 -outfile pops.first.fst # here we don't provide prior files, so we do not correct for non-independece, but we use the output file as a first guess for the fst
Rscript --vanilla --slave -e 'source("ngsTools/getMultiFST.R"); getMultiFST(filein="pops.first.fst", fileout="pops.global.fst", from_known=FALSE)' # this will compute a global fst and used it as a first guess for all sites; .R script is easily changeable for other purposes (e.g. same fst for regions fo 10-20-50Kbp)
ngsTools/bin/ngsFST -postfiles pop1.sfs.ml.norm pop2.sfs.ml.norm -priorfiles pop1.sfs.ml pop2.sfs.ml -nind 20 20 -nsites 100000 -outfile pops.corrected.fst -fstfile fst.global.fst -K 0


### ngsCovar

Program to compute the expected sample covariance matrix from genotype posterior probabilities. It receives in input genotype posterior probabilities (from angsd2.04 -doGeno 64). It can receive in input also posterior probabilities of sample allele frequencies (from angsd0.204 and sfstools) for computing the probability of each site to be variant.

Run with no arguments for help.

Examples:
ngsTools/bin/ngsCovar -probfile pop.geno -outfile pop.covar -nind 40 -nsites 100000 -block_size 20000 -call 0 -sfsfile pop.sfs.ml.norm
ngsTools/bin/ngsCovar -probfile pop.geno -outfile pop.covar -nind 40 -nsites 100000 -block_size 20000 -call 1
ngsTools/bin/ngsCovar -probfile pop.geno -outfile pop.covar -nind 40 -nsites 100000 -block_size 20000 -call 1 -minmaf 0.05

### nsgSim

Program to simulate NGS data for up to 3 populations setting an inbreeding coefficient or a FST. It outputs true genotypes, reads and genotype likelihoods.

Run with no arguments for help.

Example:
ngsTools/bin/ngsSim -outfiles pop -npop 2 -nind 20 20 -nsites 100000 -depth 4 -pvar 0.10 -F 0.3 0.3

### ngs2dSFS

Program to estimate 2D-SFS from posterior probabilities of sample allele frequencies (from angsd0.204 and sfstools).

Run with no arguments for help.

Example:
ngsTools/bin/get2DSFS -postfiles pop.sfs.ml.norm pop.sfs.ml.norm -outfile spectrum.txt -relative 1 -nind 20 20 -nsites 100000 -block_size 20000

### Misc utilities

getMultiFST.R

This script converts the output of ngsFST and compute multiple-site FST and rewrite the file with this new values of FST (to be used as -firstfile)

ngsCovar.R

This script plots some PCA figures from the output of ngsCovar.



