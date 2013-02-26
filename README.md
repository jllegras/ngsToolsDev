Set of programs/scripts to analyze NGS data for multi-population genetics analyses.

### INSTALL

See INSTALL file on how to link to repository and compile all programs. Note that you need GSL library.

### INPUT FILES

These programs receive in input files produced by the software ANGSD.

A typical pipeline can be the following, assuming we have genotype likelihoods data for one pop in 'sim' format (e.g. from nsgSim). We assume 40 individuals.

This computes genotype posterior probabilities (.geno) as well as estimates of minor allele frequencies (.maf):
```angsd0.204/angsd.g++ -sim1 pop.glf.gz -nInd 40 -doGeno 64 -doMaf 2 -outfiles all -doMajorMinor 1 -doGlf 3```

This computes sample allele frequency posterior probabilities assuming no prior (.sfs):
```angsd0.204/angsd.g++ -outfiles pop -realSFS 1 -sim1 pop.glf.gz -nInd 40```

This estimates the overall SFS (.sfs.ml) using an ML approach:
```angsd0.204/misc/optimSFS.gcc -binput pop.sfs -nChr 80 -outnames pop.sfs```

This compute sample allele frequency posterior probabilities (.sfs.ml.norm):
```angsd0.204/misc/sfstools.g++ -nChr 80 -priorFile pop.sfs.ml -sfsFile pop.sfs -dumpBinary 1 > pop.sfs.ml.norm```


### ngsFST

Program to estimate FST from NGS data. It computes expected genetic variance components and estimate FST from those.
In input it receives posterior probabilities of sample allele frequencies (from angsd0.204 and sfstools) for each population. It may receive also a 2D-SFS as a prior and in this case in gets in input posterior probabilities with uniform prior (angsd with -realSFS 1 only). Additionally it can use a corrected product of marginal spectra as prior. In this case it receives in input posterior probabilities of sample allele frequencies (from angsd0.204 and sfstools) and also marginal SFS.

Output is a tab-separated text file. Each row is a site. Columns are: EA, EAB, FACT, (EA/EAB)+FACT, pvar; where EA is the expectation of genetic variance between populations, EAB is the expectation of the total genetic variance, FACT is the correcting factor for the ratio of expectations, (EA/EAB)+FACT is the per-site FST value, pvar is the probability for the site of being variable.

Run with no arguments for help. Please note that populations must have the exact same number of sites.

Examples:
```ngsTools/bin/ngsFST -postfiles pop1.sfs pop2.sfs -priorfile spectrum.txt -nind 20 20 -nsites 100000 -block_size 20000 -outfile pops.fst``` # using a 2D-SFS as prior, estimated using ngs2dSFS
```ngsTools/bin/ngsFST -postfiles pop1.sfs.ml.norm pop2.sfs.ml.norm -nind 20 20 -nsites 100000 -block_size 20000 -outfile pops.first.fst``` # here we don't provide prior files, so we do not correct for non-independece, but we use the output file as a first guess for the fst
```Rscript --vanilla --slave -e 'source("ngsTools/getMultiFST.R"); getMultiFST(filein="pops.first.fst", fileout="pops.global.fst", from_known=FALSE)'``` # this will compute a global fst and used it as a first guess for all sites; .R script is easily changeable for other purposes (e.g. same fst for regions fo 10-20-50Kbp)
```ngsTools/bin/ngsFST -postfiles pop1.sfs.ml.norm pop2.sfs.ml.norm -priorfiles pop1.sfs.ml pop2.sfs.ml -nind 20 20 -nsites 100000 -outfile pops.corrected.fst -fstfile fst.global.fst -K 0```

Parameters:

-postfiles: .sfs files with posterior probabilities of sample allele frequencies for each population

-fstfile: file with first guesses of FST for each site obtained running the program once, check getMultiFST.R

-priorfile: 2D-SFS to be used as a prior; you can use ngs2DSFS with parameter -relative set to 1

-outfile: name of the output file

-nind: number of individuals for each population

-nsites: total number of sites; in case you want to analyze a subset of sites this is the upper limit

-K: if set to 0: automatic setting of weighting function, otherwise lambda=1/(K*FST)

-verbose: level of verbosity, if 0 suppress all messages

-block_size: to be memory efficient, set this number as the number of sites you want to analyze at each chunk

-nsums: how many terms for the correction of the expectation of the ratio; set to 1

-firstbase: in case you want to analyze a subset of your sites this is the lower limit

-isfold: boolean, is your data folded or not?

### ngsCovar

Program to compute the expected sample covariance matrix from genotype posterior probabilities. It receives in input genotype posterior probabilities (from angsd2.04 -doGeno 64). It can receive in input also posterior probabilities of sample allele frequencies (from angsd0.204 and sfstools) for computing the probability of each site to be variant.

Run with no arguments for help.

Examples:
```ngsTools/bin/ngsCovar -probfile pop.geno -outfile pop.covar -nind 40 -nsites 100000 -block_size 20000 -call 0 -sfsfile pop.sfs.ml.norm```
```ngsTools/bin/ngsCovar -probfile pop.geno -outfile pop.covar -nind 40 -nsites 100000 -block_size 20000 -call 1```
```ngsTools/bin/ngsCovar -probfile pop.geno -outfile pop.covar -nind 40 -nsites 100000 -block_size 20000 -call 1 -minmaf 0.05```

Parameters:

-probfile: file with genotype posterior probabilities

-sfsfile: file with sample allele frequency posterior probabilities, to be used to compute probability of sites of being variable

-nind: how many individuals

-nsites: how many sites are in the data or the upper limit in case you want to analyze a subset

-outfile: name of the output file

-norm: if 1 matrix is normalized by sqrt(p(1-p)) as in Patterson et al 2006

-minmaf: ignore sites below this threhsold of minor allele frequency

-block_size: memory efficiency, number of sites for each chunk

-offset: lower limit of sites in case you want to analyze a subset

-call: call genotypes based on the maximum posterior probability

-verbose: level of verbosity

### nsgSim

Program to simulate NGS data for up to 3 populations setting an inbreeding coefficient or a FST. It outputs true genotypes, reads and genotype likelihoods.

Run with no arguments for help.

Example:
```ngsTools/bin/ngsSim -outfiles pop -npop 2 -nind 20 20 -nsites 100000 -depth 4 -pvar 0.10 -F 0.3 0.3```

Parameters:

-outfiles: prefix for output files

-npop: number of populations

-nind: number of individuals for each population

-nsites: number of sites

-errate: sequencing error rate

-depth: mean sequencing depth

-pvar: probability to each site is variable in the population

-mfreq: minimum population frequency

-F: FST value(s) in case of 2/3 populations, inbreeding coefficient in case of 1 population

-model: set to for 0=fixed errate or 1=variable errate

-simpleRand: boolean, set to 1 for quick random number generator

-seed: random number

-base_freq: background allele frequencies for A,C,G,T [0.25 0.25 0.25 0.25]

-multi_depth: Simulate uneven covered individuals. -multi_depth 6 10: first 10 individuals have 6X while the rest is as -depth

### ngs2dSFS

Program to estimate 2D-SFS from posterior probabilities of sample allele frequencies (from angsd0.505 and sfstools).

Run with no arguments for help. Please note that populations must have the exact same number of sites.

Example:
```ngsTools/bin/ngs2dSFS -postfiles pop.sfs.ml.norm pop.sfs.ml.norm -outfile spectrum.txt -relative 1 -nind 20 20 -nsites 100000 -block_size 20000```

Parameters:

-postfiles: file with sample allele frequency posterior probabilities for each population

-outfile: name of output file

-nind: number of individuals per population

-nsites: number of sites, or upper limit in case of analyzing a subset

-block_size: memory efficiency, number of sites for each chunk

-offset: lower limit in case of analyzing a subset

-maxlike: if 1 compute the MLE of joint allele frequency and sum across sites, if 0 it computes the sum of the products of likelihoods

-relative: boolean, if 1 number are relative frequencies from 0 to 1 which sum up 1; if 0 numbers are absolute counts of sites having a specific joint allele frequency

### Misc utilities

getMultiFST.R

This script converts the output of ngsFST and compute multiple-site FST and rewrite the file with this new values of FST (to be used as -firstfile)

ngsCovar.R

This script plots some PCA figures from the output of ngsCovar.



