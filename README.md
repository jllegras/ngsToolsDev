Set of programs/scripts to analyze NGS data for multi-population genetic analysis.

### INSTALL

See INSTALL file on how to link to repository and compile all programs. Note that you need GSL library (actually this requirement will be soon removed).

### INPUT FILES

The programs receive in input files produced by software ANGSD.

A typical pipeline can be the following. Assuming we have genotype likelihoods data for one pop in 'sim1' format (e.g. from nsgSim). We assume 40 individuals.
We assume to use a version of ANGSD 0.505 or higher. Please check ANGSD web site for other accepted genotype likelihood formats.

First we compute genotype posterior probabilities (.geno) as well as estimates of minor allele frequencies (.maf):
``angsd -sim1 pop.glf.gz -nInd 40 -doGeno 32 -doPost 1 -doMaf 2 -out pops.geno -doMajorMinor 1``

We then compute sample allele frequency posterior probabilities assuming no prior (.sfs, values will be in log format):
``angsd -sim1 pop.glf.gz -nInd 40 -realSFS 1 -out pops``

and we estimate the overall SFS (.sfs.ml) using an ML approach (Nielsen et al. 2012):
``misc/optimSFS.gcc -binput pops.sfs -nChr 80 -nThreads 10``

Finally we compute sample allele frequency posterior probabilities using the SFS as a prior (.sfs.ml.norm, values will not be in log format anymore):
``misc/sfstools.g++ -sfsFile pops.sfs -nChr 80 -priorFile pops.sfs.ml -dumpBinary 1 > pops.norm``

Please note that if your data is folded you should use option -fold 1 at step -realSFS 1 and the set -nChr equal to -nInd.

### ngsFST

Program to estimate FST from NGS data. It computes expected genetic variance components and estimate FST from those.
In input it receives posterior probabilities of sample allele frequencies for each population (ANGSD + sfstools). It may receive also a 2D-SFS as a prior and in this case it gets in input posterior probabilities with uniform prior (ANGSD with -realSFS 1 only, do not run sfstools and set -islog 1). You can give also 2 marginal spectra as prior. 

Output is a tab-separated text file. Each row is a site. Columns are: EA, EAB, FACT, (EA/EAB)+FACT, pvar; where EA is the expectation of genetic variance between populations, EAB is the expectation of the total genetic variance, FACT is the correcting factor for the ratio of expectations, (EA/EAB)+FACT is the per-site FST value, pvar is the probability for the site of being variable.

Run with no arguments for help. Please note that populations must have the exact same, and corresponding, number of sites.

Examples:
``ngsTools/bin/ngsFST -postfiles pop1.sfs pop2.sfs -priorfile spectrum2D.txt -nind 20 20 -nsites 100000 -block_size 20000 -outfile pops.fst -islog 1`` # using a 2D-SFS as a prior, estimated using ngs2dSFS
``ngsTools/bin/ngsFST -postfiles pop1.sfs pop2.sfs -priorfiles spectrum1.txt spectrum2.txt -nind 20 20 -nsites 100000 -block_size 20000 -outfile pops.fst -islog 1`` # using marginal spectra as priors, estimated using optimSFS
``ngsTools/bin/ngsFST -postfiles pop1.sfs.ml.norm pop2.sfs.ml.norm -nind 20 20 -nsites 100000 -block_size 20000 -outfile pops.fst -islog 1`` # here we don't provide prior files, so we directly provide posterior probabilities (ANGSD+sfstools), and therefore we do not correct for non-independence;

Parameters:

-postfiles: .sfs files with posterior probabilities of sample allele frequencies for each population (with or without running sfstools)

-priorfile: 2D-SFS to be used as a prior; you can use ngs2DSFS with parameter -relative set to 1

-priorfile2: marginal spectra to be used as a prior; you can use optimSFS in ANGSD

-outfile: name of the output file

-nind: number of individuals for each population

-nsites: total number of sites; in case you want to analyze a subset of sites this is the upper limit

-verbose: level of verbosity, if 0 suppress all messages

-block_size: to be memory efficient, set this number as the number of sites you want to analyze at each chunk

-firstbase: in case you want to analyze a subset of your sites this is the lower limit

-isfold: boolean, is your data folded or not?

-islog: boolean, are postfiles in -log (from -realSFS 1 only, required if 2D-SFS is given)? If you use sfstools then -islog 1

### ngsCovar

Program to compute the expected correlation matrix between individuals from genotype posterior probabilities. It receives in input genotype posterior probabilities (from angsd -doGeno 32). It can receive in input also posterior probabilities of sample allele frequencies (from angsd and sfstools) for computing the probability of each site to be variant.

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

-norm: if 0 no normalization, if 1 matrix is normalized by (p(1-p)) as in Patterson et al 2006, if 2 normalization is 2p(1-p)

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

This script converts the output of ngsFST and compute multiple-site FST and rewrite the file with this new values of FST

ngsCovar.R

This script plots some PCA figures from the output of ngsCovar.

GetSubSFS

This program extracts a subset of .sfs files.

GetSubSim

This program extracts a subset of individuals or sites after running ngsSim.
