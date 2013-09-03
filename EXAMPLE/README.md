
## EXAMPLES

These are several examples on how to generate and analyse NGS data using ngsTools.

# Configuration
Change the path to software according to your installation directory.

    NGSTOOLS=/home/mfumagalli/Documents/Software/ngsTools    
    ANGSD=/home/mfumagalli/Documents/Software/angsd0.543    
  
## Principal Component Analysis (PCA)

# Data simulation 

Simulation settings:
* 3 populations, where the closest populations (2nd and 3rd) differ from an average FST of 0.1 and these differ from the other one (1t) from an average FST of 0.3;
* 10, 8, 6 individuals at population 1, 2, and 3, respectively;
* 1,000 independet sites, all polymorphic in the population (for ease of computational time);
* 4X mean sequencing depth;
* 0.01 mean sequencing error rate, variable across sites;
* minor options: no inbreeding, minimum frequency in the population of 0.005, output files with prefix testA.

Under these conditions, we simulate NGS data with this command line:

    $NGSTOOLS/bin/ngsSim -npop 3 -nind 10 8 6 -nsites 1000 -errate 0.01 -depth 4 -pvar 1 -mfreq 0.005 -F 0.3 0.1 -model 1 -outfiles testA    

This command will generate several files for the whole population (testA), as well as for each population (testA1, testA2, testA3). Here is a description of each file:
* .geno: true genotypes; each line is a site and each tab-separated column is an individual; genotypes are coded as 0-3, where 0 is the ancestral state;
* .reads.txt: sequencing reads counts; first 24 rows (with 24 number of individuals) represent the counts of reads from base coded as 0 (1st column) to base coded as 3 (4th column), and so on for each site; therefore data has dimension 24,000x4;
* .frq: true Site Frequency Spectrum, tab-separated, from 0 to 48 frequencies of derived alleles;
* .args: arguments used in the simulation;
* .par: simulation settings at each site (on rows): drawn ancestral population allele frequency, inbreeding coefficient, first value of FST, second value of FST;
* .seq.gz: sequencing reads in ACGT format, with A always representing the ancestral state;
* .glf.gz: genotype likelihoods in binary format; values are in double format; first 24 rows show all possible 10 genotype likelihoods in case of biallelic sites, and so on for each site; the order of genotypes is: 0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8=GT, 9=TT;

For PCA we are only interested in data from the whole population (testA), so we will ignore all testA1, test2, and test3 file.

    rm testA1* testA2* testA3*    

# Posterior probabilities of genotypes and sample allele frequencies

We use ANGSD to compute genotype posterior probabilities. Please refer to http://popgen.dk/angsd for more details on ANGSD parameters. Here is a possible command line:

    $ANGSD/angsd -sim1 testA.glf.gz -nInd 24 -doGeno 32 -doPost 1 -doMaf 2 -out testA.geno -doMajorMinor 1    

which will generate these fies:
* testA.geno.arg: parameters used;
* testA.geno.mafs: estimates of minor allele frequencies;
* testA.geno.geno: genotype posterior probabilities in binary format; values are in double format; each row show likelihoods for 2 genotypes (e.g. AA AG GG) assuming a minor and a major allele;

In case we want to use weighting scheme on each site rather than calling SNPs, we need to calculate posterior probabilities of sample allele frequencies. These commands will generate these files in ANGSD:

    $ANGSD/angsd -sim1 testA.glf.gz -nInd 24 -realSFS 1 -out testA.rf    
    $ANGSD/misc/optimSFS -binput testA.rf.sfs -nChr 48 -nThreads 10    
    $ANGSD/misc/sfstools -sfsFile testA.rf.sfs -nChr 48 -priorFile testA.rf.sfs.ml -dumpBinary 1 > testA.rf.sfs.norm    

testA.rf.sfs.norm file is in binary format; values are in double format; in case of unfolded data, each row has 49 values (2N+1 with N individuals) representing the posterior probability of having a certain derived allele frequency; in case of folded data, each row will have N+1 values.

# Covariance matrix

We use ngsCovar to estimate a covariance matrix between pairs of individuals. This can be achieved in different ways. For low coverage sequencing depth, we recommend to use `-norm 0` option which disables normalization proposed in Patterson et al. PLoS Genetics (2006).
The first way is to compute an approximation of the posterior of the covariance matrix, by weighting each site by its probability of being variable, as proposed in Fumagalli et al. Genetics (2013). Here is to command:

     $NGSTOOLS/bin/ngsCovar -probfile testA.geno.geno -outfile testA.covar1 -nind 24 -nsites 1000 -call 0 -sfsfile testA.rf.sfs.norm -norm 0     

Alternatively, one can remove non-variable or low-frequency sites with the option `-minmaf`:
    
    $NGSTOOLS/bin/ngsCovar -probfile testA.geno.geno -outfile testA.covar2 -nind 24 -nsites 1000 -call 0 -minmaf 0.05    

Finally, in case of high sequencing depth, one can call genotypes as the genotype with the highest posterior probability: 

    $NGSTOOLS/bin/ngsCovar -probfile testA.geno.geno -outfile testA.covar3 -nind 24 -nsites 100 -call 1 -minmaf 0.05     

These commands will produce text files with a symmetric covariance matrix MxM, with M individuals.

# PCA plot

We can use a simple R script to perform an eigenvalue decomposition of the covariance matrix and plot the PCA results. First, let's create a dummy plink cluster file.

    Rscript --vanilla --slave -e 'write.table(cbind(seq(1,24),rep(1,24),c(rep("A",10),rep("B",8),rep("C",6))), row.names=F, col.names=c("FID","IID","CLUSTER"), sep="\t", file="testA.clst", quote=F)'    

Assuming we want to plot the first 2 PCA components from our previous results, the command would be:

    Rscript --vanilla --slave $NGSTOOLS/plotPCA.R -i testA.covar1 -c 1-2 -a testA.clst -o testA.pca.eps     

Please note that you need 'ggplot2' and 'optparse' R libraries installed. This script will output the explained genetic variance for each component and save as output the PCA plot. As a proof, PCA plot from called genotypes (testA.covar3) shows a less clear clustering of population. This can be seen by this command:

    Rscript --vanilla --slave $NGSTOOLS/plotPCA.R testA.covar3 1 2 testA.pca.call.eps 3 10 8 6     

Much more complex PCA plots can be obtained using the script plotPCAadv.R.

## SUMMARY STATISTICS

# Data simulation

Simulation settings:
* 2 populations, differentiated by an average FST of 0.2 (please note that in the arguments you need to pass this value twice);
* 10 and 8 individuals per population respectively;
* 10,000 independet sites, all polymorphic in the population (for ease of computational time);
* 4X mean sequencing depth;
* 0.01 mean sequencing error rate, variable across sites;
* minor options: no inbreeding, minimum frequency in the population of 0.005, output files with prefix testB.

Under these conditions, we simulate NGS data with this command line:

    $NGSTOOLS/bin/ngsSim -npop 2 -nind 10 8 -nsites 10000 -errate 0.01 -depth 4 -pvar 1 -mfreq 0.005 -F 0.2 0.2 -model 1 -outfiles testB    

    rm testB.args testB.frq testB.geno testB.glf.gz testB.par testB.reads.txt testB.seq.gz

Please refer to the previous example for a description of each file produced.

# Posterior probabilities of sample allele frequencies

We use ANGSD to compute sample allele frequency posterior probabilities for each population separately. Please refer to http://popgen.dk/angsd for more details on ANGSD parameters. Here is the command line:

    $ANGSD/angsd -sim1 testB1.glf.gz -nInd 10 -realSFS 1 -out testB1.rf
    $ANGSD/misc/optimSFS -binput testB1.rf.sfs -nChr 20 -nThreads 10
    $ANGSD/misc/sfstools -sfsFile testB1.rf.sfs -nChr 20 -priorFile testB1.rf.sfs.ml -dumpBinary 1 > testB1.rf.sfs.norm    

    $ANGSD/angsd -sim1 testB2.glf.gz -nInd 8 -realSFS 1 -out testB2.rf
    $ANGSD/misc/optimSFS -binput testB2.rf.sfs -nChr 16 -nThreads 10
    $ANGSD/misc/sfstools -sfsFile testB2.rf.sfs -nChr 16 -priorFile testB2.rf.sfs.ml -dumpBinary 1 > testB2.rf.sfs.norm    

testB1.rf.sfs.norm and testB2.rf.sfs.norm filee are in binary format; values are in double format; in case of unfolded data, each row has 49 values (2N+1 with N individuals) representing the posterior probability of having a certain derived allele frequency; in case of folded data, each row will have N+1 values.

# Summary statistics
 
We use ngsStat to compute expectations of some basic statistics of the data, specifically the number of segregating sites, the expected heterozygosity, and the number of fixed differences between populations. We also want to compute these quantities in non-overlapping sliding windows of 100 sites each, by using the options `-iswin` and `-block_size`. Here is the command to achieve this goal:

     $NGSTOOLS/bin/ngsStat -npop 2 -postfiles testB1.rf.sfs.norm testB2.rf.sfs.norm -nsites 10000 -iswin 1 -nind 10 8 -outfile testB.stat -isfold 0 -islog 0 -block_size 100    

This produces a file with these values, for each window: start, end, number of variable sites in pop 1, expected heterozygosity in pop 1, number of variable sites in pop 2, expected heterozygosity in pop 2, number of fixed differences between populations. Values can be plot by a simple R script, either for both populations or only one:

     Rscript --vanilla --slave $NGSTOOLS/plotSS.R testB.stat testB.stat.eps 2 pop1 pop2   
     Rscript --vanilla --slave $NGSTOOLS/plotSS.R testB.stat testB.stat.pop1.eps 1 pop1    

## FST

We use the same data simulated and file generated for the previous analysis. In case of low coverage sequencing depth, an improvement in the estimation accuracy can be achieved by using the joint-SFS as a prior for the sample allele frequency posterior distributions, as shown in Fumagalli et al. Genetics (2013). Therefore we first estimate a joint-SFS for our pair of populations:

    $NGSTOOLS/bin/ngs2dSFS -postfiles testB1.rf.sfs.norm testB2.rf.sfs.norm -outfile testB.joint.spec -relative 1 -nind 10 8 -nsites 10000 -maxlike 0    

Alternatively, as proposed in the original formulation of the method, testB1.rf.sfs and testB2.rf.sfs can be used with `-islog 1`. Another solution, for medium to high coverage data and in case of a large number of sites, would be to set `-maxlike 1`. 
You can plot the joint SFS with this simple script:

    Rscript --vanilla --slave $NGSTOOLS/plot2dSFS.R testB.joint.spec testB.joint.spec.eps pop1 pop2    

Then we calculate method-of-moments estimator of FST, at each site, with the following command:

    $NGSTOOLS/bin/ngsFST -postfiles testB1.rf.sfs testB2.rf.sfs -priorfile testB.joint.spec -nind 10 8 -nsites 10000 -outfile testB.fst -islog 1    

Note that here we use .rf.sfs instead of .rf.sfs.norm since we are using .joint.spec as a prior. In case the latter is not available (e.g. folded data or inbreeding, see further below for these examples), we can run these alternative commands (which give the same result):

    $NGSTOOLS/bin/ngsFST -postfiles testB1.rf.sfs testB2.rf.sfs -priorfiles testB1.rf.sfs.ml testB2.rf.sfs.ml -nind 10 8 -nsites 10000 -outfile testB.fst2 -islog 1    
    $NGSTOOLS/bin/ngsFST -postfiles testB1.rf.sfs.norm testB2.rf.sfs.norm -nind 10 8 -nsites 10000 -outfile testB.fst3 -islog 0    

We can calculate and print the overall FST, as well as plot FST in sliding windows, using this simple R script:

    Rscript --vanilla --slave $NGSTOOLS/plotFST.R testB.fst testB.fst.eps 100 50


## FOLDED DATA

In many cases, ancestral allelic status is unknown and analyses are carried out using folded allele frequencies. As an illustration, we show how we can compute summary statistics even in case of folded data. We use data from pop1 from previous analysis.
As usual, we use ANGSD to compute sample allele frequency posterior probabilities:

    $ANGSD/angsd -sim1 testB1.glf.gz -nInd 10 -realSFS 1 -out testC.rf -fold 1
    $ANGSD/misc/optimSFS -binput testC.rf.sfs -nChr 10 -nThreads 10
    $ANGSD/misc/sfstools -sfsFile testC.rf.sfs -nChr 10 -priorFile testC.rf.sfs.ml -dumpBinary 1 > testC.rf.sfs.norm

We now calculate summary statistics:

    $NGSTOOLS/bin/ngsStat -npop 1 -postfiles testC.rf.sfs.norm -nsites 10000 -iswin 1 -nind 10 -outfile testC.stat -isfold 1 -islog 0 -block_size 100    

Note that these values are very similar, as expected, to the one retrieved using the unfolded data (testB.stat). Differences derive from a different prior used to compute posterior probabilities.

Currently, it is not possible to estimate a joint-SFS from folded data. To compute FST from folded data, please use one of the two alternatives provided in the example above. PCA can be performed with folded data by adding `-isfold 1` in ngsCovar if we choose to weight each site by its probability of being variable.

## INBREEDING

NGSF=/home/mfumagalli/Documents/Software/ngsF

In case of data with inbreeding, almost all analyses can be carried out in the same fashion. Main difference is in how we use ANGSD and ngsF to estimate inbreeding coefficients and incorporate them into the analyses.

Simulation settings:
* 1 population, with mean inbreeding coefficient per individual of 0.3;
* 20 individuals;
* 1,000 independent sites, all polymorphic in the population (for ease of computational time);
* 4X mean sequencing depth;
* 0.01 mean sequencing error rate, variable across sites;
* minor options: minimum frequency in the population of 0.005, output files with prefix testD.

    $NGSTOOLS/bin/ngsSim -npop 1 -nind 20 -nsites 1000 -errate 0.01 -depth 4 -pvar 1 -mfreq 0.005 -F 0.3 -model 1 -outfiles testD

We estimate inbreeding coefficients and incorporate them into the calculation of posterior probabilities. Please note that simulated files need to be converted and SNP calling is required.

    $ANGSD/angsd -sim1 testD.glf.gz -nInd 20 -doGlf 3 -doSNP 1 -doMaf 2 -minLRT 15 -out testD.geno -doMajorMinor 1
    NS=`wc -l testD.geno.mafs | cut --field=1 --delimiter=' '`
    OFF=1
    NS=$((NS-OFF))
    $NGSF/ngsF -n_ind 20 -n_sites $NS -glf testD.geno.glf -verbose 0 -out testD.indF

    $ANGSD/angsd_inbreed -sim1 testD.glf.gz -nInd 20 -doGeno 32 -doPost 1 -doMaf 2 -out testD.geno -doMajorMinor 1 -indF testD.indF    
    $ANGSD/angsd_inbreed -sim1 testD.glf.gz -nInd 20 -realSFS 2 -out testD.rf -indF testD.indF -doMaf 2 -doMajorMinor 1

We can use these files to perform analyses described above. We can calculate the covariance matrix with:

    $NGSTOOLS/bin/ngsCovar -probfile testD.geno.geno -outfile testD.covar -nind 20 -nsites $NS -call 0 -sfsfile testD.rf.sfs -norm 0

and summary statistics, for instance for the whole region, with:

     $NGSTOOLS/bin/ngsStat -npop 1 -postfiles testD.rf.sfs -nsites $NS -iswin 1 -nind 20 -outfile testD.stat -isfold 0 -islog 0 -block_size $NS


An estimate of the Site Frequency Spectrum can be retrieved with these commands, if folded:

    N_IND=20
    cat testD.rf.sfs | hexdump -v -e "$((N_IND+1))/8 \"%.10g\t\"\"\n\"" | perl -na -e '$sum=0; $sum+=exp($_) for @F; for $i (0..$#F){$frq[$i]+=exp($F[$i])/$sum}; END{$tsum+=$_ for @frq; $_/=$tsum for @frq; print join("\t",@frq)."\n"}' > testD.sfs_sum

if unfolded:

    N_IND=40
    cat testD.rf.sfs | hexdump -v -e "$((N_IND+1))/8 \"%.10g\t\"\"\n\"" | perl -na -e '$sum=0; $sum+=exp($_) for @F; for $i (0..$#F){$frq[$i]+=exp($F[$i])/$sum}; END{$tsum+=$_ for @frq; $_/=$tsum for @frq; print join("\t",@frq)."\n"}' > testD.sfs_sum

Please note that FST estimation with inbreeding cannot use a joint-SFS as a prior, and therefore alternative methods proposed, like in the case of folded data, should be used.









