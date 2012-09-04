
// this is the final program that accepts both single priors and 2D-SFS and does correction or not

// mavaf : computes expected FST from a and b and correct for the ratio of a and b; it works in blocks of sites so it is memory efficient, recieves in input a file of previous iteration of fst and compute the exponential weithing function and correct the post probs

// ultimate version: receives in input sfstools (already normalized)

#include <cstdio> //for stderr,stdout
#include <cstdlib> //for atoi
#include <sys/stat.h> //for getting file attributes
#include <cstring> //for str operations
#include <vector> 
#include <math.h> // for exponential and log
// include templates, functions
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include "mavaf_ultimate_fast.hpp"

// fast version is the one optimizing only one time for equal values of firts guess of fst

// it receives in input only post probs files
    
  // input:
  //	files with posterior probabilities for each site for each population
  //	name of output file
  //	nr of indiv for each pop
  //	nr of sites // how many sites you want to analyze
  //    verbose level ?
  //    number of sums of moments ? use 1
  //    is output a text or binary? // to do...
  //	block_size : how many sites per block?
  //	an optional file of first guess of fst (from mavaf itself)
  //	firstbase: offset , start reading from that base
  //	isfold: is data folded?

  // output:
  //	a, ab, correction factor, fst valus, pvar [all for each site]
  // (temporarly the output file is a text file)
 
  // then a locus-estimate could be simply Reynolds-style so: sum(a)/sum(a+b) without correction

 // to compile: g++ mavaf_ultimate_fast.cpp -Wall -o mavaf_ultimate_fast -lgsl -lgslcblas -lm -lz -O0 -g

 // usage
 // e.g ./mavaf -postfiles POPS1.sfs POPS2.sfs -nind 10 20 -nsites 20000 -outfile out -verbose 0 -nsums 1 -block_size 1000 -firstbase 0 -fstfile myfst.txt

int main (int argc, char *argv[]) {
  
   if (argc==1) {
    info();
   return 0;   
  }

  /// DECLARE AND INITIALIZE VARIABLES
  
  char *sfsfile1=NULL; // posterior probabilities
  char *sfsfile2=NULL;
  char *fstfile=NULL; // first guess of fst
  char *priorfile1=NULL; // posterior probabilities
  char *priorfile2=NULL;

  FILE *outpost;
  char *outfile=NULL;
  char *foutpost=NULL;
  
  int argPos = 1, increment = 0, nind = 0, nind1 = 0, nind2 = 0, nsites = 0, verbose = 0, nsums = 1, block_size = 10000, K=0, isfold=0, firstbase=0;

  /// READ AND ASSIGN INPUT PARAMETERS
  
   while (argPos<argc) {
    increment = 0;
    if(strcmp(argv[argPos],"-postfiles")==0) {
      sfsfile1 = argv[argPos+1];
      sfsfile2 = argv[argPos+2];
      increment = 1;
    }
    else if(strcmp(argv[argPos],"-fstfile")==0) {
      fstfile = argv[argPos+1];
    }
    else if(strcmp(argv[argPos],"-priorfiles")==0) {
      priorfile1 = argv[argPos+1];
      priorfile2 = argv[argPos+2];
      increment = 1;
    }     
    else if(strcmp(argv[argPos],"-outfile")==0) outfile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-nind")==0) {
      nind1 = atoi(argv[argPos+1]);
      nind2 = atoi(argv[argPos+2]);
      nind = nind1 + nind2;
      increment = 1;
    }
    else if(strcmp(argv[argPos],"-nsites")==0) nsites = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-K")==0) K = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-verbose")==0) verbose = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-block_size")==0) block_size = atoi(argv[argPos+1]);    
    else if(strcmp(argv[argPos],"-nsums")==0) nsums = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-firstbase")==0) firstbase = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-isfold")==0) isfold = atoi(argv[argPos+1]);
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0; // terminate
    }
    argPos = argPos + 2 + increment;
  }
  
  /// CHECK INPUT
  if((sfsfile1 == NULL) & (sfsfile2 == NULL) ) {
    fprintf(stderr,"\nMust supply -postfiles.\n");
    info();
    return 0;
  }
  if(outfile == NULL) {
    fprintf(stderr,"\nMust supply -outfile.\n");
    info();
    return 0;
  }
  if((priorfile1 == NULL) & (fstfile==NULL) & (K==0)) {
    fprintf(stderr,"\nPerhaps you forgot to supply -priofiles when using an automatic setting of lambda?\n");
    info();
    return 0;
  }

  /// OUTPUT
  // print input arguments
  fprintf(stderr,"\t->Using args: -nind %d -nind1 %d -nind2 %d -nsites %d -postfiles %s %s -outfile %s -verbose %d -nsums %d \n", nind, nind1, nind2, nsites, sfsfile1, sfsfile2, foutpost, verbose, nsums);
  // prepare output file (temporarly is a text file)
  foutpost = append(outfile, ".txt");
  fprintf(stderr,"\t->Dumping file: %s\n", foutpost);
  outpost = getFILE(foutpost, "w");

  // READ PRIORS (if provided)
  array<double> prior1;
  array<double> prior2;
  if (priorfile1 != NULL) {
    if (verbose==1) fprintf(stderr, "\nAdding priors...");
    prior1 = readArray(priorfile1, nind1, isfold);
    prior2 = readArray(priorfile2, nind2, isfold);
    for (int i=0; i<prior1.x; i++) { 
      if(prior1.data[i]<0.000001) prior1.data[i]=0.000001;
    }
    for (int i=0; i<prior2.x; i++) { 
      if(prior2.data[i]<0.000001) prior2.data[i]=0.000001;
    }
  }

  /// GET POSITIONS OF BLOCKS
  array<int> start; array<int> end; 
  start=getStart(nsites, firstbase, block_size);
  end=getEnd(nsites, firstbase, block_size);
 
  /// ITERATE OVER EACH BLOCK
  int nwin = start.x;
  for (int n=0; n<nwin; n++) {

    fprintf(stderr, "Block %d out of %d from %d to %d\n", n, (nwin-1), start.data[n], end.data[n]);
  
    // READ POSTERIOR PROBABILITIES FILES
    matrix<double> post1;
    matrix<double> post2;
    post1 = readFileSub(sfsfile1, nind1, start.data[n], end.data[n], isfold);
    post2 = readFileSub(sfsfile2, nind2, start.data[n], end.data[n], isfold);

    // IF NOT FST FILE PROVIDED
    if ((fstfile == NULL)) {

      // COMPUTE FST
      if (verbose==1) fprintf(stderr,"Computing FST and no first guess provided.\n");

       if (isfold) {
         computeVarReyFold(post1, post2, verbose, outpost, nsums);
       } else {
         computeVarRey(post1, post2, verbose, outpost, nsums);
       }

     } else {
     // IF FST FILE IS INDEED PROVIDED
     if (verbose==1) fprintf(stderr,"Computing FST and first guess provided.\n");

       // read file
       array<double> firstfst;
//       fprintf(stderr, "\n read fst file...");
       firstfst=readFSTsub(fstfile, nsites, start.data[n], end.data[n]);
 //      fprintf(stderr, "...done!");

       // COMPUTE LAMBDAS for only this subset of sites (and norm for post probs! not prior, a change here).

//       fprintf(stderr, "\n first fst: %d %f %f\n", firstfst.x, firstfst.data[0], firstfst.data[99]);
 //      writearray(firstfst, stderr);
 //      printf("\n");  

       array <double> sublam;
       sublam = getLambdas(firstfst, prior1, prior2, K, verbose, isfold);

       //if (verbose==1) fprintf(stderr, "\n end get lambdas %d: %f %f %f %f\n", sublam.x, sublam.data[0], sublam.data[1], sublam.data[2], sublam.data[sublam.x-1]);
       //writearray(sublam, stderr);

       if (isfold) {
         computeVarRey2Fold(post1, post2, verbose, outpost, nsums, sublam);
       } else {
         computeVarRey2(post1, post2, verbose, outpost, nsums, sublam);
       }
 
       delete [] firstfst.data;
       delete [] sublam.data;

    }    

    cleanup(post1);
    cleanup(post2);
    
  } // end blocks iterations

  delete [] prior1.data;
  delete [] prior2.data;

  delete [] start.data;
  delete [] end.data;

  free(foutpost);
   
  return 0;

} // end main





