
// this is the final program that accepts both single priors and 2D-SFS and does correction or not

// it computes expected FST from a and b and correct for the ratio of a and b; it works in blocks of sites so it is memory efficient, it may receive in input a file of previous iteration of fst and compute the exponential weithing function and correct the post probs

/// NOTE

// if you give as input a prior2Dfile then postfile are from -realSFS 1
// if you give no priorfile I assume you have run sfstools and you do not want to use weighting correcting function
// if you give 2 priorfiles I assume you want to use weighting correcting function and you have run sfstools

// in ouput a text file, each row is a site and columsn are tab separated: a, a+b, FACT, theta, pvar

#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cstring>
#include <vector> 
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include "ngsFST.hpp"
 
// to compile: g++ ngsFST.cpp -Wall -o bin/ngsFST -lgsl -lgslcblas -lm -lz -O0 -g

int main (int argc, char *argv[]) {
  
   if (argc==1) {
    info();
   return 0;   
  }

  /// DECLARE AND INITIALIZE VARIABLES
  
  char *sfsfile1=NULL; // posterior probabilities
  char *sfsfile2=NULL;
  char *fstfile=NULL; // first guess of fst
  char *priorfile1=NULL; // priors (needed for weighting function only)
  char *priorfile2=NULL;
  char *priorfile12=NULL; // joint prior, it is 2D-SFS

  FILE *outpost;
  char *outfile=NULL;
  char *foutpost=NULL;
  
  int argPos = 1, increment = 0, nind = 0, nind1 = 0, nind2 = 0, nsites = 0, verbose = 0, nsums = 1, block_size = 0, K=0, isfold=0, firstbase=0;

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
    else if(strcmp(argv[argPos],"-priorfile")==0) {
      priorfile12 = argv[argPos+1];
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
    fprintf(stderr,"\nIf you want to correct the marginal product then perhaps you forgot to supply -priofiles and -fstfile when using an automatic setting of lambda? Ignore this message if this is the first iteration assuming independence.\n");
  }
  if((priorfile1 != NULL) & (priorfile12!=NULL)) {
    fprintf(stderr,"\nYou should give either -priorfiles or -priorfile, otherwise I don't know if you want to use a 2D-SFS or the corrected product of marginal spectra as prior.\n");
    info();
    return 0;
  }
  if((fstfile != NULL) & (priorfile12!=NULL)) {
    fprintf(stderr,"\nIf you give -fstfile I assume you want to use the correction for marginal spectra. So why are you giving -priorfile too? You should only give -priorfiles eventually if K=0.\n");
    info();
    return 0;
  }
  if((priorfile12!=NULL) & (isfold)) {
    fprintf(stderr,"\nRemember that the 2D-SFS given as a prior should be n1*2;n2*2 dimensions even if it is folded (non possible configurations should be set to 0");
   }

  /// OUTPUT
  foutpost = append(outfile, "");
  fprintf(stderr,"\t->Dumping file: %s\n", foutpost);
  outpost = getFILE(foutpost, "w");

  // print input arguments // UPDATE THIS !!!
  fprintf(stderr,"\t->Using some of these args: -nind %d -nind1 %d -nind2 %d -nsites %d -postfiles %s %s -priorfiles %s %s -priorfile %s -fstfile %s -outfile %s -verbose %d -nsums %d -offset %d -K %d\n", nind, nind1, nind2, nsites, sfsfile1, sfsfile2, priorfile1, priorfile2, priorfile12, fstfile, foutpost, verbose, nsums, firstbase, K);

  // READ PRIORS (if provided)
  // marginal spectra
  array<double> prior1;
  array<double> prior2;
  if (priorfile1 != NULL) {
    if (verbose==1) fprintf(stderr, "\nAdding priors...");
    prior1 = readArray(priorfile1, nind1, isfold);
    prior2 = readArray(priorfile2, nind2, isfold);
  }
  // 2D-SFS
  matrix<double> prior12;
  if ((priorfile12==NULL)==0) {
    if (verbose==1) fprintf(stderr, "\nAdding 2D prior...");
    prior12 = readPrior12(priorfile12, nind1*2+1, nind2*2+1); // same dimensions even if it is folded
    if (verbose==2) {
      fprintf(stderr, "\nPrior 2d:\n");
      writematrix(prior12, stderr);
    }
    //// the difference with this prior is that I don't add the prior, but I add the prior directly at computeFST step
  }


  /// GET POSITIONS OF BLOCKS
  if (block_size>(nsites-firstbase+1)) block_size=(nsites-firstbase+1);
  if (block_size==0) block_size=nsites-firstbase+1;
  array<int> start; array<int> end;
  start=getStart(nsites, firstbase, block_size);
  end=getEnd(nsites, firstbase, block_size);
  int nwin = start.x;

  if (verbose==1) fprintf(stderr, "\n num win %d win0 is %d %d", nwin, start.data[0], end.data[0]);

  /// ITERATE OVER EACH BLOCK
  for (int n=0; n<nwin; n++) {

    fprintf(stderr, "Block %d out of %d from %d to %d\n", n, (nwin-1), start.data[n], end.data[n]);
  
    // READ POSTERIOR PROBABILITIES FILES
    matrix<double> post1;
    matrix<double> post2;
    post1 = readFileSub(sfsfile1, nind1, start.data[n], end.data[n], isfold);
    post2 = readFileSub(sfsfile2, nind2, start.data[n], end.data[n], isfold);

    if (priorfile12!=NULL) {
      // if from -realSFS 1 they are -log
      normSFS(post1, 1); // 2nd argument is islog
      normSFS(post2, 1);
    }

    // IF NOT FST FILE PROVIDED
    if ((fstfile == NULL)) {
       if (verbose==1) fprintf(stderr,"Computing FST and no first guess provided.\n");
       if (priorfile12==NULL) {
         if (isfold) {
           computeVarReyFold(post1, post2, verbose, outpost, nsums);
         } else {
           computeVarRey(post1, post2, verbose, outpost, nsums);
         }
       } else {
         if (verbose==1) fprintf(stderr,"Using 2D-SFS as prior. You didn't run sfstools, right??? Use only -realSFS 1.\n");
         computeVarRey12New(post1, post2, verbose, outpost, nsums, prior12);
       }
     } else {
     // IF FST FILE IS INDEED PROVIDED
     if (verbose==1) fprintf(stderr,"Computing FST and first guess provided.\n");
       array<double> firstfst;
       firstfst=readFSTsub(fstfile, nsites, start.data[n], end.data[n]);
       array <double> sublam;
       sublam = getLambdas(firstfst, prior1, prior2, K, verbose, isfold);
       if (verbose==1) fprintf(stderr,"Computed lambdas.\n");
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
 
  delete [] start.data;
  delete [] end.data;

  if ((priorfile12==NULL)==0) cleanup(prior12);

  if ((priorfile1==NULL)==0) delete [] prior1.data;
  if ((priorfile2==NULL)==0) delete [] prior2.data;

  fclose(outpost);
  free(foutpost);
  
  return 0;

} // end main





