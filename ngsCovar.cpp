#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cstring>
#include <vector>
#include <math.h>
#include "ngsCovar.hpp"

// to compile: g++ -Wall -O0 -g ngsCovar.cpp -o ngsCovar

// implement offset

int main (int argc, char *argv[]) {
  
  // CALL HELP FUNCTION
  if (argc==1) {
    info();
    return 0; // terminate
  }
 
  // DECLARE
 
  // possible inputs
  char *estfile=NULL; // estimated genotypes with probs (from angsd -doGeno 64)
  char *sfsfile=NULL; // sfs probs (-realSFS 1 + optimSFS = sfstools [already normalized and not log]); this is a major change 
  
  FILE *outest;
  char *outfiles=NULL;
  char *foutest=NULL;
  
  int argPos = 1, increment = 0, nind = 0, nsites = 0, debug = 0, norm = 0, block_size = 10000, call=0, offset=0, maxgeno=0;
  double esites = 0.0, minmaf = 0.0;
 
  /// READ AND ASSIGN INPUT PARAMETERS
  
  while (argPos<argc) {
    increment = 0; // leave this in case of future implementation of more parameters per argument
    if(strcmp(argv[argPos],"-probfile")==0) 
      estfile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-sfsfile")==0) 
      sfsfile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-nind")==0) 
      nind = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-nsites")==0) 
      nsites = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-outfile")==0) 
      outfiles = argv[argPos+1];
    else if(strcmp(argv[argPos],"-norm")==0) 
      norm = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-minmaf")==0) 
      minmaf = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-block_size")==0)
      block_size = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-offset")==0)
      offset = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-call")==0)
      call = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-verbose")==0) 
      debug = atoi(argv[argPos+1]);
    else { // input is not a valid one 
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0; // terminate
    }
    argPos = argPos + 2 + increment;
  } // end while all inputs
  
  // check if there is the input file
  if ((estfile == NULL)) {
    fprintf(stderr,"\nMust supply -probfile.\n");
    info();
    return 0;
  }

  // check if there is the sfs file
  if ((sfsfile != NULL) & (minmaf>0)) {
    fprintf(stderr,"\n-sfsfile and -minmaf are in conflict! If use use -sfsfile to weight each file there is no need (theoretically) to filter out sites with low maf. In this case -minmaf won't be applied and the program terminates here.\n");
    info();
    return 0;
  }
  
  // check if there is the output file
  if(outfiles == NULL) {
    fprintf(stderr,"\nMust supply -outfile.\n");
    info();
    return 0;
  }

  // if block_size longer than nsites
  if (block_size>nsites) block_size=nsites;

  // prepare output files
  foutest = append(outfiles, "");
  
  // print input arguments
  fprintf(stderr,"\t->Using args: -nind %d -nsites %d -probfile %s -sfsfile %s -outfile %s -verbose %d -norm %d -minmaf %f -block_size %d -call %d -offset %d\n", nind, nsites, estfile, sfsfile, foutest, debug, norm, minmaf, block_size, call, offset); 
 
  // BLOCKS
  int nwin = (nsites/block_size);
  /// GET POSITIONS OF BLOCKS
  array<int> start; array<int> end;
  start=getStart(nsites, offset, block_size);
  end=getEnd(nsites, offset, block_size);
  int maxlen=end.data[0]-start.data[0]+1; // len for each win, it will never be greater than this
 
  // prepare out
  fprintf(stderr,"\t->Dumping file: %s\n", foutest);
  outest = getFILE(foutest, "w");

  // initialize covariance matrix  
  matrix<double> covar;
  double **cdata = new double*[nind];
  for(int i=0;i<nind;i++){
    double *ctmp = new double[nind];
    cdata[i]= ctmp;
  }
  covar.x = nind;
  covar.y = nind;
  covar.data = cdata;
  for (int i=0;i<nind;i++) {
    for (int j=0;j<nind;j++) {
      covar.data[i][j]=0.0;
    }
  }

  // initialize pvar
  array<double> pvar;
  double *temp2 = new double [maxlen]; // init to the size of the largest window
  for (int i=0; i<maxlen; i++)
   temp2[i]=0.0;
  pvar.x=maxlen;
  pvar.data=temp2;

  // ITERATING for each block
  for (int n=0; n<nwin; n++) {

    matrix<double> sfs;
    matrix<double> esti;
    array<double> pp;

    fprintf(stderr, "Block %d out of %d from %d to %d\n", n, (nwin-1), start.data[n], end.data[n]);

    // compute esti
    if (debug==1) fprintf(stderr, "\nGetting esti...");
    esti = readEstiSub(estfile, nind, start.data[n], end.data[n]);
    if (debug==2) writematrix(esti, stderr);
    if (debug==1) fprintf(stderr, ": %d %d , %f %f", esti.x, esti.y, esti.data[0][0], esti.data[1][1]);  

    // IF CALL GENOTYPES, set max prob to 1
    if (call) {
      if (debug==1) fprintf(stderr, "\nCalling genotypes...");
      for (int i=0; i<esti.x; i++) {
        maxgeno=maxposarr(esti, i);
        for (int j=0; j<3; j++) esti.data[i][j]=0.0;
        esti.data[i][maxgeno]=1.0;
      }
    }

    /// COMPUTE ALLELE FREQUENCIES
    if (debug==1) fprintf(stderr, "\nGetting mu...");
    pp = getAlleFreq(esti);
    if (debug==3) writearray(pp, stderr);
    // this is NOT strictly the alle freq, is mu in Patterson notation, thus it's the average number of derived alleles per individual; the allele freq p is simply mu/2
    if (debug==2) {
      for (int i=0; i<pp.x; i++)
       if (pp.data[i]<0) fprintf(stderr, "\n %d %f", i, pp.data[i]);
    }
    if (debug==1) fprintf(stderr, ": %d, %f %f ", pp.x, pp.data[0], pp.data[4]);
 
    /// COMPUTE COVARIANCE, here simply update
    if (debug==1) fprintf(stderr, "\nUpdating covar...");
    if (sfsfile!=NULL) {
      // read sfs
      fprintf(stderr, "...weighting...");
      sfs = readFileSub(sfsfile, nind, start.data[n], end.data[n], 0);
      if (debug==1) fprintf(stderr, "\nGot  sfs: %d %d, e.g. %f %f", sfs.x, sfs.y, sfs.data[0][0], sfs.data[1][1]);
      if (debug==1) fprintf(stderr, "\nGetting pvar...");
      getPvar(sfs, pvar);
      if (debug==1) fprintf(stderr, ": %d , %f %f", pvar.x, pvar.data[0], pvar.data[10]);
      cleanup(sfs);
      if (debug==1) fprintf(stderr, "\nUpdating covar...");
      calcCovarUpProb(esti, pp, norm, covar, pvar);
      if (debug==1) fprintf(stderr, "\nUpdating esite...");
      // divide by expected nr of segr sites -1 in case of pvar
      for (int i=0;i<pvar.x;i++) esites = esites + pvar.data[i];
      if (debug==1) fprintf(stderr, ": %f", esites);
    } else {
      if (debug==1) fprintf(stderr, "\n I am using this minmaf %f ", minmaf);
      fprintf(stderr, "...no weighting..."); 
      double tmp_eff_nsites=calcCovarUp(esti, pp, norm, covar, minmaf);
      esites=esites+tmp_eff_nsites;
    }
       
    cleanup(esti);

    delete [] pp.data;
  
  } // end for n in nwin block

  fprintf(stderr, "\n(Exp/eff) nr sites: %f\n", esites);

  delete [] pvar.data;

  // divide
  for (int i=0;i<nind;i++) {
    for (int j=0;j<nind;j++) {
      covar.data[i][j]=covar.data[i][j]/(esites-1);
    }
  }

  writematrix(covar, outest);  
  
  // free
  cleanup(covar);
  free(foutest);
    
  return 0;
  
} // end main

