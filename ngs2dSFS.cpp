
// a program to estimate 2D-SFS from genotype posterior probabilities

#include <cstdio> //for stderr,stdout
#include <cstdlib> //for atoi
#include <sys/stat.h> //for getting file attributes
#include <cstring> //for str operations
#include <vector> 
#include <math.h> // for exponential and log
#include "ngs2dSFS.hpp" // include templates, functions

// input is the output of sfstools

// to compile: g++ -Wall -O0 -g ngs2dSFS.cpp -o ngs2dSFS

int main (int argc, char *argv[]) {
  
  if (argc==1) {
    printf("\n// this is to estimate 2D-SFS from geno post probabilities (from sfstools) // input: -postfiles, -outfile, -relative, -nind, -nsites, -block_size, -offset\n");
    return 0;    
  }

  /// DECLARE AND INITIALIZE VARIABLES
  
  char *sfsfile1=NULL;
  char *sfsfile2=NULL;
  
  FILE *outpost;
  char *outfile=NULL;
  char *foutpost=NULL;
  
  int argPos = 1;
  int increment = 0;
  int nind1 = 0, nind2 = 0, nsites = 0;
  int block_size = 10000;
  int firstbase = 0;
  int relative = 1;
  int maxlike=1;

  // READ AND ASSIGN INPUT PARAMETERS
  
   while (argPos<argc) {
    increment = 0;
    if(strcmp(argv[argPos],"-postfiles")==0) {
      sfsfile1 = argv[argPos+1];
      sfsfile2 = argv[argPos+2];
      increment = 1;
    }
    else if(strcmp(argv[argPos],"-outfile")==0) outfile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-nind")==0) {
      nind1 = atoi(argv[argPos+1]);
      nind2 = atoi(argv[argPos+2]);
      increment = 1;
    }
    else if(strcmp(argv[argPos],"-nsites")==0) nsites = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-block_size")==0) block_size = atoi(argv[argPos+1]);    
    else if(strcmp(argv[argPos],"-offset")==0) firstbase = atoi(argv[argPos+1]); 
    else if(strcmp(argv[argPos],"-maxlike")==0) maxlike = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-relative")==0) relative = atoi(argv[argPos+1]);
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      return 0;
    }
  
    argPos = argPos + 2 + increment;
  
  } // end while all inputs

  // check if there is the input file
  if((sfsfile1 == NULL) & (sfsfile2 == NULL) ) {
    fprintf(stderr,"\nMust supply -postfiles.\n");
    return 0;
  }
  // check if there is the output file
  if(outfile == NULL) {
    fprintf(stderr,"\nMust supply -outfile.\n");
    return 0;
  }

  // prepare output file
  foutpost = append(outfile, "");
  outpost = getFILE(foutpost, "w");

  // output
  matrix<double> spec;
  spec.x=(2*nind1)+1;
  spec.y=(2*nind2)+1;

  double **data = new double*[spec.x];
  for (int i=0; i<spec.x; i++) {
    double *tmp = new double[spec.y];
    for (int j=0; j<spec.y; j++) {
      tmp[j] = 0.0;
     }
    data[i] = tmp;
  }
  spec.data = data;


  /// READ POSTERIOR PROBABILITY FILES

  // BLOCKS

  if (block_size>(nsites-firstbase)) block_size=nsites;

  int nwin = ((nsites-firstbase)/block_size);
  array<int> start, end;
    
  if ( ((nsites-firstbase) % block_size)!=0) {
    start.x=nwin+1; end.x=nwin+1;
  } else {
    start.x=nwin, end.x=nwin;
  }
  int *tStart= new int [nwin]; 
  int *tEnd= new int [nwin]; 
    
  for (int i=0; i<nwin; i++) {
    tEnd[i]=(i+1)*block_size-1;
    tStart[i]=(i)*block_size;
  }
  // if there is rest
  if ( (nsites % block_size)!=0) {
    tEnd[nwin]=nsites-1;
    tStart[nwin]=tStart[(nwin-1)]+block_size;
    nwin=nwin+1;
  }
  // if you dont start from beginning
  if (firstbase>0) {
    for (int i=0; i<nwin; i++) {
      tEnd[i]=tEnd[i]+firstbase;
      tStart[i]=tStart[i]+firstbase;
    }
  }
  start.data=tStart;
  end.data=tEnd;

  matrix<double> post1;
  matrix<double> post2; 

  // for each block
  for (int n=0; n<nwin; n++) {
      
    fprintf(stderr, "win %d out of %d from %d to %d", n, (nwin-1), start.data[n], end.data[n]);
    post1 = readFileSub(sfsfile1, nind1, start.data[n], end.data[n], 0);
    post2 = readFileSub(sfsfile2, nind2, start.data[n], end.data[n], 0);
  
    // COMPUTE SFS
    sumSpectrum(spec, post1, post2, maxlike);

    // clean
    cleanup(post1);
    cleanup(post2);

  } // end for in nwin

  delete [] start.data;
  delete [] end.data;

  // if relative divide by nsites
  if (relative) {
    for (int i=0; i<spec.x; i++) {
      for (int j=0; j<spec.y; j++) {
        spec.data[i][j]=spec.data[i][j] / (nsites-firstbase);
      }
    }
  }

  // write
  writematrix(spec, outpost);
  cleanup(spec);
  
  // close file
  free(foutpost);
   
  return 0;

} // end main
