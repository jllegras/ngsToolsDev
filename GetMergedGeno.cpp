
// in input it receives 2 .geno files
// in output it generates a new merged .geno file (binary)

#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cstring>
#include <vector>
#include <math.h>

#include <string> //for str operations
using namespace std;

// a general matrix style structure
template <typename T>
struct matrix{
  int x;
  int y;
  T** data;
};

template <typename T>
struct array{
  int x;
  T* data;
};

template <typename T>
T *collapse(std::vector<T> &v){
  T *tmp = new T[v.size()];
  for(int i=0;i<v.size();i++)
    tmp[i] = v[i];
  return tmp;
}

//function to cleanup our generic matrix structure
template <typename T>
void cleanup(matrix<T> &m){//using a reference to avoid copying the data
  for(int i=0;i<m.x;i++)
    delete [] m.data[i];
  delete [] m.data;
}

// get the filesize of a file
size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}

// find out if a file exists
int fexists(const char* str) {
  struct stat buffer ;
  return (stat(str, &buffer )==0 );
}

// a nice wrapper for getting files
FILE *getFILE(const char*fname,const char* mode) {
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;
  if(writeFile&&fexists(fname)){//DRAGON
    fprintf(stderr,"\t-> File exists: %s exiting...\n",fname);
    exit(0);
  }
  FILE *fp;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t->Error opening FILE handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

// read a file into a matrix
matrix<double> readFile(char *fname, int nInd, int nSites) {

  FILE *fp = getFILE(fname,"rb");
  size_t filesize =fsize(fname);

  if((filesize %(sizeof(double)*nInd*3*nSites)) ) {
      fprintf(stderr,"\n\t-> Possible error,binaryfiles might be broken\n");
      exit(-1);
  }

  double **data = new double*[nSites*nInd];
  fseek(fp, 0, SEEK_SET);
  
  for(int i=0; i<(nSites*nInd); i++) {
    double *tmp = new double[3];
    fread(tmp,sizeof(double),3,fp);
    data[i]= tmp;
  }

  fclose(fp);
  matrix<double> ret;
  ret.x = nSites*nInd;
  ret.y = 3;
  ret.data = data;
  return ret;

}

int main (int argc, char *argv[]) {

  char *infile1=NULL;
  char *infile2=NULL;

  char *outfile;
  char *foufile=NULL;

  int argPos = 1, nind = 0, nind1=0, nind2=0, nsites = 0, increment=0, verbose=0, len=0;

  // CHECK HERE FOR USAGE:

  while (argPos<argc) {
    increment = 0;
    if(strcmp(argv[argPos],"-infiles")==0) {
      infile1 = argv[argPos+1];
      infile2 = argv[argPos+2];
      increment = 1;
    } else if(strcmp(argv[argPos],"-nind")==0) {
      nind1 = atoi(argv[argPos+1]);
      nind2 = atoi(argv[argPos+2]);
      increment = 1;
    } else if(strcmp(argv[argPos],"-nsites")==0)
      nsites = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-outfile")==0)
      outfile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-verbose")==0)
      verbose = atoi(argv[argPos+1]);
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      return 0;
    }
    argPos = argPos + 2 + increment;
  } 

  nind=nind1+nind2;

  matrix<double> geno1;
  geno1 = readFile(infile1, nind1, nsites);
  if (verbose) fprintf(stderr, "Dim input %d , %d; example %f %f %f\n", geno1.x, geno1.y, geno1.data[0][0], geno1.data[0][1], geno1.data[0][2]);

  matrix<double> geno2;
  geno2 = readFile(infile2, nind2, nsites);
  if (verbose) fprintf(stderr, "Dim input %d , %d; example %f %f %f\n", geno2.x, geno2.y, geno2.data[0][0], geno2.data[0][1], geno2.data[0][2]);

  // initialize
  int new_nrow=nsites*nind; 
  matrix<double> new_geno;
  double **cdata = new double*[new_nrow];
  for(int i=0;i<new_nrow;i++){
    double *ctmp = new double[geno1.y];
    cdata[i]= ctmp;
  }
  new_geno.x=new_nrow;
  new_geno.y=geno1.y;
  new_geno.data = cdata;
  for (int i=0; i<new_geno.x; i++) {
    for (int j=0; j<new_geno.y; j++) {
      new_geno.data[i][j]=0.0;
    }
  }
  if (verbose) fprintf(stderr, "Dim output %d , %d; example %f %f %f\n", new_geno.x, new_geno.y, new_geno.data[0][0], new_geno.data[0][1], new_geno.data[0][2]);

  // merge
  int conta=0;

  for (int s=0; s<nsites; s++) {

    for (int i=(s*nind1); i<((s*nind1)+(nind1)); i++) {
      for (int j=0; j<geno1.y; j++) {
        new_geno.data[conta][j]=geno1.data[i][j];
      }
      conta++;
    }

    for (int i=(s*nind2); i<((s*nind2)+(nind2)); i++) {
      for (int j=0; j<geno2.y; j++) {
        new_geno.data[conta][j]=geno2.data[i][j];
      }
      conta++;
    }

  }

  if (verbose) fprintf(stderr, "Dim output %d , %d; example %f %f %f\n", new_geno.x, new_geno.y, new_geno.data[0][0], new_geno.data[0][1], new_geno.data[0][2]);

  if ( conta != (new_geno.x)) {
    fprintf(stderr, "Possible error, the counts and the new geno are not the same %d %d\n", conta, new_geno.x);
    exit(-1);
  }

  FILE *fp = fopen(outfile,"wb");

  for (int i=0; i<new_geno.x; i++)
    fwrite(new_geno.data[i], sizeof(double), new_geno.y, fp);

  fclose(fp);
  //fclose(infile);
  //fclose(posfile);

} // main


