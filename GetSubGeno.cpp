
// in input it receives a file .txt lines / sites to be kept (1 based)
// in output it generates a new .geno file (binary) with only the sites retained

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

array<int> readArray(const char *fname, int len) {
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if(filesize==0){
    fprintf(stderr,"file:%s looks empty\n",fname);
    exit(0);
  }
  array<int> ret;
  int *tmp = new int[len];
  char *buf = new char[filesize];
  fread(buf,sizeof(char),filesize,fp);
  tmp[0] = atoi(strtok(buf,"\t \n"));
  for(int i=1;i<(len);i++)
    tmp[i] = atoi(strtok(NULL,"\t \n"));
  fclose(fp);
  ret.x = len;
  ret.data = tmp;
  return ret;
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

  char *infile=NULL;
  char *posfile=NULL; 

  char *outfile;
  char *foufile=NULL;

  char *genoquality=NULL;

  int argPos = 1, nind = 0, nsites = 0, increment=0, verbose=0, len=0;

  while (argPos<argc) {
    increment = 0;
    if(strcmp(argv[argPos],"-infile")==0)
      infile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-posfile")==0)
      posfile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-nind")==0)
      nind = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-nsites")==0)
      nsites = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-len")==0)
      len = atoi(argv[argPos+1]);
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

  matrix<double> geno;
  geno = readFile(infile, nind, nsites);
  if (verbose) fprintf(stderr, "Dim input %d , %d; example %f %f\n", geno.x, geno.y, geno.data[0][0], geno.data[1][1]);

  // posfile is a file newline-separated with numbers of lines/sites to keep
  array<int> pos;
  pos = readArray(posfile, len);
  if (verbose) fprintf(stderr, "Dim pos %d; example %d %d \n", pos.x, pos.data[0], pos.data[1]);

  // initialize
  int new_nrow=len*nind; 
  matrix<double> new_geno;
  double **cdata = new double*[new_nrow];
  for(int i=0;i<new_nrow;i++){
    double *ctmp = new double[geno.y];
    cdata[i]= ctmp;
  }
  new_geno.x=new_nrow;
  new_geno.y=geno.y;
  new_geno.data = cdata;
  if (verbose) fprintf(stderr, "Dim output %d , %d", new_geno.x, new_geno.y);
  for (int i=0; i<new_geno.x; i++) {
    for (int j=0; j<new_geno.y; j++) {
      new_geno.data[i][j]=0.0;
    }
  }
  if (verbose) fprintf(stderr, "\nDim output %d , %d; example %f %f;", new_geno.x, new_geno.y, new_geno.data[0][0], new_geno.data[1][1]);

  // fill
  // pos is 1 based
  int row=0;
  for (int s=0; s<pos.x; s++) {
    row=pos.data[s];
    for (int i=((row-1)*nind); i<(((row-1)*nind)+(nind-1)); i++) {
      for (int j=0; j<geno.y; j++) {
        new_geno.data[s][j]=geno.data[pos.data[i]][j];
      }
    }
  }
  
  if (verbose) fprintf(stderr, "\nDim output %d , %d; example %f %f;", new_geno.x, new_geno.y, new_geno.data[0][0], new_geno.data[1][1]);

  FILE *fp = fopen(outfile,"wb");

  for (int i=0; i<new_geno.x; i++)
    fwrite(new_geno.data[i], sizeof(double), new_geno.y, fp);

  fclose(fp);
  //fclose(infile);
  //fclose(posfile);

} // main


