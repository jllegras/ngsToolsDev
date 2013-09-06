
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
    tmp[i] = (atoi(strtok(NULL,"\t \n")));
  fclose(fp);
  ret.x = len;
  ret.data = tmp;
  return ret;
}

// read a file into a matrix but only for a specific subsets of positions (0-based notation)
matrix<double> readFileSub(char *fname, int nInd, int start, int end, int isfold) {
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if (isfold==0) {
    if((filesize %(sizeof(double)*(2*nInd+1)) )) {
      fprintf(stderr,"\n\t-> Possible error,binaryfiles might be broken\n");
      exit(-1);
    }
  } else {
    if((filesize %(sizeof(double)*(nInd+1)) )) {
      fprintf(stderr,"\n\t-> Possible error,binaryfiles might be broken\n");
      exit(-1);
    }
  }
  int nsites = end-start;
  double **data = new double*[nsites];
  if (isfold) {
          fseek(fp, sizeof(double)*(nInd+1)*start, SEEK_SET);
  } else {
          fseek(fp, sizeof(double)*(2*nInd+1)*start, SEEK_SET);
  }
  if (isfold) {
    for(int i=0; i<nsites; i++) {
      double *tmp = new double[nInd+1];
      fread(tmp,sizeof(double),nInd+1,fp);
      data[i]= tmp;
    }
  } else {
    for(int i=0; i<nsites; i++) {
      double *tmp = new double[2*nInd+1];
      fread(tmp,sizeof(double),2*nInd+1,fp);
      data[i]= tmp;
    }
  }
  fclose(fp);
  matrix<double> ret;
  ret.x = nsites;
  if (isfold) {
    ret.y = nInd+1;
  } else {
    ret.y = 2*nInd+1;
  }
  ret.data = data;
  return ret;
}

int main (int argc, char *argv[]) {

  char *infile=NULL;
  char *posfile=NULL; 

  char *outfile;
  char *foufile=NULL;

  int argPos = 1, nind = 0, nsites = 0, isfold=0, increment=0, verbose=0, len=0;

  // CHECK HERE FOR USAGE

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
    else if(strcmp(argv[argPos],"-isfold")==0) // whether the .sfs. is folded or not 
      isfold = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-verbose")==0)
      verbose = atoi(argv[argPos+1]);
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      return 0;
    }
    argPos = argPos + 2 + increment;
  } 

  //outfile = getFILE(foutest, "wb");

  matrix<double> sfs;
  sfs = readFileSub(infile, nind, 0, nsites, isfold);
  if (verbose) fprintf(stderr, "Dim input %d , %d; example %f %f\n", sfs.x, sfs.y, sfs.data[0][0], sfs.data[1][1]);

  array<int> pos;
  pos = readArray(posfile, len);
  if (verbose) fprintf(stderr, "Dim pos %d; example %d %d \n", pos.x, pos.data[0], pos.data[pos.x-1]);

  matrix<double> new_sfs;
  double **cdata = new double*[pos.x];
  for(int i=0;i<pos.x;i++){
    double *ctmp = new double[sfs.y];
    cdata[i]= ctmp;
  }
  new_sfs.x=pos.x;
  new_sfs.y=sfs.y;
  new_sfs.data = cdata;
  fprintf(stderr, "Dim output %d , %d", new_sfs.x, new_sfs.y);
  for (int i=0; i<new_sfs.x; i++) {
    for (int j=0; j<new_sfs.y; j++) {
      new_sfs.data[i][j]=0.0;
    }
  }
  if(verbose) fprintf(stderr, "\nDim output %d , %d; example %f %f;", new_sfs.x, new_sfs.y, new_sfs.data[0][0], new_sfs.data[1][1]);

  for (int i=0; i<new_sfs.x; i++) {
      for (int j=0; j<new_sfs.y; j++) {
        //fprintf(stderr, "%d %f\t", j, sfs.data[pos.data[i]][j]);
        new_sfs.data[i][j]=sfs.data[(pos.data[i]-1)][j]; // -1 because it is 1-based to 0-based
      }
  }

  if (verbose) fprintf(stderr, "\nDim output %d , %d; example %f %f;", new_sfs.x, new_sfs.y, new_sfs.data[0][0], new_sfs.data[1][1]);

  FILE *fp = fopen(outfile,"wb");

  for (int i=0; i<new_sfs.x; i++)
    fwrite(new_sfs.data[i], sizeof(double), new_sfs.y, fp);

  fclose(fp);
  //fclose(infile);
  //fclose(posfile);

} // main


