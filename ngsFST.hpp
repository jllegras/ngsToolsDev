
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>


/// TEMPLATES

// a general matrix style structure
template <typename T>
struct matrix{
  int x;
  int y;
  T** data;
};

// a general array style structure
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

// struct for optimization (gsl_function data type)
typedef struct {
  double k;
  double z;
  double F;
  array<double> p1;
  array<double> p2;
} my_f_params;

/// FUNCTIONS

array<int> getStart(int nsites, int firstbase, int block_size) {
  int nwin = ((nsites-firstbase)/block_size);
  array<int> start;
  if ( ((nsites-firstbase) % block_size)!=0) {
    start.x=nwin+1;
  } else {
    start.x=nwin;
  }
  int *tStart= new int [nwin]; 
  for (int i=0; i<nwin; i++) {
    tStart[i]=(i)*block_size;
  }
  // if there is rest
  if ( (nsites % block_size)!=0) {
    tStart[nwin]=tStart[(nwin-1)]+block_size;
    nwin=nwin+1;
  }
  // if you dont start from beginning
  if (firstbase>0) {
    for (int i=0; i<nwin; i++) {
      tStart[i]=tStart[i]+firstbase;
    }
  }
  start.data=tStart;
  return start;
}

array<int> getEnd(int nsites, int firstbase, int block_size) {
  int nwin = ((nsites-firstbase)/block_size);
  array<int> end;
  if ( ((nsites-firstbase) % block_size)!=0) {
    end.x=nwin+1;
  } else {
    end.x=nwin;
  }
  int *tEnd= new int [nwin]; 
  for (int i=0; i<nwin; i++) {
    tEnd[i]=(i+1)*block_size-1;
  }
  // if there is rest
  if ( (nsites % block_size)!=0) {
    tEnd[nwin]=nsites-1;
    nwin=nwin+1;
  }
  // if you dont start from beginning
  if (firstbase>0) {
    for (int i=0; i<nwin; i++) {
      tEnd[i]=tEnd[i]+firstbase;
    }
  }
  end.data=tEnd;
  return end;
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

// read a file of prior into an array
array<double> readArray(const char *fname, int nInd, int isfold) {
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if(filesize==0){
    fprintf(stderr,"file:%s looks empty\n",fname);
    exit(0);
  }
  array<double> ret;
  if (isfold) {
    double *tmp = new double[nInd+1];
    char *buf = new char[filesize];
    fread(buf,sizeof(char),filesize,fp);
    tmp[0] = atof(strtok(buf,"\t \n"));
    for(int i=1;i<(nInd+1);i++)
      tmp[i] = atof(strtok(NULL,"\t \n"));
    fclose(fp);
    ret.x = nInd+1;
    ret.data = tmp;
   } else {
    double *tmp = new double[2*nInd+1];
    char *buf = new char[filesize];
    fread(buf,sizeof(char),filesize,fp);
    tmp[0] = atof(strtok(buf,"\t \n"));
    for(int i=1;i<(2*nInd+1);i++)
      tmp[i] = atof(strtok(NULL,"\t \n"));
    fclose(fp);
    ret.x = 2*nInd+1;
    ret.data = tmp;  
   }
  return ret;
}

matrix<double> readPrior12(const char *fname, int nrow, int ncol) {
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if(filesize==0){
    fprintf(stderr,"file:%s looks empty\n",fname);
    exit(0);
  }
  double *tmp = new double[nrow*ncol];
  char *buf = new char[filesize];
  fread(buf,sizeof(char),filesize,fp);
  tmp[0] = atof(strtok(buf,"\t \n"));
  for(int i=1;i<(nrow*ncol);i++)
    tmp[i] = atof(strtok(NULL,"\t \n"));
  fclose(fp);
  array<double> allvalues;
  allvalues.x = nrow*ncol;
  allvalues.data = tmp;
  int index=0;
  double **data = new double*[nrow];
  for(int i=0;i<nrow;i++){
    double *tmp2 = new double[ncol];
    for(int k=0;k<ncol;k++) {
      tmp2[k]=allvalues.data[index];
      index=index+1;
    }
    data[i]= tmp2;
  }
  matrix<double> ret;
  ret.x=nrow;
  ret.y=ncol;
  ret.data=data;
  delete [] allvalues.data;
  return ret;
}


// read a file of posterior probabilities into a matrix but only for a specific subsets of positions (0-based notation)
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
  int nsites = end-start+1;  
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

// read text output of fst values into an array (again for a subset of values)
array<double> readFSTsub(const char *fname, int nsites, int start, int end) {
  // nsites is how many sites you want
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if(filesize==0){
    fprintf(stderr,"file:%s looks empty\n",fname);
    exit(0);
  }
  double *tmp = new double[nsites*5];
  char *buf = new char[filesize];
  fread(buf,sizeof(char),filesize,fp);
  tmp[0] = atof(strtok(buf,"\t \n"));
  for(int i=1;i<(nsites*5);i++)
    tmp[i] = atof(strtok(NULL,"\t \n"));

  fclose(fp);
  array<double> allvalues;
  allvalues.x = nsites*5;
  allvalues.data = tmp;
  // get only column with fst for only the sites you want
  int index=(-2);
  int tt=0;
  double *tmp2 = new double[end-start+1];
  for(int i=start;i<=end;i++) {
    index=index+5;
    tmp2[tt] = allvalues.data[index];
    tt=tt+1;
  }
  delete [] allvalues.data;
  array<double> ret;
  ret.x = tt;
  ret.data = tmp2;
  return ret;
}

// write an array of doubles into a file
void writearray(array<double> &m,FILE *fp) {
  for(int i=0;i<m.x;i++)
    fprintf(fp,"%f\t",m.data[i]);
  fprintf(fp,"\n");
}

// write an array of ints into a file
void writearrayInt(array<int> &m,FILE *fp) {
  for(int i=0;i<m.x;i++)
    fprintf(fp,"%d\t",m.data[i]);
  fprintf(fp,"\n");
  
}

// write a matrix of doubles into a file
void writematrix(matrix<double> &m,FILE *fp) {
  for(int i=0;i<m.x;i++){
    for(int j=0;j<m.y;j++)
      fprintf(fp,"%f\t",m.data[i][j]);
    fprintf(fp,"\n");
  }
}

// write a matrix of ints into a file
void writematrixInt(matrix<int> &m,FILE *fp) {
  for(int i=0;i<m.x;i++){
    for(int j=0;j<m.y;j++)
      fprintf(fp,"%d\t",m.data[i][j]);
    fprintf(fp,"\n");
  }
}

// to append names
char *append(const char* a,const char *b){
  char *c =(char *) malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}

// print help
void info() {
  fprintf(stderr,"\t -> Required args:\n");
  fprintf(stderr,"\t\t-postfiles\tFiles with posterior probabilities\n");
  fprintf(stderr,"\t\t-outfile\tOutput file\n");
  fprintf(stderr,"\t\t-nind\tNumber of individuals per population\n");
  fprintf(stderr,"\t\t-nsites\tNumber sites (how many sites you want to investigate)\n");
  fprintf(stderr,"\t\t-verbose\tVerbose? [0]\n");
  fprintf(stderr,"\t\t-nsums\tnr of sums of moments? [0]\n");
  fprintf(stderr,"\t\t-firstbase\twhere to start in the file [0]\n");
}

/// 

// compute a and ab (using short-cut formulas) given allele frequencies and sample sizes
array<double> calcAB(int s1, int s2, int n1, int n2, int debug) {
  if (debug) fprintf(stderr, "\nInside calcAB");
  double alfa1 = 0.0, alfa2 = 0.0, p1 = 0.0, p2 = 0.0, a = 0.0, ab = 0.0;
  array<double> res; res.x = 2;
  double *values= new double [2];
  if (debug) fprintf(stderr, "%f\t%f\t%f\t%f\t%f\t%f\n", alfa1, alfa2, p1,p2,a,ab);
  p1 = static_cast<double>(s1) / (2*n1);
  p2 = static_cast<double>(s2) / (2*n2);
  if (debug) fprintf(stderr, "%f\t%f\n", p1, p2);
  alfa1 = 2*p1*(1-p1);  
  alfa2 = 2*p2*(1-p2);  
  if (debug) fprintf(stderr, "%f\t%f\n", alfa1, alfa2);
  a = ((p1-p2)*(p1-p2)) - (((n1+n2) * (n1 * alfa1 + n2 * alfa2)) / (4*n1*n2*(n1+n2-1)));
  ab = ((p1-p2)*(p1-p2)) + ( ( (4*n1*n2 - n1 - n2) * (n1 * alfa1 + n2 * alfa2) ) / ((4*n1*n2)*(n1+n2-1)) );
  values[0]=a;
  values[1]=ab;
  if (debug) fprintf(stderr, "%f\t%f\n", a, ab);
  res.data = values;
  return res;
}

// optimizing function gsl version
double myfunc_gsl (double lambda, void * params) {
  my_f_params * p = (my_f_params*) params;
  double k = (p->k);
  double z = (p->z);
  double F = (p->F);
  array<double> p1 = (p->p1);
  array<double> p2 = (p->p2);
  double somma=0.0;
  double dp=0.0;
  array<double> fmat;
  double mu = 0.0;
  mu = 1/lambda;
  for (int i=0; i<=k; i++) {
    for (int j=0; j<=z; j++) {
   	  fmat=calcAB(i, j, k/2, z/2, 0);
      if (isnan(fmat.data[0]/fmat.data[1])==0) {
       	dp= fabs( (static_cast<double>(i)/k ) - (static_cast<double>(j)/z ) );
        somma = somma + (fmat.data[0]/fmat.data[1])*p1.data[i]*p2.data[j]*gsl_cdf_exponential_Q(dp+0.0001,mu);
      }
      delete [] fmat.data;
    }
  }
  return fabs(somma-F);
}

// optimizing function if folded gsl version
double myfuncFold_gsl (double lambda, void * params) {
  my_f_params * p = (my_f_params*) params;
  double k = (p->k);
  double z = (p->z);
  double F = (p->F);
  array<double> p1 = (p->p1);
  array<double> p2 = (p->p2);
  double somma=0.0;
  double dp=0.0;
  array<double> fmat;
  double mu = 0.0;
  mu = 1/lambda;
  for (int i=0; i<=k; i++) {
    for (int j=0; j<=z; j++) {
   	  fmat=calcAB(i, j, k, z, 0);
      if (isnan(fmat.data[0]/fmat.data[1])==0) {
       	dp= fabs( (static_cast<double>(i)/(k*2) ) - (static_cast<double>(j)/(z*2) ) );
        somma = somma + (fmat.data[0]/fmat.data[1])*p1.data[i]*p2.data[j]*gsl_cdf_exponential_Q(dp+0.0001,mu);
      }
      delete [] fmat.data;
    }
  }
  return fabs(somma-F);
}

// optimizing function
double myfunc(double lambda, int k, int z, double F, array<double> p1, array<double> p2) {
  double somma=0.0;
  double dp=0.0;
  array<double> fmat;
  double mu = 0.0;
  mu = 1/lambda;
  for (int i=0; i<=k; i++) {
    for (int j=0; j<=z; j++) {
   	  fmat=calcAB(i, j, k/2, z/2, 0);
      if (isnan(fmat.data[0]/fmat.data[1])==0) {
       	dp= fabs( (static_cast<double>(i)/k ) - (static_cast<double>(j)/z ) );
        somma = somma + (fmat.data[0]/fmat.data[1])*p1.data[i]*p2.data[j]*gsl_cdf_exponential_Q(dp+0.0001,mu);
      }
      delete [] fmat.data;
    }
  }
  return fabs(somma-F);
}

// optimizing function if folded
double myfuncFold(double lambda, int k, int z, double F, array<double> p1, array<double> p2) {
  double somma=0.0;
  double dp=0.0;
  array<double> fmat;
  double mu = 0.0;
  mu = 1/lambda;
  for (int i=0; i<=k; i++) {
    for (int j=0; j<=z; j++) {
  	  fmat=calcAB(i, j, k, z, 0);
      if (isnan(fmat.data[0]/fmat.data[1])==0) {
       	dp= fabs( (static_cast<double>(i)/(k*2) ) - (static_cast<double>(j)/(z*2) ) );
        somma = somma + (fmat.data[0]/fmat.data[1])*p1.data[i]*p2.data[j]*gsl_cdf_exponential_Q(dp+0.0001,mu);
      }
      delete [] fmat.data;
    }
  }
  return fabs(somma-F);
}

// compute lambdas
array<double> getLambdas(array<double> myfst, array<double> prob1, array<double> prob2, int K, int verbose, int isfold) {

//    fprintf(stderr, "\ninside get lamb; fst are %d \n", myfst.x);
//    writearray(myfst, stdout);

//    fprintf(stderr, "\nprob1");
//    writearray(prob1, stdout);

//    fprintf(stderr, "\nprob2");
//    writearray(prob2, stdout);

    array <double> lambdas;
    lambdas.x=myfst.x;
    double *ltmp = new double[myfst.x];
    for (int i=0; i<myfst.x; i++) {
      ltmp[i]=0.0;
    }
    lambdas.data=ltmp;
    

    int nind1 = (prob1.x-1)/2;
    int nind2 = (prob2.x-1)/2;

    // initialize values for optimization
    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    double m = 1; // starting values
    double a = 0.0001, b = 30; // interval
    if (verbose) fprintf(stderr, "\nOptimizing to get lambdas...");
    // my grid
    double grid[35] = {0.00015,0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 11.0, 14.0, 17.0, 20.0, 22.0, 25.0, 29.0};
    double tfm[35]= {}; // f(grid)
    double minu=100.00; // init
    int iminu=999; // init

    int alreadydone=(-1); // if -1 no duplicates, otherwise it is the index of last duplicate

    // for each site
    if (K==0) {

      for (int i=0; i<myfst.x; i++) {

        alreadydone=-1;
        //fprintf(stderr, "\n site %d : %d", i, alreadydone);
        if(i>0) {
          for (int j=0; j<i; j++) {
            if (myfst.data[i]==myfst.data[j]) {
              alreadydone=j;
              //fprintf(stderr, "\ncheck site %d at %f ? %f : %d", i, myfst.data[i], myfst.data[j], alreadydone);
              break;
            }
          }
        }

        if (alreadydone==(-1)) {

//      fprintf(stderr, "\n!!! site %d out of %d with fst %f", i, myfst.x, myfst.data[i]);

        // re-init
        a=0.0001;
        b=30;
        m=1;
        // check if fst out of range
        if (myfst.data[i]<(0.0001))  {
          m = 30; }
        else if (myfst.data[i]>0.3) {
          m = 0.0001; }
        else {
         // a mix of grid-search and brent algorithm
         // first a grid-search for 
         for (int j=0; j<35; j++) {
//           printf("\n %f", grid[j]);
           if (isfold) {
             tfm[j]=myfuncFold(grid[j], nind1, nind2, myfst.data[i], prob1, prob2);
           } else {
             tfm[j]=myfunc(grid[j], nind1*2, nind2*2, myfst.data[i], prob2, prob2);
           }
//           printf("\n %f:%f", grid[j], tfm[j]);
         }
       // what's the minimum?
       minu=100.00;
       iminu=999;
       for (int j=0; j<35; j++) {
         if (tfm[j]<minu) {
           minu=tfm[j];
           iminu=j;
         }
       }
//       printf("\n %f:%d", minu, iminu);
       // reset a,b,m
       m = grid[iminu]; // starting values
       // if it is one of the bounds then set miniumu as 0.0001 or 30
       if (iminu==0) {
         m=0.0001; }
       else if (iminu==34) {
         m=30; 
       } else {
        // do brent
        a=grid[iminu-1];
        b=grid[iminu+1];
        my_f_params params = { nind1*2, nind2*2, myfst.data[i], prob1, prob2};
        gsl_function F;
        if (isfold) {
          F.function = &myfuncFold_gsl;
        } else {
          F.function = &myfunc_gsl;        	
        }
        F.params = &params;

        T = gsl_min_fminimizer_brent;
        s = gsl_min_fminimizer_alloc (T);
        gsl_min_fminimizer_set (s, &F, m, a, b);

        do {
          iter++;
          status = gsl_min_fminimizer_iterate (s);

          if (status == GSL_SUCCESS) {
//            printf ("Converged:\n");
//            printf ("\t %5d [%.7f, %.7f] %.7f %.7f\n", iter, a, b, m, b - a);
          }
        } // end do
        while (status == GSL_CONTINUE && iter < max_iter);

        gsl_min_fminimizer_free (s);

      } // end if do brent

    } // end if it is not a bound

    lambdas.data[i]=m;

    //if (i==0) fprintf(stderr, "\n ended at: %f and ltmp is %f",m, lambdas.data[i]);

    } else { // end if not already done

    lambdas.data[i]=lambdas.data[alreadydone];
    //fprintf(stderr, ".. lam %f is %f", lambdas.data[i], lambdas.data[alreadydone]);

    }

    //sleep(2);

    } // end for every sites

    //lambdas.x=myfst.x;
    //lambdas.data = ltmp;

  } else { // use fixed lambdas based on 1/(K*FST)
   for (int i=0; i<myfst.x; i++) {
 //    fprintf(stdout, "\n %d %f ", i, myfst.data[i]);
     lambdas.data[i]=1/(K*myfst.data[i]);

   }

//   lambdas.x=myfst.x;
//   lambdas.data = ltmp;

  }
  // now I have my lambdas
  // now for each site (post probs per site) compute the corrected post probs
  return lambdas;
} 

// compute a and ab estimate for all possible combinations of sample size, then weight by their prob, but no correction from a first guess of fst
void computeVarRey(matrix<double> &m1, matrix<double> &m2, int verbose, FILE *fname, int nsums) {
  // m1 and m2 are post probs
  /// DEFINITION (DECLARATION AND INITIALIZATION)
  if (verbose==2) fprintf(stderr, "\nInside computeVarRey");
  int n1 = 0, n2 = 0; // sample sizes
  int nsites = 0; // nsites
  // estimates nsites and pop sizes from matrices
  n1 = (m1.y-1)/2; // nind1
  n2 = (m2.y-1)/2; // nind2
  nsites = m1.x; // nsites from file
  double VAR = 0.0, COVAR = 0.0, FACT = 0.0;
  matrix<double> A;
  matrix<double> AB;
  array<double> temp;
  for (int s=0; s<nsites; s++) {
    if (verbose==2) fprintf(stderr, "\t s %d", s);
    // for each possible value of freq 1 and freq 2 compute the FST, so compute A, AB, VAR, COVAR (see Price paper for its meaning)
    A.x=AB.x=(n1*2)+1;
    A.y=AB.y=(n2*2)+1;
    // FIRST CYCLE: get expected A and AB and retain matrices of A and AB
    double **dataA = new double*[(n1*2)+1];
    double **dataAB = new double*[(n1*2)+1];
    // get also the probability of site being variable
    double pvar = 0.0;
    pvar = 1 - m1.data[s][0]*m2.data[s][0] - m1.data[s][2*n1]*m2.data[s][2*n2];
    if (verbose==2) fprintf(stderr, "\t first cycle");
    double EA = 0.0, EAB = 0.0;
    for (int i=0; i<(2*n1+1); i++) {
      if (verbose==2) fprintf(stderr, "\ti%d",i);
      double *bufA = new double[(n2*2)+1];
      double *bufAB = new double[(n2*2)+1];
      for (int j=0; j<(2*n2+1); j++) {
       if (verbose==2) fprintf(stderr, "\tj%d",j);
        temp = calcAB(i, j, n1, n2, 0);
        bufA[j]=temp.data[0];
        bufAB[j]=temp.data[1];
        EA = EA + temp.data[0]*m1.data[s][i]*m2.data[s][j];
        EAB = EAB + temp.data[1]*m1.data[s][i]*m2.data[s][j];
        delete [] temp.data;
      } // end for in j
      dataA[i]=bufA;
      dataAB[i]=bufAB;
    } // end for in i
    A.data=dataA;
    AB.data=dataAB;
    // SECOND CYCLE: get VAR and COVAR, and then the correcting FACTor, according to number of sums to retain
    VAR = 0.0, COVAR = 0.0, FACT = 0.0;
    if ((verbose==4) & (s==0)) fprintf(stderr, "\t second cycle %d", nsums);
    for (int q=1; q<=nsums; q++) {
      VAR = 0.0; COVAR = 0.0;
      for (int i=0; i<(2*n1+1); i++) {
        for (int j=0; j<(2*n2+1); j++) {
          VAR = VAR + pow((AB.data[i][j]-EAB), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j];
          COVAR = COVAR + pow((AB.data[i][j]-EAB)*(A.data[i][j]-EA), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j];
        } // end for in j
      } // end for in i
      FACT = FACT + pow( (-1.0), static_cast <double> (q)) *( (EA*VAR + COVAR) / pow(EAB, static_cast <double> (q+1)));
      if ((verbose==4) & (s==0)) fprintf(stderr, "\n q %d v %f c %f f %f",q,VAR,COVAR,FACT);
    } // end for in nsums
    // print results
    fprintf(fname, "%f\t%f\t%f\t%f\t%f\n", EA, EAB, FACT, (EA/EAB)+FACT, pvar);
    cleanup(A);
    cleanup(AB);
  } // end for s in nsites
} // end


// compute a and ab estimate for all possible combinations of sample size, then weight by their prob, , but no correction from a first guess of fst, here if folded!
void computeVarReyFold(matrix<double> &m1, matrix<double> &m2, int verbose, FILE *fname, int nsums) {
  // m1 and m2 are post probs
  /// DEFINITION (DECLARATION AND INITIALIZATION)
  if (verbose==2) fprintf(stderr, "\nInside computeVarRey");
  int n1 = 0, n2 = 0; // sample sizes
  int nsites = 0; // nsites
  // estimates nsites and pop sizes from matrices
  n1 = (m1.y-1); // nind1
  n2 = (m2.y-1); // nind2
  nsites = m1.x; // nsites from file
  double VAR = 0.0, COVAR = 0.0, FACT = 0.0;
  matrix<double> A;
  matrix<double> AB;
  array<double> temp;
  for (int s=0; s<nsites; s++) {
    if (verbose==2) fprintf(stderr, "\t s %d", s);
    // for each possible value of freq 1 and freq 2 compute the FST, so compute A, AB, VAR, COVAR (see Price paper for its meaning)
    A.x=AB.x=(n1)+1;
    A.y=AB.y=(n2)+1;
    // FIRST CYCLE: get expected A and AB and retain matrices of A and AB
    double **dataA = new double*[(n1)+1];
    double **dataAB = new double*[(n1)+1];
    // get also the probability of site being variable
    double pvar = 0.0;
    pvar = 1 - m1.data[s][0]*m2.data[s][0] - m1.data[s][2*n1]*m2.data[s][2*n2];
    if (verbose==2) fprintf(stderr, "\t first cycle");
    double EA = 0.0, EAB = 0.0;
    for (int i=0; i<(n1+1); i++) {
      if (verbose==2) fprintf(stderr, "\ti%d",i);
      double *bufA = new double[(n2)+1];
      double *bufAB = new double[(n2)+1];
      for (int j=0; j<(n2+1); j++) {
       if (verbose==2) fprintf(stderr, "\tj%d",j);
        temp = calcAB(i, j, n1, n2, 0);
        bufA[j]=temp.data[0];
        bufAB[j]=temp.data[1];
        EA = EA + temp.data[0]*m1.data[s][i]*m2.data[s][j];
        EAB = EAB + temp.data[1]*m1.data[s][i]*m2.data[s][j];
        delete [] temp.data;
      } // end for in j
      dataA[i]=bufA;
      dataAB[i]=bufAB;
    } // end for in i
    A.data=dataA;
    AB.data=dataAB;
    // SECOND CYCLE: get VAR and COVAR, and then the correcting FACTor, according to number of sums to retain
    VAR = 0.0, COVAR = 0.0, FACT = 0.0;
    if ((verbose==4) & (s==0)) fprintf(stderr, "\t second cycle %d", nsums);
    for (int q=1; q<=nsums; q++) {
      VAR = 0.0; COVAR = 0.0;
      for (int i=0; i<(n1+1); i++) {
        for (int j=0; j<(n2+1); j++) {
          VAR = VAR + pow((AB.data[i][j]-EAB), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j];
          COVAR = COVAR + pow((AB.data[i][j]-EAB)*(A.data[i][j]-EA), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j];
        } // end for in j
      } // end for in i
      FACT = FACT + pow( (-1.0), static_cast <double> (q)) *( (EA*VAR + COVAR) / pow(EAB, static_cast <double> (q+1)));
      if ((verbose==4) & (s==0)) fprintf(stderr, "\n q %d v %f c %f f %f",q,VAR,COVAR,FACT);
    } // end for in nsums
    // print results
    fprintf(fname, "%f\t%f\t%f\t%f\t%f\n", EA, EAB, FACT, (EA/EAB)+FACT, pvar);
    cleanup(A);
    cleanup(AB);
  } // end for s in nsites
} // end

// compute weights normalizing over posterior probabilities
matrix<double> getWeights(int n1, int n2, matrix<double> post1, matrix<double> post2, int s, double herelam) {
  int k = n1*2+1;
  int z = n2*2+1;
  double dp=0.0;
  double mu=1/herelam;
  matrix<double> weights;
  // compute weights from exponential and from the computed lambda
  double **data = new double*[k];
  for(int i=0;i<k;i++) {
    double *tmp = new double[z];
    for (int j=0; j<z; j++) {
      dp= fabs( (static_cast<double>(i)/(k-1) ) - (static_cast<double>(j)/(z-1) ) );
      tmp[j]= gsl_cdf_exponential_Q(dp+0.0001, mu);
    }
    data[i]=tmp;
  }
  weights.x=k;
  weights.y=z;
  weights.data=data;
  // normalize
  double somma=0.0;
  for(int i=0;i<k;i++) {
    for (int j=0; j<z; j++) {
      somma=somma+post1.data[s][i]*post2.data[s][j]*weights.data[i][j];
    }
  }
  for(int i=0;i<k;i++) {
    for (int j=0; j<z; j++) {
      weights.data[i][j]=weights.data[i][j]/somma;
    }
  }
  // now they sum up to 1
  return weights;
}

// compute weights normalizing over posterior probabilities if folded
matrix<double> getWeightsFold(int n1, int n2, matrix<double> post1, matrix<double> post2, int s, double herelam) {
  int k = n1+1;
  int z = n2+1;
  double dp=0.0;
  double mu=1/herelam;
  matrix<double> weights;
  // compute weights from exponential and from the computed lambda
  double **data = new double*[k];
  for(int i=0;i<k;i++) {
    double *tmp = new double[z];
    for (int j=0; j<z; j++) {
      dp= fabs( (static_cast<double>(i)/(k-1) ) - (static_cast<double>(j)/(z-1) ) );
      tmp[j]= gsl_cdf_exponential_Q(dp+0.0001, mu);
    }
    data[i]=tmp;
  }
  weights.x=k;
  weights.y=z;
  weights.data=data;
  // normalize
  double somma=0.0;
  for(int i=0;i<k;i++) {
    for (int j=0; j<z; j++) {
      somma=somma+post1.data[s][i]*post2.data[s][j]*weights.data[i][j];
    }
  }
  for(int i=0;i<k;i++) {
    for (int j=0; j<z; j++) {
      weights.data[i][j]=weights.data[i][j]/somma;
    }
  }
  // now they sum up to 1
  return weights;
}

// compute a and ab estimate for all possible combinations of sample size, then weight by their prob, and correction from a first guess of fst
void computeVarRey2(matrix<double> &m1, matrix<double> &m2, int verbose, FILE *fname, int nsums, array<double> lambdas) {
  // m1 and m2 are post probs
  /// DEFINITION (DECLARATION AND INITIALIZATION)
  if (verbose==20) fprintf(stderr, "\nInside computeVarRey");
  int n1 = 0, n2 = 0; // sample sizes
  int nsites = 0; // nsites
  // estimates nsites and pop sizes from matrices
  n1 = (m1.y-1)/2; // nind1
  n2 = (m2.y-1)/2; // nind2
  nsites = m1.x; // nsites from file
  double VAR = 0.0, COVAR = 0.0, FACT = 0.0;
  matrix<double> A;
  matrix<double> AB;
  matrix <double> wpp;
  if (verbose==20) fprintf(stderr, ";init...");
  array<double> temp;
  for (int s=0; s<nsites; s++) {
    // for each possible value of freq 1 and freq 2 compute the FST, so compute A, AB, VAR, COVAR (see Price paper for its meaning)
    A.x=AB.x=(n1*2)+1;
    A.y=AB.y=(n2*2)+1;
    // FIRST CYCLE: get expected A and AB and retain matrices of A and AB
    double **dataA = new double*[(n1*2)+1];
    double **dataAB = new double*[(n1*2)+1];
    // get the weights for each site
    wpp = getWeights(n1, n2, m1, m2, s, lambdas.data[s]);
    // get also the probability of site being variable
    double pvar = 0.0;
    pvar = 1 - m1.data[s][0]*m2.data[s][0]*wpp.data[0][0] - m1.data[s][2*n1]*m2.data[s][2*n2]*wpp.data[2*n1][2*n2];
    if (verbose==20) fprintf(stderr, "\t first cycle");
    double EA = 0.0, EAB = 0.0;
    for (int i=0; i<(2*n1+1); i++) {
      if (verbose==20) fprintf(stderr, "\ti%d",i);
      double *bufA = new double[(n2*2)+1];
      double *bufAB = new double[(n2*2)+1];
      for (int j=0; j<(2*n2+1); j++) {
       if (verbose==20) fprintf(stderr, "\tj%d",j);
        temp = calcAB(i, j, n1, n2, 0);
        if (verbose==20) fprintf(stderr, "\t a ab %f %f",temp.data[0],temp.data[1]);
        bufA[j]=temp.data[0];
        bufAB[j]=temp.data[1];
        EA = EA + temp.data[0]*m1.data[s][i]*m2.data[s][j]*wpp.data[i][j];
        EAB = EAB + temp.data[1]*m1.data[s][i]*m2.data[s][j]*wpp.data[i][j];
        delete [] temp.data;
      } // end for in j
      dataA[i]=bufA;
      dataAB[i]=bufAB;
    } // end for in i
    A.data=dataA;
    AB.data=dataAB;
    // SECOND CYCLE: get VAR and COVAR, and then the correcting FACTor, according to number of sums to retain
    VAR = 0.0, COVAR = 0.0, FACT = 0.0;
    if ((verbose==4) & (s==0)) fprintf(stderr, "\t second cycle %d", nsums);
    for (int q=1; q<=nsums; q++) {
      VAR = 0.0; COVAR = 0.0;
      for (int i=0; i<(2*n1+1); i++) {
        for (int j=0; j<(2*n2+1); j++) {
          VAR = VAR + pow((AB.data[i][j]-EAB), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j]*wpp.data[i][j];
          COVAR = COVAR + pow((AB.data[i][j]-EAB)*(A.data[i][j]-EA), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j]*wpp.data[i][j];
        } // end for in j
      } // end for in i
      FACT = FACT + pow( (-1.0), static_cast <double> (q)) *( (EA*VAR + COVAR) / pow(EAB, static_cast <double> (q+1)));
      if ((verbose==4) & (s==0)) fprintf(stderr, "\n q %d v %f c %f f %f",q,VAR,COVAR,FACT);
    } // end for in nsums
    // print results
    fprintf(fname, "%f\t%f\t%f\t%f\t%f\n", EA, EAB, FACT, (EA/EAB)+FACT, pvar);
    cleanup(A);
    cleanup(AB);
    cleanup(wpp);
  } // end for s in nsites
} // end

// compute a and ab estimate for all possible combinations of sample size, then weight by their prob, and correction from a first guess of fst, if folded
void computeVarRey2Fold(matrix<double> &m1, matrix<double> &m2, int verbose, FILE *fname, int nsums, array<double> lambdas) {
  // m1 and m2 are post probs
  /// DEFINITION (DECLARATION AND INITIALIZATION)
  if (verbose==20) fprintf(stderr, "\nInside computeVarRey");
  int n1 = 0, n2 = 0; // sample sizes
  int nsites = 0; // nsites
  // estimates nsites and pop sizes from matrices
  n1 = (m1.y)/2; // nind1
  n2 = (m2.y)/2; // nind2
  nsites = m1.x; // nsites from file
  double VAR = 0.0, COVAR = 0.0, FACT = 0.0;
  matrix<double> A;
  matrix<double> AB;
  matrix <double> wpp;
  if (verbose==20) fprintf(stderr, ";init...");
  array<double> temp;
  for (int s=0; s<nsites; s++) {
    if (verbose==20) fprintf(stderr, "\t s %d", s);
    // for each possible value of freq 1 and freq 2 compute the FST, so compute A, AB, VAR, COVAR (see Price paper for its meaning)
    A.x=AB.x=(n1)+1;
    A.y=AB.y=(n2)+1;
    // FIRST CYCLE: get expected A and AB and retain matrices of A and AB
    double **dataA = new double*[(n1)+1];
    double **dataAB = new double*[(n1)+1];
    // get the weights for each site
    wpp = getWeightsFold(n1, n2, m1, m2, s, lambdas.data[s]);
    // get also the probability of site being variable
    double pvar = 0.0;
    pvar = 1 - m1.data[s][0]*m2.data[s][0]*wpp.data[0][0] - m1.data[s][2*n1]*m2.data[s][2*n2]*wpp.data[2*n1][2*n2];
    if (verbose==20) fprintf(stderr, "\t first cycle");
    double EA = 0.0, EAB = 0.0;
    for (int i=0; i<(n1+1); i++) {
      if (verbose==20) fprintf(stderr, "\ti%d",i);
      double *bufA = new double[(n2)+1];
      double *bufAB = new double[(n2)+1];
      for (int j=0; j<(n2+1); j++) {
       if (verbose==20) fprintf(stderr, "\tj%d",j);
        temp = calcAB(i, j, n1, n2, 0);
        if (verbose==20) fprintf(stderr, "\t a ab %f %f",temp.data[0],temp.data[1]);
        bufA[j]=temp.data[0];
        bufAB[j]=temp.data[1];
        EA = EA + temp.data[0]*m1.data[s][i]*m2.data[s][j]*wpp.data[i][j];
        EAB = EAB + temp.data[1]*m1.data[s][i]*m2.data[s][j]*wpp.data[i][j];
        delete [] temp.data;
      } // end for in j
      dataA[i]=bufA;
      dataAB[i]=bufAB;
    } // end for in i
    A.data=dataA;
    AB.data=dataAB;
    // SECOND CYCLE: get VAR and COVAR, and then the correcting FACTor, according to number of sums to retain
    VAR = 0.0, COVAR = 0.0, FACT = 0.0;
    if ((verbose==4) & (s==0)) fprintf(stderr, "\t second cycle %d", nsums);
    for (int q=1; q<=nsums; q++) {
      VAR = 0.0; COVAR = 0.0;
      for (int i=0; i<(n1+1); i++) {
        for (int j=0; j<(n2+1); j++) {
          VAR = VAR + pow((AB.data[i][j]-EAB), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j]*wpp.data[i][j];
          COVAR = COVAR + pow((AB.data[i][j]-EAB)*(A.data[i][j]-EA), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j]*wpp.data[i][j];
        } // end for in j
      } // end for in i
      FACT = FACT + pow( (-1.0), static_cast <double> (q)) *( (EA*VAR + COVAR) / pow(EAB, static_cast <double> (q+1)));
      if ((verbose==4) & (s==0)) fprintf(stderr, "\n q %d v %f c %f f %f",q,VAR,COVAR,FACT);
    } // end for in nsums
    // print results
    fprintf(fname, "%f\t%f\t%f\t%f\t%f\n", EA, EAB, FACT, (EA/EAB)+FACT, pvar);
    cleanup(A);
    cleanup(AB);
    cleanup(wpp);
  } // end for s in nsites
} // end

// compute a and ab estimate for all possible combinations of sample size, then weight by their prob12 computed from post1, pos12 and prior12 (normalize it)
void computeVarRey12New(matrix<double> &m1, matrix<double> &m2, int verbose, FILE *fname, int nsums, matrix<double> &p12) {

  int n1 = 0, n2 = 0; // sample sizes
  int nsites = 0; // nsites

  // estimates nsites and pop sizes from matrices
  n1 = (m1.y-1)/2; // nind1
  n2 = (m2.y-1)/2; // nind2
  nsites = m1.x; // nsites from file

  double VAR = 0.0, COVAR = 0.0, FACT = 0.0;

  for (int s=0; s<nsites; s++) {

    // m1 and m2 are post probs, p12 is the 2d sfs
    matrix<double> m12;
    m12.x=m1.y;
    m12.y=m2.y;
    double **ddata = new double*[m12.x];
      for(int i=0;i<m12.x;i++){
        double *dtmp = new double[m12.y];
        for(int k=0;k<m12.y;k++) {
          dtmp[k]=0.0;
        }
        ddata[i]= dtmp;
      }
    m12.data=ddata;

    // multiply
    for (int j=0; j<m12.x; j++) {
      for (int i=0;i<m12.y; i++) {
        m12.data[j][i] = m1.data[s][j]*m2.data[s][i]*p12.data[j][i];
      }
    }
    
    // for each possible value of freq 1 and freq 2 compute the FST, so compute A, AB, VAR, COVAR (see Price paper for its meaning)
    matrix<double> A;
    matrix<double> AB;

    A.x=AB.x=(n1*2)+1;
    A.y=AB.y=(n2*2)+1;

    // FIRST CYCLE: get expected A and AB and retain matrices of A and AB
    double **dataA = new double*[(n1*2)+1];
    double **dataAB = new double*[(n1*2)+1];

    // get also the probability of site being variable
    double pvar = 0.0;
    pvar = 1 - m12.data[0][0] - m12.data[2*n1][2*n2];

    if (verbose==2) fprintf(stderr, "\t first cycle");

    double EA = 0.0, EAB = 0.0;

    for (int i=0; i<(2*n1+1); i++) {

      if (verbose==2) fprintf(stderr, "\ti%d",i);

      double *bufA = new double[(n2*2)+1];
      double *bufAB = new double[(n2*2)+1];

      for (int j=0; j<(2*n2+1); j++) {

        array<double> temp;
        temp = calcAB(i, j, n1, n2, 0);
        bufA[j]=temp.data[0];
        bufAB[j]=temp.data[1];

        EA = EA + temp.data[0]*m12.data[i][j];
        EAB = EAB + temp.data[1]*m12.data[i][j];

        delete [] temp.data;

      } // end for in j

      dataA[i]=bufA;
      dataAB[i]=bufAB;

    } // end for in i

    A.data=dataA;
    AB.data=dataAB;

    // SECOND CYCLE: get VAR and COVAR, and then the correcting FACTor, according to number of sums to retain
    VAR = 0.0, COVAR = 0.0, FACT = 0.0;

    if ((verbose==4) & (s==0)) fprintf(stderr, "\t second cycle %d", nsums);

    for (int q=1; q<=nsums; q++) {

      VAR = 0.0; COVAR = 0.0;
      for (int i=0; i<(2*n1+1); i++) {
        for (int j=0; j<(2*n2+1); j++) {

          VAR = VAR + pow((AB.data[i][j]-EAB), static_cast <double> (q) )*m12.data[i][j];

          COVAR = COVAR + pow((AB.data[i][j]-EAB)*(A.data[i][j]-EA), static_cast <double> (q) ) *m12.data[i][j];

        } // end for in j
      } // end for in i

      FACT = FACT + pow( (-1.0), static_cast <double> (q)) *( (EA*VAR + COVAR) / pow(EAB, static_cast <double> (q+1)));

      if ((verbose==4) & (s==0)) fprintf(stderr, "\n q %d v %f c %f f %f",q,VAR,COVAR,FACT);

    } // end for in nsums

    // print results
    fprintf(fname, "%f\t%f\t%f\t%f\t%f\n", EA, EAB, FACT, (EA/EAB)+FACT, pvar);

    cleanup(A);
    cleanup(AB);
    cleanup(m12);

  } // end for s in nsites



} // end


