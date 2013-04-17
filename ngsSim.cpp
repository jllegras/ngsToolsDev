// program to simulate NGS data for up 3 populations

// author Rasmus Nielsen with contributions from Thorfinn Kornelliussen and Matteo Fumagalli

// to compile: g++ ngsSim.cpp -o ngsSim -lm -lz -O0

/*Note 0=A, 1=C, 2=G, 3=T*/

#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <cstring>
#include <zlib.h>
#include <sys/stat.h>
#define PI 3.141592654 // for Poisson distribution

#include "rbeta.cpp" // random generator of values beta-distributed

// define some constants or global variables

double minfreq, myConst;

int dumpBinary = 1; // output as a binary file?
int simpleRand = 0; // which uniform function?

// declare functions

int static SetSeed(int seed) {
  int static z_rndu;
  z_rndu = 170*(seed%178) + 137;
  return z_rndu;
}

/*U(0,1): AS 183: Appl. Stat. 31:188-190
Wichmann BA & Hill ID.  1982.  An efficient and portable pseudo-random number generator.  Appl. Stat. 31:188-190 x, y, z are any numbers in the range 1-30000.  Integer operation up to 30323 required.
Suggested to me by Z. Yang who also provided me with the source code used here. */
double uniform_rasmus(int seed=0) {
  static int x_rndu=11, y_rndu=23;
  static int z_rndu;
  if (seed!=0)
    z_rndu = SetSeed(seed);
  else
    z_rndu=137;
  double r;
  x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
  y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
  z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
  if (x_rndu<0) x_rndu+=30269;
  if (y_rndu<0) y_rndu+=30307;
  if (z_rndu<0) z_rndu+=30323;
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
  return (r-(int)r);
}

double uniform_thorfinn() {
  return((double)rand() / (double)RAND_MAX);
}

// which uniform function to use
double uniform(int seed=0) {
  double sampled;
  if(seed==0)
    sampled = uniform_thorfinn();
  else
    sampled = uniform_rasmus(seed);
  return sampled;
}

// function to check for file existence, using stat.h.
int fexists(const char* str){ ///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

// functions to get gz files
gzFile getGz(const char*fname,const char* mode){
  if(fexists(fname)){
    fprintf(stderr,"\t->File exists: %s exiting...\n",fname);
    exit(0);
  }
  gzFile fp;
  if(NULL==(fp=gzopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening gzFile handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

FILE *getFile(const char*fname,const char *mode){
  if(strncmp(mode,"r",1) && fexists(fname)){
    fprintf(stderr,"\t-> File exists: %s exiting... Terminate\n",fname);
    exit(0);
  }
  FILE *fp=NULL;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening File handle for file:%s\n",fname);
    exit(0);
  }
  return fp;
}

// gamma distribution
double gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947,-86.50532032942, 24.01409824083,-1.231739572450, 0.1208650973866e-2,-0.5395239385e-5};

  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -=(x+0.5)*log(tmp);
  ser=1.00000000019015;
  for (j=0; j<5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

// Poisson distribution
double Poisson(double xm) {
  double gammln(double xx);
  static double sq,alxm,g,oldm=(-1.0);
  double em, t, y;
  
  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em=-1;
    t=1.0;
    do {
      ++em;
      t *=uniform();
    } while (t>g);
  } 
  else {
    if (xm!=oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
	y=tan(PI*uniform());
	em=sq*y+xm;
      } while (em< 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (uniform()>t);
  }
  return em;
}

// generate sequencing errors by picking a different base
int basepick_with_errors(double errate, int inbase) {
  // from error rate and a starting base
  int outbase;
  if (uniform()<errate){ // if error occurs then change base
    while ((outbase=(floor(4*uniform()))) == inbase); // then take a different base
    return outbase;
  }
  else return inbase;
}

// pick a genotype, but with errors for each allele of the genotype[2]
int pick_a_base(double errate, int genotype[2])	{
  if (uniform()<0.5) 
    return basepick_with_errors(errate, genotype[0]);
  else 
    return basepick_with_errors(errate, genotype[1]);
}

// calculate the likelihood for each genotype configuration
void calclike(int base, double errate, double *like)	{
  /*0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
  double  prob;
  int i, j, k=0;
  for (i=0; i<4; i++){  // all 4 possible alleles for 1st base of genotype
    for (j=i; j<4; j++){ // all 4 possible alleles for 2nd base of genotype 
      if (base==i) 
	prob = (1.0-errate)/2.0; // if right base
      else 
	prob = errate/6.0; // if wrong base
      if (base==j) 
	prob = prob +(1.0-errate)/2.0;
      else 
	prob = prob + errate/6.0;
       if (prob <= 0.0) 
	like[k] = -1000000.0;
      else 
	like[k] = like[k] + log(prob); // add all log(prob) as likelihood
      k++;
    }
  }
}

// translate from int code to letters
char int_to_base(int b) {
  if (b==0) return 'A';
  else if (b==1) return 'C';
  else if (b==2) return 'G';
  else return 'T';
}

// compute and print results into files
int print_ind_site(double errate, double meandepth, int genotype[2], gzFile resultfile, gzFile glffile, FILE *fname) {
  
  int i, b, numreads;
  double like[10];
  
  int debug = 0;
  
  int ireads[4]; // sum of 0 1 2 3 reads for each individual
  for (int j=0; j<4; j++)
    ireads[j]=0;
  
  numreads = Poisson(meandepth);  // mumber of reads, poisson distributed with lambda=meandepth
  char res[numreads]; // char alleles for all the reads

  for (i=0; i<10; i++) like[i] = 0.0; // initialize GL values for all 10 possible genotypes
  
  // compute like
  for (i=0; i<numreads; i++) {
    b = pick_a_base(errate, genotype); // pick a base including errors

    // add ireads
    for (int j=0; j<4; j++)
     if (j==b) ireads[j]++;

    // debugging
    if (debug) {
      for (int j=0; j<4; j++)
	fprintf(stderr, "\nn %d count reads j: %d is : %d", i, j, ireads[j]);
    }  
       
    res[i] = int_to_base(b); 
    calclike(b, errate, like); // compute likelihood
  }

  // set maximum likelihoos to 0 and scale the rest
  int max_pos=-1;
  float max_val = -99.0;
  for (i=0; i<10; i++) {
    if (like[i]>max_val) {
      max_val=like[i];
      max_pos=i;
    }
  }
  for (i=0; i<10; i++) like[i] = like[i]-max_val;
  like[max_pos]=0;

  fprintf(fname, "%d\t%d\t%d\t%d\n", ireads[0], ireads[1], ireads[2], ireads[3]);

  // write into files
  if(dumpBinary) { // if binary 
    gzwrite(resultfile, res, numreads);
    gzwrite(glffile, like, sizeof(double)*10); // print likelihoods values
  } else { // text
    fprintf(stderr,"textout not implemented\n"); // not yet...
    exit(0);
  }

  return numreads; // output is the number of reads

}

// compute and print results into files when you have 2 populations
int print_ind_site2(double errate, double meandepth, int genotype[2], gzFile resultfile, gzFile glffile, FILE *fname, gzFile resultfile2, gzFile glffile2, FILE *fname2, int first) {
  int i, b, numreads;
  double like[10];
  int ireads[4]; // sum of 0 1 2 3 reads for each individual
  for (int j=0; j<4; j++)
    ireads[j]=0;
  numreads = Poisson(meandepth);  // mumber of reads, poisson distributed with lambda=meandepth
  char res[numreads]; // char alleles for all the reads
  for (i=0; i<10; i++) like[i] = 0.0; // initialize GL values for all 10 possible genotypes
  // compute like
  for (i=0; i<numreads; i++) {
    b = pick_a_base(errate, genotype); // pick a base including errors
    // add ireads
    for (int j=0; j<4; j++)
     if (j==b) ireads[j]++;
    res[i] = int_to_base(b); 
    calclike(b, errate, like); // compute likelihood
  }
  fprintf(fname, "%d\t%d\t%d\t%d\n", ireads[0], ireads[1], ireads[2], ireads[3]);
  fprintf(fname2, "%d\t%d\t%d\t%d\n", ireads[0], ireads[1], ireads[2], ireads[3]);
  // write into files  
  if(dumpBinary) { // if binary
   if (first==1) {
     char sep[1]={'\n'};
     gzwrite(resultfile, sep, 1);
   }
    gzwrite(resultfile, res, numreads);
    gzwrite(glffile, like, sizeof(double)*10); // print likelihoods values
    gzwrite(resultfile2, res, numreads);
    gzwrite(glffile2, like, sizeof(double)*10); // print likelihoods values
  } else { // text
    fprintf(stderr,"textout not implemented\n"); // not yet...
    exit(0);
  }
  return numreads; // output is the number of reads
}

// randomly pick a base from prior basefreq
int pick_base_from_prior(double basefreq[4])	{
  int i = 0;
  double U, p;
  
  U = uniform();
  p = basefreq[0];
  while (U>p)
    {
      i++;
      p = p + basefreq[i];
    }
  return i;
}

// simulate sfs for ancestral population
double simfreq() {
  return exp(uniform()*myConst-myConst);
}

// simulate sfs for 2 subpopulations using Balding-Nichols distribution
double simfreqBN(double F, double p_anc) {
  // FST values for the 2 subpops
  // p_anc: ancestral population allele frequency, drawn from a truncated exponential distribution (some other authors use a uniform distribution instead)
  return rbeta( ((1-F)/F)*p_anc, ((1-F)/F)*(1-p_anc)) ;
}

// to append names
char *append(const char* a,const char *b){
  char *c =(char *) malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}

// help printout
void info() {
  fprintf(stderr,"\t -> Required arg:\n\t\t-outfiles PREFIX\t PREFIX.seq PREFIX.glf PREFIX.frq PREFIX.arg\n");
  fprintf(stderr,"\t -> Optional arg:\n\t\t-npop\tNumber of populations. This MUST be set before -nind [1]\n");
  fprintf(stderr,"\t\t-nind\tNumber of diploid individuals for each population [10]\n");
  fprintf(stderr,"\t\t-nsites\tNumber of sites [500000]\n");
  fprintf(stderr,"\t\t-errate\tThe sequencing error rate [0.0075]\n");
  fprintf(stderr,"\t\t-depth\tMean sequencing depth [5]\n");
  fprintf(stderr,"\t\t-pvar\tProbability that a site is variable in the population [0.015]\n");
  fprintf(stderr,"\t\t-mfreq\tMinimum population frequency [0.0005]\n");
  fprintf(stderr,"\t\t-F\tFST value of 1st and 2nd split [0.4 0.1] OR inbreeding value/file in case of 1 pop [0]\n");
  fprintf(stderr,"\t\t-model\t0=fixed errate 1=variable errate [1]\n");
  fprintf(stderr,"\t\t-simpleRand\tboolean [1]\n");
  fprintf(stderr,"\t\t-seed\trandom number [0]\n");
  fprintf(stderr,"\t\t-base_freq\tBackground allele frequencies for A,C,G,T [0.25 0.25 0.25 0.25]\n");
  fprintf(stderr,"\t\t-multi_depth\tSimulate uneven covered individuals. -multi_depth 6 10: first 10 individuals have 6X while the rest is as -depth.[0 0]\n");
}

///		///
/// MAIN	///
///		///

int main(int argc, char *argv[]) { // read input parameters

  if(argc==1) { // if no argument (call of the program is considered as 1)
    info(); // return info
    return 0; // terminate
  }

  /// define and initialize the variables (with default values)
  
  int i=0, j=0, k=0, b1=0, b2=0, var=0, nsites = 500000, nind = 10, npop=1, model=1, nind1=0, nind2=0, nind3=0, increment=0, seed=0, multi_depth_coverage=0, multi_depth_nind=0;
  static int genotype[2], genotype1[2], genotype2[2], genotype3[2]; // array /matrix of genotypes for all pops
  double pfreq=0.0, pfreq1=0.0, pfreq2=0.0, pfreq3=0.0, pfreqB=0.0, pvar= 0.015, meandepth = 5, errate = 0.0075, F=0.0, F1=0.0, F2=0.0, minfreq=0.0001;
  double basefreq[4] = {0.25, 0.25, 0.25, 0.25}; // background frequencies
  double* indF; //per individual F
 
  int debug = 0; // change to 1 for debugging
  
  // for debug
  static int basecheck[4], basecheck1[4], basecheck2[4], basecheck3[4];

  //filehandles and their associated names:
  // resultfile: reads
  // glffile: genotype likelihoods per site
  // freqfile: whole SFS
  // argfile: input parameters
  //  gzFile glffile, resultfile, glffile1, glffile2, resultfile1, resultfile2, ireadsfile, ireadsfile1, ireadsfile2;

  gzFile glffile, resultfile, glffile1, glffile2, glffile3, resultfile1, resultfile2, resultfile3;
  FILE  *freqfile, *argfile, *genofile, *freqfile1, *freqfile2, *freqfile3, *genofile1, *genofile2, *genofile3, *freqfile12, *ireadsfile, *ireadsfile1, *ireadsfile2, *ireadsfile3, *parfile;
  char *fGlf=NULL,*fFreq=NULL,*fSeq=NULL,*fArg=NULL,*fGeno=NULL,*fireads=NULL;
  char *fGlf1=NULL,*fFreq1=NULL,*fSeq1=NULL,*fGeno1=NULL,*fireads1=NULL;
  char *fGlf2=NULL,*fFreq2=NULL,*fSeq2=NULL,*fGeno2=NULL,*fireads2=NULL;
  char *fGlf3=NULL,*fFreq3=NULL,*fSeq3=NULL,*fGeno3=NULL,*fireads3=NULL;
  char *fFreq12=NULL;
  char *outfiles=NULL;
  char *fparfile=NULL;
  
  // read input parameters and assign to values variables (if defined)
  int argPos=1; // we checked before if there were no inputs
  increment=0; // increment over the input parameters (one input may have more values than 1)
  
  while (argPos<argc) { // for all inputs

    // increment of 2 is for only 1 value input, otherwise is larger 
    increment=0; // increment over the input parameters (one input may have more values than 1)

    if(strcmp(argv[argPos],"-npop")==0) npop = atoi(argv[argPos+1]); 

    else if(strcmp(argv[argPos],"-nind")==0) {
      increment = increment + npop -1; // argPos will be then 1+1=2 if 2 pops
      if (npop==1) {
	nind = atoi(argv[argPos+1]);
	nind1 = nind2 = 0;
      } else if (npop==2) {
	nind1 = atoi(argv[argPos+1]);
	nind2 = atoi(argv[argPos+2]);
	nind = nind1+nind2;
      } else if (npop==3) {
	nind1 = atoi(argv[argPos+1]);
	nind2 = atoi(argv[argPos+2]);
	nind3 = atoi(argv[argPos+3]);
	nind = nind1+nind2+nind3;
      } else {
	printf("\n-npop must be 1 or 2 or 3. Terminate.\n");
	return 0;
      } // end if npop
    }

    else if(strcmp(argv[argPos],"-F")==0) {
      increment = increment + npop -1; // argPos will be then 1+1=2 if 2 pops

      if (npop==1) {
	char *str;
	F = strtod(argv[argPos+1], &str); // inbreeding
	if(strcmp(str,"")!=0)
	  F=argPos+1;
	else if(F<0 || F>1) {
	  printf("error in F (%f); must be [0,1]\n", F); 
	  exit(-1);
	}
      } else if (npop==2) {
	F1 = atof(argv[argPos+1]); // FST
	F2 = atof(argv[argPos+2]); // FST (redundant)
	F = 0.0; // no inbreeding"
      } else if (npop==3) {
	F1 = atof(argv[argPos+1]); // this is FST first split
	F2 = atof(argv[argPos+2]); // FST second split
	F = 0.0; // no inbreeding"
        increment = increment -1;
      } // end if npop

    }

    // argPos+1 because you count the -flag too (argPos is 1!)
    // argPos+1+npop-1 == argPos+npop; es. 2 pops, increment is 3 for the next inpute read 

    else if(strcmp(argv[argPos],"-errate")==0)   errate  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-depth")==0) meandepth  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-pvar")==0)  pvar = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-mfreq")==0)  minfreq = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-outfiles")==0) outfiles = (argv[argPos+1]);
    else if(strcmp(argv[argPos],"-nsites")==0)  nsites = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-seed")==0)  seed = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-multi_depth")==0) {
      multi_depth_coverage  = atoi(argv[argPos+1]);
      multi_depth_nind  = atoi(argv[argPos+2]);
      increment = increment + 1;
      }
    else if(strcmp(argv[argPos],"-model")==0)  model = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-simpleRand")==0) simpleRand = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-base_freq")==0) {
     increment = increment + 3;
     basefreq[0] = atof(argv[argPos+1]); basefreq[1] = atof(argv[argPos+2]); basefreq[2] = atof(argv[argPos+3]); basefreq[3] = atof(argv[argPos+4]); 
    }

    else { // input is not a valid one 
     
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0; // terminate
      
    }
    
    argPos = argPos + 2 + increment; // increment next value (+2 because un count the input flags too)
  
  } // end while all inputs
  
  // check if there is the output file
  if(outfiles == NULL) {
    fprintf(stderr,"\nMust supply -outfiles. Terminate\n");
    return 0;
  }
  
  // Initialize indF
  indF = (double*) malloc(nind*sizeof(double));
  // Read indF file
  if(F > 1) {
    int cnt = 0;
    char buf[100];
    FILE* indF_f = getFile(argv[int(F)], "r");
    while( fgets(buf,100,indF_f) )
      indF[cnt++] = atof(buf);
    fclose(indF_f);
    F = indF[0];
  } else {
    for (int i = 0; i < nind; i++)
      indF[i] = F;
  }
  
  // check if multi_depth_nind is higher than nind and if you have multiple populations
  if(multi_depth_nind>nind) {
    fprintf(stderr,"\nmulti_depth individuals must be lower than total number of individuals -npop. Terminate\n");
    return 0;
  }
  if(multi_depth_nind>0 && npop>1) {
    fprintf(stderr,"\nmulti_depth supported only if simulating 1 population. Terminate\n");
    return 0;
  }



  // to adjust the exponential function according to the lowest allele frequency detectable
  myConst = -log(minfreq);
  
  // print input arguments
  fprintf(stderr,"\t->Using args: -npop %d -nind %d -nind1 %d -nind2 %d -nind3 %d -errate %f -depth %f -pvar %f -mfreq %f -nsites %d -F %f -F1 %f -F2 %f -model %d -simpleRand %d -seed %d -base_freq %f %f %f %f \n", npop, nind, nind1, nind2, nind3, errate, meandepth, pvar, minfreq, nsites, F, F1, F2, model, simpleRand, seed, basefreq[0], basefreq[1], basefreq[2], basefreq[3]);
  
  /// output files
  fArg = append(outfiles,".args");
  fGlf = append(outfiles,".glf.gz");
  fireads = append(outfiles,".reads.txt");
  fFreq = append(outfiles,".frq");
  fSeq = append(outfiles,".seq.gz");
  fGeno = append(outfiles,".geno");
  fparfile = append(outfiles, ".par");
  // print a message
  fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile: %s\ttruefreq: %s args:%s geno:%s reads:%s par:%s\n",fSeq,fGlf,fFreq,fArg,fGeno,fireads,fparfile);

  // get files
  argfile=getFile(fArg,"w");
  resultfile=getGz(fSeq,"w");
  glffile = getGz(fGlf,"w");
  ireadsfile = getFile(fireads,"w");
  freqfile= getFile(fFreq,"w"); 
  genofile =getFile(fGeno,"w");
  parfile = getFile(fparfile, "w");

  // eventually if 2 pops, prepare other files and print messages
  if (npop==2) {
    
    fGlf1 = append(outfiles,"1.glf.gz");
    fireads1 = append(outfiles,"1.reads.txt");
    fFreq1 = append(outfiles,"1.frq");
    fSeq1 = append(outfiles,"1.seq.gz");
    fGeno1 = append(outfiles,"1.geno");
    fGlf2 = append(outfiles,"2.glf.gz");
    fireads2 = append(outfiles,"2.reads.txt");
    fFreq2 = append(outfiles,"2.frq");
    fSeq2 = append(outfiles,"2.seq.gz");
    fGeno2 = append(outfiles,"2.geno");
    fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile1: %s\ttruefreq1: %s geno1:%s reads1: %s\n", fSeq1, fGlf1, fFreq1, fGeno1,fireads1);
    fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile2: %s\ttruefreq2: %s geno2:%s reads2: %s\n", fSeq2, fGlf2, fFreq2, fGeno2,fireads2);

    resultfile1=getGz(fSeq1, "w");
    glffile1 = getGz(fGlf1, "w");
    ireadsfile1 = getFile(fireads1, "w");
    freqfile1= getFile(fFreq1, "w"); 
    genofile1 =getFile(fGeno1, "w");
    resultfile2=getGz(fSeq2, "w");
    glffile2 = getGz(fGlf2, "w");
    ireadsfile2 = getFile(fireads2, "w");
    freqfile2= getFile(fFreq2, "w"); 
    genofile2 =getFile(fGeno2, "w");

    // joint-SFS
    fFreq12 = append(outfiles,"12.frq");    
    freqfile12 = getFile(fFreq12, "w"); 
  
  }

 // eventually if 3 pops, prepare other files and print messages
  if (npop==3) {
    
    fGlf1 = append(outfiles,"1.glf.gz");
    fireads1 = append(outfiles,"1.reads.txt");
    fFreq1 = append(outfiles,"1.frq");
    fSeq1 = append(outfiles,"1.seq.gz");
    fGeno1 = append(outfiles,"1.geno");
    fGlf2 = append(outfiles,"2.glf.gz");
    fireads2 = append(outfiles,"2.reads.txt");
    fFreq2 = append(outfiles,"2.frq");
    fSeq2 = append(outfiles,"2.seq.gz");
    fGeno2 = append(outfiles,"2.geno");
    fGlf3 = append(outfiles,"3.glf.gz");
    fireads3 = append(outfiles,"3.reads.txt");
    fFreq3 = append(outfiles,"3.frq");
    fSeq3 = append(outfiles,"3.seq.gz");
    fGeno3 = append(outfiles,"3.geno");

    fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile1: %s\ttruefreq1: %s geno1:%s reads1: %s\n", fSeq1, fGlf1, fFreq1, fGeno1,fireads1);
    fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile2: %s\ttruefreq2: %s geno2:%s reads2: %s\n", fSeq2, fGlf2, fFreq2, fGeno2,fireads2);
    fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile3: %s\ttruefreq3: %s geno3:%s reads3: %s\n", fSeq3, fGlf3, fFreq3, fGeno3,fireads3);

    resultfile1=getGz(fSeq1, "w");
    glffile1 = getGz(fGlf1, "w");
    ireadsfile1 = getFile(fireads1, "w");
    freqfile1= getFile(fFreq1, "w"); 
    genofile1 =getFile(fGeno1, "w");

    resultfile2=getGz(fSeq2, "w");
    glffile2 = getGz(fGlf2, "w");
    ireadsfile2 = getFile(fireads2, "w");
    freqfile2= getFile(fFreq2, "w"); 
    genofile2 =getFile(fGeno2, "w");

    resultfile3=getGz(fSeq3, "w");
    glffile3 = getGz(fGlf3, "w");
    ireadsfile3 = getFile(fireads3, "w");
    freqfile3= getFile(fFreq3, "w"); 
    genofile3 =getFile(fGeno3, "w");


    // joint-SFS
    fFreq12 = append(outfiles,"12.frq"); // maybe useless here?   
    freqfile12 = getFile(fFreq12, "w"); 
  
  }
  

  // write args file
  fprintf(argfile,"\t->Using args: -npop %d -nind %d -nind1 %d -nind2 %d -errate %f -depth %f -pvar %f -mfreq %f -nsites %d -F %f -F1 %f -F2 %f -model %d -simpleRand %d -seed %d -base_freq %f %f %f %f -multi_depth %d %d\n", npop, nind, nind1, nind2, errate, meandepth, pvar, minfreq, nsites, F, F1, F2, model, simpleRand, seed, basefreq[0], basefreq[1], basefreq[2], basefreq[3], multi_depth_coverage, multi_depth_nind); 

  /// COMPUTE
  
  // initialize SFS (whole locus), as the unfolded spectrum (ancestral and derived)
  // as arrays
  int freqspec[nind*2+1], freqspec1[nind1*2+1], freqspec2[nind2*2+1], freqspec3[nind3*2+1];
  int freqspec12[nind1*2+1][nind2*2+1]; // joint-SFS, matrix

  for (i=0; i<nind*2+1; i++) freqspec[i]=0; // initialize sfs (from 0 to 2*nind)
  for (i=0; i<nind1*2+1; i++) freqspec1[i]=0;
  for (i=0; i<nind2*2+1; i++) freqspec2[i]=0;
  for (i=0; i<nind3*2+1; i++) freqspec3[i]=0;
  for (i=0; i<nind1*2+1; i++) {
    for (j=0; j<nind2*2+1; j++) {
      freqspec12[i][j]=0;
    }
  } // end initialize joint-SFS fo first 2 pops

  /// FOR EACH SITE
  
  for (i=0; i<nsites; i++) {
    
    /*debug code*/
    for (k=0; k<4; k++) basecheck[k]=basecheck1[k]=basecheck2[k]=basecheck3[k]=0; // initialize base checks      
    
    // basechecks are arrays of 4 elements counting the nr of alleles 0,1,2,3
    
    /// test if site is variable or not
    
    if (uniform(seed) >= pvar) { // if it is NOT variable... 

      var=0; // not variable (boolean)
      // assign alleles to genotype
      genotype[0]=genotype[1]=0; // all ancestral A

      /*debug code*/
      basecheck[genotype[0]] = 2*nind; // all individuals are monorphic for ancestral A, basechecks are arrays of 4 elements counting the nr of alleles 0,1,2,3
      if (npop==2) {
	genotype1[0]=genotype1[1]=genotype2[0]=genotype2[1]=0;
	basecheck1[genotype1[0]] = 2*nind1; // all individuals are monomorphic
	basecheck2[genotype2[0]] = 2*nind2; // all individuals are monomorphic
      }
      if (npop==3) {
	genotype1[0]=genotype1[1]=genotype2[0]=genotype2[1]=genotype3[0]=genotype3[1]=0;
	basecheck1[genotype1[0]] = 2*nind1; // all individuals are monomorphic
	basecheck2[genotype2[0]] = 2*nind2; // all individuals are monomorphic
	basecheck3[genotype3[0]] = 2*nind3; // all individuals are monomorphic
      }      

      fprintf(parfile, "%f\t%f\t%f\t%f\n", 0.0, F, F1, F2);

    } else { // if it IS variable
      
      var = 1; // TRUE 
      
      // choose the 2 alleles
      b2 = 0; //changed such that reference is always A, for convenience
      while ((b1=pick_base_from_prior(basefreq))==b2); // take second allele, different from the first A
      
      // we will assign genotypes on the basis of the 2 alleles later, not here
      // here it is just a check if the site is variable and which alleles to take.
	
      // simulate population allele frequency (or ancestral if 2/3 subpops)
      pfreq=simfreq(); // if site is not variable, you don't need to compute pfreq
      
      fprintf(parfile, "%f\t%f\t%f\t%f\n", pfreq, F, F1, F2);

    } // end test if it is variable or not (pvar)
    
    /// WHOLE POPULATIONS (or the only one if -npop=1)
 
    // now for each individual, assigned genotypes at each individual as couples of alleles/bases
   if (npop==1) {
    for (j=0; j<nind; j++) {

      if (var==1) {

	if (uniform(seed)>=indF[j]) { //no inbreeding case
	  for (k=0; k<2; k++) {
	    if (uniform(seed)<=pfreq) 
	      genotype[k] = b1;
	    else genotype[k] = b2; 
	  }
	} else { //inbreeding case
	  if (uniform(seed)<=pfreq) {
	    genotype[0] = b1;
	    genotype[1] = b1;
	  } else {
	    genotype[0] = b2; 
	    genotype[1] = b2;
	  }
	} // end if uniform<=> F
	
	/*debug code*/ 
	// increment the count for each allele type on the basis of just assigned genotype
	basecheck[genotype[0]]++; basecheck[genotype[1]]++;
	
      } // end assignment of genotypes of site is variable
      
      // write genotypes in the output file (append it)
      int has_reads =0;
      // compute and print likelihoods
      if (debug) fprintf(stderr,"\nStart writing reads for whole sample; indiv %d", j);      
      if (model==0) {
        if (multi_depth_coverage>0 && j<multi_depth_nind) {
          has_reads=print_ind_site(errate, multi_depth_coverage, genotype, resultfile, glffile, ireadsfile);
        } else {        
	  has_reads=print_ind_site(errate, meandepth, genotype, resultfile, glffile, ireadsfile);
        }
      } else {
        if (multi_depth_coverage>0 && j<multi_depth_nind) {
          has_reads=print_ind_site(2*errate*uniform(seed), multi_depth_coverage, genotype, resultfile, glffile, ireadsfile);
        } else {
	has_reads=print_ind_site(2*errate*uniform(seed), meandepth, genotype, resultfile, glffile, ireadsfile);
        }
      }
      fprintf(genofile,"%d %d\t",genotype[0],genotype[1]);      
      
      if (debug) fprintf(stderr,"\nEnd writing reads for whole sample");      
      
      // now write binary files
      if (j<nind-1) {
	if(dumpBinary) {
	  char sep[1]={'\t'};
	  gzwrite(resultfile,sep,1);
	} else {
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);
	  //fprintf(glffile,",");
	  // fprintf(resultfile,", ");
	}
      } // if nind

    } // end for j in nind, all individuals
   } // if npop==1
    
    /// MULTI POPULATIONS
    
    if (npop==2) {

      if (debug) fprintf(stderr,"\n\n START 2 POPS");      
      
      // then assign subpop freq based on time-split? theory of pop genet
  
      // use now balding-nichols and write one set of files for each population
      // simulate distinct pops allele freq from Balding-Nickols distribution, given an ancestral population allele frequency and an FST value
      // 2 independent draws
      // pfreq already simulated!
      pfreq1=simfreqBN(F1, pfreq);
      pfreq2=simfreqBN(F2, pfreq);
     
      /// 1st pop (same as the whole)
      for (j=0; j<nind1; j++) {
	if (var==1) {
	  if (uniform(seed)>=F) { // this will alwasy be TRUE now 
	    for (k=0; k<2; k++) {
	      if (uniform(seed)<=pfreq1) 
		genotype1[k] = b1;
	    else genotype1[k] = b2; 
	  }
	} else {
	  if (uniform(seed)<=pfreq1	) {
	    genotype1[0] = b1;
	    genotype1[1] = b1;
	  } else {
	    genotype1[0] = b2; 
	    genotype1[1] = b2;
	  }
	}	
	basecheck1[genotype1[0]]++; basecheck1[genotype1[1]]++;
	basecheck[genotype1[0]]++; basecheck[genotype1[1]]++;
	}
      fprintf(genofile1,"%d %d\t",genotype1[0],genotype1[1]);     
      // write also into whole genotype file
      fprintf(genofile,"%d %d\t",genotype1[0],genotype1[1]);
      int has_reads1 =0, has_reads=0;
      if (debug) fprintf(stderr,"\nStart writing reads for 1 sample");      
      if (model==0) {
	has_reads1=has_reads=print_ind_site2(errate, meandepth, genotype1, resultfile, glffile, ireadsfile, resultfile1, glffile1, ireadsfile1, 0);
      } else {
        has_reads1=has_reads=print_ind_site2(2*errate*uniform(seed), meandepth, genotype1, resultfile, glffile, ireadsfile, resultfile1, glffile1, ireadsfile1, 0);
      }
      if (debug) fprintf(stderr,"\nEnd writing reads for 1 sample");      
      if (j<nind1-1) {
	if(dumpBinary) {
	  char sep[1]={'\t'};
	  gzwrite(resultfile1, sep, 1);
	  gzwrite(resultfile, sep, 1); // write also into whole genotype file
	} else {
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);
	}
      }
      }

      /// 2nd pop (same as above)
      for (j=0; j<nind2; j++) {
	if (var==1) {
	  if (uniform(seed)>F) { // this will always be TRUE now
	    for (k=0; k<2; k++) {
	      if (uniform(seed)<=pfreq2) 
		genotype2[k] = b1;
	    else genotype2[k] = b2; 
	  }
	} else {
	  if (uniform(seed)<=pfreq2	) {
	    genotype2[0] = b1;
	    genotype2[1] = b1;
	  } else {
	    genotype2[0] = b2; 
	    genotype2[1] = b2;
	  }
	}	
	basecheck2[genotype2[0]]++; basecheck2[genotype2[1]]++;
	basecheck[genotype2[0]]++; basecheck[genotype2[1]]++;
      }

      fprintf(genofile2,"%d %d\t",genotype2[0],genotype2[1]);      
      fprintf(genofile,"%d %d\t",genotype2[0],genotype2[1]);      
      int has_reads2 =0, has_reads=0;
      if (debug) fprintf(stderr,"\nStart writing reads for 2 sample");      
      if (model==0) {
	has_reads2=has_reads=print_ind_site2(errate, meandepth, genotype2, resultfile, glffile, ireadsfile, resultfile2, glffile2, ireadsfile2, (j+1));
      } else {
	has_reads2=has_reads=print_ind_site2(2*errate*uniform(seed), meandepth, genotype2, resultfile, glffile, ireadsfile, resultfile2, glffile2, ireadsfile2, (j+1));
      }
      if (debug) fprintf(stderr,"\nEnd writing reads for 2 sample");    

      if (j<nind2-1) {
	if(dumpBinary) {
	  char sep[1]={'\t'};
	  gzwrite(resultfile2, sep, 1);
	  gzwrite(resultfile, sep, 1);
	} else {
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);

	}
      }
    }
    
    } // end if multi populations npop==2


    // 3 POPULATIONS

    if (npop==3) {

      if (debug) fprintf(stderr,"\n\n START 3 POPS");      
      

      // use now balding-nichols and write one set of files for each population
      // simulate distinct pops allele freq from Balding-Nickols distribution, given an ancestral population allele frequency and an FST value

      // 2 independent draws
      // pfreq already simulated!
      pfreq1=simfreqBN(F1, pfreq);
      pfreqB=simfreqBN(F1, pfreq);
      // in case of 3 pops then second draw is the second ancestra/
      pfreq2=simfreqBN(F2, pfreqB);
      pfreq3=simfreqBN(F2, pfreqB);


      /// 1st pop (same as the whole)
      for (j=0; j<nind1; j++) {
	if (var==1) {
	  if (uniform(seed)>=F) { // this will alwasy be TRUE now 
	    for (k=0; k<2; k++) {
	      if (uniform(seed)<=pfreq1) 
		genotype1[k] = b1;
	    else genotype1[k] = b2; 
	  }
	} else {
	  if (uniform(seed)<=pfreq1	) {
	    genotype1[0] = b1;
	    genotype1[1] = b1;
	  } else {
	    genotype1[0] = b2; 
	    genotype1[1] = b2;
	  }
	}	
	basecheck1[genotype1[0]]++; basecheck1[genotype1[1]]++;
	basecheck[genotype1[0]]++; basecheck[genotype1[1]]++;
	}
      fprintf(genofile1,"%d %d\t",genotype1[0],genotype1[1]);     
      // write also into whole genotype file
      fprintf(genofile,"%d %d\t",genotype1[0],genotype1[1]);
      int has_reads1 =0, has_reads=0;
      if (debug) fprintf(stderr,"\nStart writing reads for 1 sample");      
      if (model==0) {
	has_reads1=has_reads=print_ind_site2(errate, meandepth, genotype1, resultfile, glffile, ireadsfile, resultfile1, glffile1, ireadsfile1, 0);
      } else {
        has_reads1=has_reads=print_ind_site2(2*errate*uniform(seed), meandepth, genotype1, resultfile, glffile, ireadsfile, resultfile1, glffile1, ireadsfile1, 0);
      }
      if (debug) fprintf(stderr,"\nEnd writing reads for 1 sample");      
      if (j<nind1-1) {
	if(dumpBinary) {
	  char sep[1]={'\t'};
	  gzwrite(resultfile1, sep, 1);
	  gzwrite(resultfile, sep, 1); // write also into whole genotype file
	} else {
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);
	}
      }
      }

      /// 2nd pop (same as above)
      for (j=0; j<nind2; j++) {
	if (var==1) {
	  if (uniform(seed)>F) { // this will always be TRUE now
	    for (k=0; k<2; k++) {
	      if (uniform(seed)<=pfreq2) 
		genotype2[k] = b1;
	    else genotype2[k] = b2; 
	  }
	} else {
	  if (uniform(seed)<=pfreq2	) {
	    genotype2[0] = b1;
	    genotype2[1] = b1;
	  } else {
	    genotype2[0] = b2; 
	    genotype2[1] = b2;
	  }
	}	
	basecheck2[genotype2[0]]++; basecheck2[genotype2[1]]++;
	basecheck[genotype2[0]]++; basecheck[genotype2[1]]++;
      }

      fprintf(genofile2,"%d %d\t",genotype2[0],genotype2[1]);      
      fprintf(genofile,"%d %d\t",genotype2[0],genotype2[1]);      
      int has_reads2 =0, has_reads=0;
      if (debug) fprintf(stderr,"\nStart writing reads for 2 sample");      
      if (model==0) {
	has_reads2=has_reads=print_ind_site2(errate, meandepth, genotype2, resultfile, glffile, ireadsfile, resultfile2, glffile2, ireadsfile2, (j+1));
      } else {
	has_reads2=has_reads=print_ind_site2(2*errate*uniform(seed), meandepth, genotype2, resultfile, glffile, ireadsfile, resultfile2, glffile2, ireadsfile2, (j+1));
      }
      if (debug) fprintf(stderr,"\nEnd writing reads for 2 sample");    

      if (j<nind2-1) {
	if(dumpBinary) {
	  char sep[1]={'\t'};
	  gzwrite(resultfile2, sep, 1);
	  gzwrite(resultfile, sep, 1);
	} else {
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);

	}
      }
    }
    
    /// 3rd pop (same as above), exactly as pop 2 but use genotype3, files3...

     for (j=0; j<nind3; j++) {
	if (var==1) {
	  if (uniform(seed)>F) { // this will always be TRUE now
	    for (k=0; k<2; k++) {
	      if (uniform(seed)<=pfreq3) 
		genotype3[k] = b1;
	    else genotype3[k] = b2; 
	  }
	} else {
	  if (uniform(seed)<=pfreq3) {
	    genotype3[0] = b1;
	    genotype3[1] = b1;
	  } else {
	    genotype3[0] = b2; 
	    genotype3[1] = b2;
	  }
	}	
	basecheck3[genotype3[0]]++; basecheck3[genotype3[1]]++;
	basecheck[genotype3[0]]++; basecheck[genotype3[1]]++;
      }

      fprintf(genofile3,"%d %d\t",genotype3[0],genotype3[1]);      
      fprintf(genofile,"%d %d\t",genotype3[0],genotype3[1]);
      
      int has_reads3 =0, has_reads=0;
      if (debug) fprintf(stderr,"\nStart writing reads for 2 sample");      
      if (model==0) {
	has_reads3=has_reads=print_ind_site2(errate, meandepth, genotype3, resultfile, glffile, ireadsfile, resultfile3, glffile3, ireadsfile3, (j+1));
      } else {
	has_reads3=has_reads=print_ind_site2(2*errate*uniform(seed), meandepth, genotype3, resultfile, glffile, ireadsfile, resultfile3, glffile3, ireadsfile3, (j+1));
      }
      if (debug) fprintf(stderr,"\nEnd writing reads for 2 sample");    

      if (j<nind3-1) {
	if(dumpBinary) {
	  char sep[1]={'\t'};
	  gzwrite(resultfile3, sep, 1);
	  gzwrite(resultfile, sep, 1);
	} else {
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);

	}
      }
    }

    } // end if multi populations npop==3

    // EOL
    fprintf(genofile,"\n");
    if (npop>1) { 
      fprintf(genofile1,"\n"); 
      fprintf(genofile2,"\n"); 
    }
    if (npop==3) { 
      fprintf(genofile3,"\n"); 
    }

    /// DEBUG CODE AND COMPUTE SFS

    if (1) { 

      int k=0, j=0, j1=0, j2=0, j3=0;
      for (k=1; k<4; k++) { // count only derived alleles
	if (basecheck[k]>j) j=basecheck[k];
      } // j is the count of the derived allele   
      if ((j+basecheck[0])!=nind*2) { // if sum of allele-counts is not equal to nind*2 (nr chromosomes)
	printf("error in freqspec calculation (%i %i %i %i)\n", basecheck[0], basecheck[1], basecheck[2], basecheck[3]); 
	exit(-1);
      }
      // then now j is the count of non-ancestral alleles (not 0)
      freqspec[j]++; // this if unfolded

      if (npop>1) {

	k=0, j1=0;
	for (k=1; k<4; k++) { // count only derived alleles
	  if (basecheck1[k]>j1) j1=basecheck1[k];
	} // j is the count of the derived allele     
	if ((j1+basecheck1[0])!=nind1*2) { // if sum of allele-counts is not equal to nind*2 (nr chromosomes)
	  printf("error in freqspec1 calculation (%i %i %i %i)\n", basecheck1[0], basecheck1[1], basecheck1[2], basecheck1[3]); 
	  exit(-1);
	}
	// then now j is the count of non-ancestral alleles (not 0)
	freqspec1[j1]++; // this if unfolded

	k=0, j2=0;
	for (k=1; k<4; k++) { // count only derived alleles
	  if (basecheck2[k]>j2) j2=basecheck2[k];
	} // j is the count of the derived allele     
	if ((j2+basecheck2[0])!=nind2*2) { // if sum of allele-counts is not equal to nind*2 (nr chromosomes)
	  printf("error in freqspec2 calculation (%i %i %i %i)\n", basecheck2[0], basecheck2[1], basecheck2[2], basecheck2[3]); 
	  exit(-1);
	}
	// then now j is the count of non-ancestral alleles (not 0)
	freqspec2[j2]++; // this if unfolded
	
	// joint-SFS
	freqspec12[j1][j2]++;
      
      if (npop==3) {

	k=0, j3=0;
	for (k=1; k<4; k++) { // count only derived alleles
	  if (basecheck3[k]>j3) j3=basecheck3[k];
	} // j is the count of the derived allele

    
	if ((j3+basecheck3[0])!=nind3*2) { // if sum of allele-counts is not equal to nind*2 (nr chromosomes)
	  printf("error in freqspec3 calculation:base (%i %i %i %i) j3 %d tot %d:\n", basecheck3[0], basecheck3[1], basecheck3[2], basecheck3[3],j3,nind3*2);

	  exit(-1);
	}
	// then now j is the count of non-ancestral alleles (not 0)
	freqspec3[j3]++; // this if unfolded
	
      } // end if npop 3

     } // end if npop>1

    } // end if (1) 

   
    /// write reads and GLF
    if (dumpBinary) {
      char sep[1] = {'\n'};
      gzwrite(resultfile, sep, 1);
    } else {
      fprintf(stderr, "non binary output disabled\n");
      exit(0);
    }
    if (npop>1) {
      // 1st
      if (dumpBinary) {
        char sep[1] = {'\n'};
        gzwrite(resultfile1, sep, 1);
      } else {
        fprintf(stderr, "non binary output disabled\n");
       exit(0);
      }
      // 2nd
      if (dumpBinary) {
        char sep[1] = {'\n'};
        gzwrite(resultfile2, sep, 1);
      } else {
        fprintf(stderr, "non binary output disabled\n");
        exit(0);
      }
    } // end if pop>1

    if (npop==3) {
      if (dumpBinary) {
        char sep[1] = {'\n'};
        gzwrite(resultfile3, sep, 1);
      } else {
        fprintf(stderr, "non binary output disabled\n");
        exit(0);
      }

    } // end dump binary if 3 pops

    // computed and written for a site

  } // end for in i each sites

  /// COMPUTE WHOLE SFS
  k=0;
  for (i=0; i<nind*2+1; i++)
    k=k+freqspec[i]; // how many alleles in the whole sample
  for (i=0; i<nind*2+1; i++) { // compute relative frequency for each allele freq bin
    fprintf(freqfile,"%f\t",(double)freqspec[i]/(double)k);
  }
  fprintf(freqfile,"\n"); // EOL
       
  // npop==2
  if (npop>1) {
    k=0;
    for (i=0; i<nind1*2+1; i++)
      k=k+freqspec1[i]; // how many alleles in the whole sample
    for (i=0; i<nind1*2+1; i++) { // compute relative frequency for each allele freq bin
      fprintf(freqfile1,"%f\t",(double)freqspec1[i]/(double)k);
    }
    fprintf(freqfile1,"\n"); // EOL

    k=0;
    for (i=0; i<nind2*2+1; i++)
    k=k+freqspec2[i]; // how many alleles in the whole sample
    for (i=0; i<nind2*2+1; i++) { // compute relative frequency for each allele freq bin
      fprintf(freqfile2,"%f\t",(double)freqspec2[i]/(double)k);
    }
    fprintf(freqfile2,"\n"); // EOL 
    
    k=0; // joint
    for (i=0; i<nind1*2+1; i++) {
      for (j=0; j<nind2*2+1; j++) {
	k=k+freqspec12[i][j];
      }
    }

    for (i=0; i<nind1*2+1; i++) {
      for (j=0; j<nind2*2+1; j++) {
	fprintf(freqfile12,"%f\t",(double)freqspec12[i][j]/(double)k);
      }
      fprintf(freqfile12,"\n"); // EOL
    }
    
  } // end if 2 pops

  // add 3
  if (npop==3) {
    k=0;
    for (i=0; i<nind3*2+1; i++)
    k=k+freqspec3[i]; // how many alleles in the whole sample
    for (i=0; i<nind3*2+1; i++) { // compute relative frequency for each allele freq bin
      fprintf(freqfile3,"%f\t",(double)freqspec3[i]/(double)k);
    }
    fprintf(freqfile3,"\n"); // EOL 
  }

  /// FREE MEMORY, FLUSH AND CLOSE
  
  free(fGlf);
  free(fireads);
  free(fFreq);
  free(fSeq);
  free(fArg);
  free(fparfile);

  gzclose(resultfile); //fclose flushed automaticly
  gzclose(glffile);
  fclose(ireadsfile);
  fclose(argfile);
  fclose(freqfile);
  fclose(genofile);
  fclose(parfile);

  if (npop>1) {

    free(fGlf1);
    free(fireads1);
    free(fFreq1);
    free(fSeq1);

    free(fGlf2);
    free(fireads2);
    free(fFreq2);
    free(fSeq2);
   
    free(fFreq12);
    
    gzclose(resultfile1); //fclose flushed automaticly
    gzclose(glffile1);
    fclose(ireadsfile1);
    fclose(freqfile1);
    fclose(genofile1);
    gzclose(resultfile2); //fclose flushed automaticly
    gzclose(glffile2);
    fclose(ireadsfile2);
    fclose(freqfile2);
    fclose(genofile2);

    fclose(freqfile12);
    
  }

  if (npop==3) {

    free(fGlf3);
    free(fireads3);
    free(fFreq3);
    free(fSeq3);

    gzclose(resultfile3); //fclose flushed automaticly
    gzclose(glffile3);
    fclose(ireadsfile3);
    fclose(freqfile3);
    fclose(genofile3);

  }

  return 0; // return value

} // end main












