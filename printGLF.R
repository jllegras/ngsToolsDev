

# Usage: Rscript -i infile.stat -o outfile.eps -n name1-name2
# Usage: Rscript -i infile.stat -o outfile.eps -n name1 # if only 1 pop

library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file'),
                    make_option(c('-n','--names'), action='store', type='character', default=1-2, help='Name(s) of population(s)'),
                    make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))




readGLF <- function(file=NULL, ncat=10, nsites=10, nind=10) {
 ff <- gzfile(file,"rb")
 m<-matrix(readBin(ff,"double",ncat*nsites*nind),ncol=ncat,byrow=TRUE)
 close(ff)
 return(m)
}



