# High Performance Computers 
_aka_ *supercomputers* or *high performance computers* (HPC)

## QSUB and R

Need to edit R script to 'behave' like a normal #! (linux command line) script.  This is easy with [getopt](http://cran.r-project.org/web/packages/getopt/index.html) package. 
```{r,eval=F}
library(getopt)
## get options
opta <- getopt(matrix(c(
  'date', 'd', 1, 'character',
  'help', 'h', 0, 'logical'
), ncol=4, byrow=TRUE))
if ( !is.null(opta$help) )
{
  prg <- commandArgs()[1];
  cat(paste("Usage: ", prg,  " --date | -d <file> :: The date to process\n", sep=""));
  q(status=1);
}
## extract value
date=opta$date 
```

```{R,eval=F}
Rscript script.R --date 20131105 
```

## Driving cluster from R

Possible to drive the cluster from within R via QSUB.  First, define the jobs:
  ```{r,eval=FALSE}
script="/path/to/Rscript.r"

write.table(
  paste(script,"--date",dates),                     
  file="process.txt",
  row.names=F,col.names=F,quote=F)

### Set up submission script
queue="devel"
nodes=120
walltime=24
```

## Write the QSUB script

```{r,eval=F}
### write qsub script to disk from R
cat(paste("
#PBS -S /bin/bash
#PBS -l select=",nodes,":ncpus=8:mpiprocs=8
#PBS -l walltime=",walltime,":00:00
#PBS -q ",queue,"

CORES=",nodes*8,"

source $HDIR/etc/environ.sh
IDIR=/nobackupp1/awilso10/mod35/
WORKLIST=$IDIR/process.txt
EXE=Rscript
LOGSTDOUT=$IDIR/log/mod35_stdout
LOGSTDERR=$IDIR/log/mod35_stderr
          
### use mpiexec to parallelize across days
mpiexec -np $CORES pxargs -a $WORKLIST -p $EXE 1> $LOGSTDOUT 2> $LOGSTDERR
",sep=""),file=paste("mod35_qsub",sep=""))

## run it!
system(paste("qsub mod35_qsub",sep=""))
```

