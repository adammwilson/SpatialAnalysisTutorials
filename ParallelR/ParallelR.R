#' ---
#' title: "Introduction to Parallel Computing with R"
#' author: "Adam M. Wilson"
#' date: "November 5, 2014"
#' output:
#'   knitrBootstrap::bootstrap_document:
#'     highlight: Magula
#'     highlight.chooser: no
#'     theme: cerulean
#'     theme.chooser: no
#'   pdf_document:
#'     toc: true
#'     toc_depth: 2
#'   md_document:
#'     variant: markdown_github
#' ---
#' 
#' Introduction to Parallel Computing with R
#' ====
#' 
#' This script is available [here](https://rawgit.com/adammwilson/SpatialAnalysisTutorials/master/ParallelR/ParallelR.html).  You can view a version of this script with commented text [here](https://rawgit.com/adammwilson/SpatialAnalysisTutorials/master/ParallelR/ParallelR.R).
#'  
#'  
## ----message=F,warning=FALSE---------------------------------------------
library(foreach)
library(doParallel)
library(knitr)
library(raster)
library(rasterVis)
library(arm)
library(coda)
library(fields)
library(dplyr)
library(ggplot2)
library(ggmcmc)

#' 
## ----echo=FALSE----------------------------------------------------------
opts_chunk$set(cache=TRUE)
# purl("ParallelR/ParallelR.Rmd","ParallelR/ParallelR.R",documentation=2)
# rmarkdown::render("ParallelR/ParallelR.Rmd", "all")

presentation_theme <- theme_grey()+theme(text = element_text(size = 25, colour = "black"))
theme_set(presentation_theme)

#' 
#' If you don't have the packages above, install them in the package manager or by running `install.packages("doParallel")`. 
#' 
#' # Introduction
#' 
#' ## Serial Computing
#' Most (legacy) software is written for serial computation:
#' 
#'   * Problem broken into discrete set of instructions
#'   * Instructions executed sequentially on a single processor
#'   
#' ![https://computing.llnl.gov/tutorials/parallel_comp/](assets/serialProblem.gif)
#' <span style="color:grey; font-size:1em;">Figure from [here](https://computing.llnl.gov/tutorials/parallel_comp/) </span>
#' 
#' ## Parallel computation
#' 
#' Parallel computing is the simultaneous use of multiple compute resources:
#' 
#'   * Problem divided into discrete parts that can be solved concurrently
#'   * Instructions from each part execute simultaneously on different processors
#'   * An overall control/coordination mechanism is employed
#' 
#' ![https://computing.llnl.gov/tutorials/parallel_comp/](assets/parallelProblem.gif)
#' <span style="color:grey; font-size:1em;">Figure from [here](https://computing.llnl.gov/tutorials/parallel_comp/) </span>
#' 
#' 
#' 
#' ## Flynn's taxonomy
#' 
#' *  *Single Instruction, Single Data (SISD)*
#'     * No parallelization
#' *  *Single Instruction, Multiple Data (SIMD)*
#'     * Run the same code/analysis on different datasets
#'     * e.g. species or climate projections
#'     * MCMC chains
#'     * each job independent
#' * *Multiple Instruction, Single Data (MISD)*
#'     * Run different code/analyses on the same data
#' * *Multiple Instruction, Multiple Data streams (MIMD)*
#'     * Run different code/analyses on different data
#'     
#' ![http://en.wikipedia.org/wiki/Flynn%27s_taxonomy](assets/SISD.png)
#' <span style="color:grey; font-size:1em;">Figure from [here](http://en.wikipedia.org/wiki/Flynn%27s_taxonomy)</span>
#' 
#' ## Our focus: *Single Instruction, Multiple Data (SIMD)*
#' 1. Parallel loops within an R script
#'     * starts on single processor
#'     * runs looped elements on multiple 'slave' processors
#'     * returns results of all iterations to the original instance
#'     * foreach, multicore, plyr, raster
#' 2. Run many instances of same R script in parallel
#'     * need another operation to combine the results
#'     * preferable for long, complex jobs
#' 
#' ### R Packages
#' There are many R packages for parallelization, check out the CRAN Task View on [High-Performance and Parallel Computing](http://cran.r-project.org/web/views/HighPerformanceComputing.html) for an overview.  For example: 
#' 
#' * [Rmpi](http://cran.r-project.org/web/packages/Rmpi/index.html): Built on MPI (Message Passing Interface), a de facto standard in parallel computing.
#' * [snow](http://cran.r-project.org/web/packages/snow/index.html):  Simple Network of Workstations can use several standards (PVM, MPI, NWS)
#' * [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) Built in R package (since v2.14.0).
#' 
#' ---------------
#' 
#' ## Foreach Package
#' In this session we'll focus on the foreach package, which has numerous advantages including:
#' 
#'   * intuitive `for()` loop-like syntax
#'   * flexibility of choosing a parallel 'backend' for laptops through to supercomputers (using multicore, parallel, snow, Rmpi, etc.)
#'   * nice options for combining output from parallelized jobs
#' 
#' ### Documentation for foreach:
#'  - [foreach](http://cran.r-project.org/web/packages/foreach/vignettes/foreach.pdf)
#'  - [Nested Loops](http://cran.r-project.org/web/packages/foreach/vignettes/nested.pdf)
#' 
#' 
#' ### Foreach _backends_
#'  - [doParallel](http://cran.r-project.org/web/packages/doParallel/index.html) best for use on multicore machines (uses `fork` on linux/mac and `snow` on windows).
#'  - [doMPI](http://cran.r-project.org/web/packages/doMPI/vignettes/doMPI.pdf): Interface to MPI (Message-Passing Interface)
#'  - [doSNOW](http://cran.r-project.org/web/packages/doSNOW/doSNOW.pdf): Simple Network of Workstations
#' 
#' 
#' # Simple examples
#' 
#' ## _Sequential_ `for()` loop
## ------------------------------------------------------------------------
x=vector()
for(i in 1:3) x[i]=i^2
x

#' 
#' 
#' 
#' ## _Sequential_ `foreach()` loop
## ------------------------------------------------------------------------
x <- foreach(i=1:3) %do% i^2
x

#' 
#' Note that `x` is a list with one element for each iterator variable (`i`).  You can also specify a function to use to combine the outputs with `.combine`.  Let's concatenate the results into a vector with `c`.
#' 
#' ## _Sequential_ `foreach()` loop with `.combine`
## ------------------------------------------------------------------------
x <- foreach(i=1:3,.combine='c') %do% i^2
x

#' 
#' ## _Sequential_ `foreach()` loop with `.combine`
## ------------------------------------------------------------------------
x <- foreach(i=1:3,.combine='rbind') %do% i^2
x

#' 
#' So far we've only used `%do%` which only uses a single processor.
#' 
#' ## _Parallel_ `foreach()` loop
#' 
#' Before running `foreach()` in parallel, you have to register a _parallel backend_ with one of the `do` functions such as `doParallel()`. On most multicore systems, the easiest backend is typically `doParallel()`. On linux and mac, it uses `fork` system call and on Windows machines it uses `snow` backend. The nice thing is it chooses automatically for the system.
#' 
## ----,message=FALSE------------------------------------------------------
# register specified number of workers
registerDoParallel(3)
# or, reserve all all available cores 
#registerDoParallel()		
# check how many cores (workers) are registered
getDoParWorkers() 	

#' 
#' > _NOTE_ It's a good idea to use n-1 cores for processing (so you can still use your computer to do other things while the analysis is running)
#' 
## ------------------------------------------------------------------------
## run the loop
x <- foreach(i=1:3, .combine='c') %dopar% i^2
x

#' 
#' 
#' ## A slightly more complicated example
#' 
#' In this section we will:
#' 
#' 1. Generate data with known parameters
#' 2. Fit multiple chains of a bayesian regression to recover those parameters
#' 3. Compare processing times for sequential vs. parallel execution
#' 
#' Make some data
## ----makedata1-----------------------------------------------------------
n <- 100000              # number of data points
x1 <- rnorm (n)          # make up x1 covariate
x2 <- rbinom (n, 1, .5)  #make up x2 covariate
b0 <- 1.8                # set intercept (beta0)
b1 <- 1.5                # set beta1
b2 <- 2                  # set beta2
y <- rbinom (n, 1, invlogit(b0+b1*x1+b2*x2))  # simulate data with noise
data=cbind.data.frame(y=y,x1=x1,x2=x2)

#' 
#' Let's look at the data:
## ------------------------------------------------------------------------
kable(head(data),row.names = F,digits = 2)

#' 
#' 
#' Now we will specify the number of chains and fit separate bayesian GLMs using [bayesglm](http://www.inside-r.org/packages/cran/arm/docs/bayesglm) in the [ARM](http://cran.r-project.org/web/packages/arm/index.html) package.
#' 
## ----fitmodelp-----------------------------------------------------------
nchains=3

ptime <- system.time({
  result <- foreach(i=1:nchains,.combine = rbind.data.frame,.packages=c("arm")) %dopar% {
  M1=bayesglm (y ~ x1 + x2, data=data,family=binomial(link="logit"),n.iter=1e8)
  ## return parameter estimates
  cbind.data.frame(chain=i,t(coefficients(M1)))
    }
  })
ptime


#' 
#' 
#' Look at `results` object containing slope and aspect from subsampled models. There is one row per sample (`1:trials`) with columns for the estimated intercept and slope for that sample.
#' 
## ------------------------------------------------------------------------
kable(result,digits = 2)

#' 
#' So we were able to perform `r nchains` independent chains in `r ptime[3]` seconds.  Let's see how long it would have taken in sequence.
#' 
## ----fitmodel2-----------------------------------------------------------
stime <- system.time({
  result <- foreach(i=1:nchains,.combine = rbind.data.frame,.packages=c("arm")) %do% {
  M1=bayesglm (y ~ x1 + x2, data=data,family=binomial(link="logit"),n.iter=1e8)  
  ## return mean estimates
  cbind.data.frame(chain=i,t(coefficients(M1)))
    }
  })
stime


#' 
#' So we were able to run `r nchains` independent chains in `r ptime[3]` seconds when using `r getDoParWorkers()` CPUs and `r stime[3]` seconds on one CPU.  That's `r round(stime[3]/ptime[3],1)`X faster for this simple example.
#' 
#' 
#' ## Things to consider
#' ### Organizing results with `.combine`
#' Typical functions are `c` to make a vector of results, `cbind` or `cbind.data.frame` to make a table where columns correspond to different jobs, and `rbind` or `rbind.data.frame` to make a table where rows correspond to different jobs.  But you can use any function that will combine output.  
#' 
#' For example, let's extract all posteriors from the model above
## ----fitmodel------------------------------------------------------------
nchains=3
registerDoParallel(3)  	

  chains <- foreach(i=1:nchains,.combine = 'mcmc.list',.packages = c("arm","coda"),.multicombine=TRUE) %dopar% {
  M1=bayesglm(y ~ x1 + x2, data=data,family=binomial(link="logit"),n.iter=1000)
  ## extract posteriors and convert to mcmc object
  mcmc(sim(M1,500)@coef)
    }

#' 
#' That makes is relatively easy to parallelize across chains and combine the results.
## ------------------------------------------------------------------------
ggs_traceplot(ggs(chains))

#' 
#' 
#' ### Writing data to disk
#' For long-running processes, you may want to consider writing results to disk _as-you-go_ rather than waiting until the end in case of a problem (power failure, single job failure, etc.).
#' 
## ----writedata-----------------------------------------------------------
## assign target directory
td=tempdir()

  result <- foreach(i=1:nchains,.combine = rbind.data.frame,.packages=c("arm")) %dopar% {
  M1=bayesglm (y ~ x1 + x2, data=data,family=binomial(link="logit"),n.iter=1000)  
  ## return mean estimates
  results=cbind.data.frame(chain=i,t(coefficients(M1)))
  ## write results to disk
  file=paste0(td,"/results_",i,".csv")
  write.csv(results,file=file)
  return(NULL)
    }

#' 
#' That will save the result of each subprocess to disk (be careful about duplicated file names!):
## ------------------------------------------------------------------------
list.files(td,pattern="results")

#' 
#' 
#' 
#' # Spatial example
#' In this section we will:
#' 
#' 1. Generate some _spatial_ data
#' 2. Perform a moving window mean for the full area and introduce a tiling scheme
#' 3. Compare processing times for sequential vs. parallel execution
#' 
#' ## Generate Spatial Data
#' 
#' A function to generate `raster` object with spatial autocorrelation.
## ------------------------------------------------------------------------
simrast=function(nx=60,ny=60,theta=10,seed=1234){
      ## create a random raster with some spatial structure
      ## Theta is the scale of an exponential decay function.  
      ## This controls degree of autocorrelation, 
      ## values close to 1 are close to random while values near nx/4 have high autocorrelation
     r=raster(nrows=ny, ncols=nx,vals=1,xmn=-nx/2, xmx=nx/2, ymn=-ny/2, ymx=ny/2)
      names(r)="z"
      # Simulate a Gaussian random field with an exponential covariance function
      set.seed(seed)  #set a seed so everyone's maps are the same
      grid=list(x=seq(xmin(r),xmax(r)-1,by=res(r)[1]),y=seq(ymin(r),ymax(r)-1,res(r)[2]))
      obj<-Exp.image.cov(grid=grid, theta=theta, setup=TRUE)
      look<- sim.rf( obj)      
      values(r)=t(look)*10
      return(r)
      }


#' 
#' To parallelize spatial data, you often need to _tile_ the data and process each tile separately. Here is a function that will take a bounding box, tile size and generate a tiling system.  If given an `overlap` term, it will also add buffers to the tiles to reduce/eliminate edge effects, though this depends on what algorithm/model you are using.
#' 
## ------------------------------------------------------------------------
tilebuilder=function(raster,size=10,overlap=NULL){
  ## get raster extents
  xmin=xmin(raster)
  xmax=xmax(raster)
  ymin=ymin(raster)
  ymax=ymax(raster)
  xmins=c(seq(xmin,xmax-size,by=size))
  ymins=c(seq(ymin,ymax-size,by=size))
  exts=expand.grid(xmin=xmins,ymin=ymins)
  exts$ymax=exts$ymin+size
  exts$xmax=exts$xmin+size
  if(!is.null(overlap)){
  #if overlapped tiles are requested, create new columns with buffered extents
    exts$yminb=exts$ymin
    exts$xminb=exts$xmin
    exts$ymaxb=exts$ymax
    exts$xmaxb=exts$xmax
    
    t1=(exts$ymin-overlap)>=ymin
    exts$yminb[t1]=exts$ymin[t1]-overlap
    t2=exts$xmin-overlap>=xmin
    exts$xminb[t2]=exts$xmin[t2]-overlap    
    t3=exts$ymax+overlap<=ymax
    exts$ymaxb[t3]=exts$ymax[t3]+overlap
    t4=exts$xmax+overlap<=xmax
    exts$xmaxb[t4]=exts$xmax[t4]+overlap  
  }
  exts$tile=1:nrow(exts)
  return(exts)
}

#' 
#' Generate a raster and tiling system using these functions.
## ----generateraster------------------------------------------------------
r=simrast(nx=2000,ny=1000,theta = 100)
jobs=tilebuilder(r,size=1000,overlap=80)

#' 
#' Plot the raster showing the grid.
## ----plotraster,fig.height=4---------------------------------------------
ggplot(jobs)+
  geom_raster(aes(x=coordinates(r)[,1],y=coordinates(r)[,2],fill = values(r)))+ 
  scale_fill_gradient(low = 'white', high = 'blue')+
  geom_rect(mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill="transparent",lty="dashed",col="darkgreen")+
  geom_rect(aes(xmin=xminb,xmax=xmaxb,ymin=yminb,ymax=ymaxb),
            fill="transparent",col="black")+
  geom_text(aes(x=(xminb+xmax)/2,y=(yminb+ymax)/2,label=tile),size=10)+
  coord_equal()+ylab("Y")+xlab("X")

#' 
#' ## Run the `focal` function from the raster package to calculate a 3x3 moving window mean.
## ------------------------------------------------------------------------
stime2=system.time({
  r_focal1=focal(r,w=matrix(1,101,101),mean,pad=T)
  })
stime2

## plot it
gplot(r_focal1)+
  geom_raster(aes(fill = value))+ 
  scale_fill_gradient(low = 'white', high = 'blue')+
  coord_equal()+ylab("Y")+xlab("X")

#' 
#' That works great (and fast) for this little example, but as the data (or the size of the window) get larger, it can become prohibitive.  
#' 
#' ## Repeat the analysis, but parallelize using the tile system.
#' 
## ------------------------------------------------------------------------
## write a function that breaks up the original raster, computes the focal mean, then puts it back together
focal_par=function(i,raster,jobs,w=matrix(1,101,101)){
  ## identify which row in jobs to process
  t_ext=jobs[i,]
  ## crop original raster to (buffered) tile
  r2=crop(raster,extent(t_ext$xminb,t_ext$xmaxb,t_ext$yminb,t_ext$ymaxb))
  ## run moving window mean over tile
  rf=focal(r2,w=w,mean,pad=T)
  ## crop to tile
  rf2=crop(rf,extent(t_ext$xmin,t_ext$xmax,t_ext$ymin,t_ext$ymax))
  ## return the object
  return(rf2)
}


registerDoParallel(2)  	

ptime2=system.time({
  r_focal=foreach(i=1:nrow(jobs),.combine=merge,.packages=c("raster")) %dopar% focal_par(i,r,jobs)
  })


#' 
#' Are they the same?
## ------------------------------------------------------------------------
identical(r_focal,r_focal1)

#' 
#' So we were able to process the data in `r ptime2[3]` seconds when using `r getDoParWorkers()` CPUs and `r stime2[3]` seconds on one CPU.  That's `r round(stime2[3]/ptime2[3],1)`X faster for this simple example.
#' 
#' > R's Raster package can automatically parallelize some functions, check out [`clusterR`](http://www.inside-r.org/packages/cran/raster/docs/endCluster)
#' 
#' 
#' ## Other useful `foreach` parameters
#' 
#'   * `.inorder` (true/false)  results combined in the same order that they were submitted?
#'   * `.errorhandling` (stop/remove/pass)
#'   * `.packages` packages to made available to sub-processes
#'   * `.export` variables to export to sub-processes
#'   
#' 
#' # Choose your method
#' 1. Run from master process (e.g. `foreach`)
#'      - easier to implement and collect results
#'      - fragile (one failure can kill it and lose results)
#'      - clumsy for *big* jobs
#' 2. Run as separate R processes via pxargs
#'      - safer for big jobs: each job totally independent
#'      - easy to re-run incomplete submissions
#'      - compatible with qsub
#'      - forces you to have a clean processing script
#'  
#' 
#' > Each task should involve computationally-intensive work.  If the tasks are very small, it can take _longer_ to run in parallel.
#' 
#' ## Further Reading
#' 
#' * [CRAN Task View: High-Performance and Parallel Computing with R](http://cran.r-project.org/web/views/HighPerformanceComputing.html)
#' * [Simple Parallel Statistical Computing in R](www.stat.uiowa.edu/~luke/talks/uiowa03.pdf)
#' * [Parallel Computing with the R Language in a Supercomputing Environment](http://download.springer.com/static/pdf/832/chp%253A10.1007%252F978-3-642-13872-0_64.pdf?auth66=1415215123_43bf0cbf5ae8f5143b7ee309ff5e3556&ext=.pdf)
