[![Codecov test coverage](https://codecov.io/gh/cmatKhan/brentlabRnaSeqTools/branch/master/graph/badge.svg)](https://codecov.io/gh/cmatKhan/brentlabRnaSeqTools?branch=master)[![R-CMD-check](https://github.com/cmatKhan/brentlabRnaSeqTools/workflows/R-CMD-check/badge.svg)](https://github.com/cmatKhan/brentlabRnaSeqTools/actions)

# Installation and updating 
The following will both install, and update if there are changes in the repository.
```
library(devtools)
# remove build_vignettes to save time
install_github("cmatKhan/brentlabRnaSeqTools", build_vignettes = TRUE, dependencies = TRUE)

# after you get the package installed, do this:
library(brentlabRnaSeqTools)

# if you think there are changes, but install_github disagrees, try using the argument force = TRUE
```
I have also installed this on my htcf cluster profile like so:
```
ml miniconda # note: this will definitely load, and most likely work as expected. But it does not come with a promise. It is a cluster module I wrote. If you have issues which you suspect to be a conda problem, I suggest that you install a version of miniconda in your home profile. It will be easier to address any conda related issues that way.

conda install -n brentlabRnaSeqTools # or whatever you want to call your env name

conda install r r-essentials libpq

$ R

> install.packages(devtools)
# YOU HAVE TO DO THIS! do not update RSQLite (as of 20210702 there is an install error in the boost/c++ package which is a dependency. You do not need to worry about this when you're installing)
> remotes::install_version("RSQLite", version = "2.2.5")
> install_github("cmatKhan/brentlabRnaSeqTools")
```
See the bamtools vignette for examples of how to use the functions to examine bam files in an Rscript that you could run with SLURM

# uninstall
```
remove.packages("brentlabRnaSeqTools")
```

# Documentation

if you used the `build_vignettes = TRUE` argument in install_github, then you can check the available vignettes like so:
```
> browseVignettes("brentlabRnaSeqTools")
```
Additionally, all functions and data variables are documented. If you type the following into your console:
```
> brentlabRnaSeqTools::
```
and hit tab, a list of functions and variables will appear. Select any, and place a question mark in the beginning to
view the documentation:
```
> ?brentlabRnaSeqTools::getMetadata
```
There is also online documentation here, though the site is so full of advertisements it is hard to use:
See [Online Documentation](https://rdrr.io/github/cmatKhan/brentlabRnaSeqTools/)

# issues  
please do post issues to the issues tab. Please include the full error code and the command/context that lead to the error

# to contribute  
1. fork the repo
2. develop in a branch
3. create a pull request for the branch


This is the featureCounts/subreads homepage. In particular, has a good example of how to make mean/variance graph with voom
http://bioinf.wehi.edu.au/RNAseqCaseStudy/

# TODOs
 - read more about packrat, add some instructions on how to use
 - update R and dependencies to R version 4

# brentlabRnaSeqTools

This is a very helpful tutorial on making an R package:  
https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html  

this github post helped with installing bioconductor dependencies (deseq2 in this case):  
https://bioinformatics.stackexchange.com/a/3375  

and this helped with installing from github:  
https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html

Finally, here is a nice package development cheatsheet (for R):  
https://rawgit.com/rstudio/cheatsheets/master/package-development.pdf
