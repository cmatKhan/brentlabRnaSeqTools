[![Codecov test coverage](https://codecov.io/gh/cmatKhan/brentlabRnaSeqTools/branch/master/graph/badge.svg)](https://codecov.io/gh/cmatKhan/brentlabRnaSeqTools?branch=master)[![R-CMD-check](https://github.com/cmatKhan/brentlabRnaSeqTools/workflows/R-CMD-check/badge.svg)](https://github.com/cmatKhan/brentlabRnaSeqTools/actions)

See [Online Documentation](https://rdrr.io/github/cmatKhan/brentlabRnaSeqTools/)

all functions documentation can be accessed in RStudio with ?functionName

# Installation  
In R, do the following:  
```
library(devtools)
# remove build_vignettes to save time
install_github("cmatKhan/brentlabRnaSeqTools", build_vignettes = TRUE)
# NOTE: you can set the argument upgrade = "ask" or upgrade = "never" to avoid upgrading other packages, though if you're in your base environment, you might as well upgrade if you have time

You might also try this, if you get a namespace error:  
install_github("cmatKhan/brentlabRnaSeqTools", force=TRUE)
.rs.restartR()

# after you get the package installed, do this:
library(brentlabRnaSeqTools)
```
# uninstall
```
remove.packages("brentlabRnaSeqTools")
```

More documentation is on the way, but for now do this in the console:  
```
> ?createNinetyMinuteInductionSet
```
Scroll to the bottom and click 'index' to get an index of avaialable functions

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
