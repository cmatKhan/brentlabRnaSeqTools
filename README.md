# brentlabRnaSeqTools

This is a very helpful tutorial on making an R package:  
https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html  

this github post helped with installing bioconductor dependencies (deseq2 in this case):  
https://bioinformatics.stackexchange.com/a/3375  

and this helped with installing from github:  
https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html

Finally, here is a nice package development cheatsheet (for R):  
https://rawgit.com/rstudio/cheatsheets/master/package-development.pdf

# Installation  
In R, do the following:  
```
library(devtools)

install_github("cmatKhan/brentlabRnaSeqTools") 
# NOTE: you can set the argument upgrade = "ask" or upgrade = "never" to avoid upgrading other packages, though if you're in your base environment, you might as well upgrade if you have time

You might also try this, if you get a namespace error:  
install_github("cmatKhan/brentlabRnaSeqTools", force=TRUE)
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
