# brentlabRnaSeqTools

This is a very helpful tutorial on making an R package:  
https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

# Installation  
In R, do the following:  

```
library(devtools)

install_github("cmatKhan/brentlabRnaSeqTools") 
# NOTE: you can set the argument upgrade = "ask" or upgrade = "never" to avoid upgrading other packages, though if you're in your base environment, you might as well upgrade if you have time

library(brentlabRnaSeqTools)
```
You will then have the functions in the R subdirectory (see the code repository). Each has a description at the top, and more documentation is on its way.
