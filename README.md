# brentlabRnaSeqTools

This is a very helpful tutorial on making an R package:  
https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

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
