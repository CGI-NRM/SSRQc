This R package can be installed with:

```
library(devtools)
devtools::install_github("CGI-NRM/SSRQc")
```

It has a set of dependencies that can be found in the file DESCRIPTION
or listed in the NAMESPACE file found in the repo. 

This package is unlikely to be useful for anyone not working with the
Swedish brown bear monitoring project.

The main functions is to take exports from genotyping results of SSR
and generate input data files suitable for matchning genotypes to
known individual genotypes. In addition there are tools here to merge
already matched samples with metadata and create files suitable for
imports to rovbase. 

If you want to try it out there is i vignette files capturing most of
the functionality of the package on small made up files delivered
with the package. 
