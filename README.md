**Chicdiff: A differential caller for Capture Hi-C data** 

Chicdiff is a pipeline for identifying differential Capture Hi-C interactions between two conditions. 

Chicdiff takes advantage of Capture Hi-C parameters learned by the Chicago pipeline, and requires that the data for each replicate of each condition be processed by Chicago first. Please refer to the [Chicago bitbucket page](http://www.bitbucket.org/chicagoTeam/chicago) for details on using this pipeline.

A paper describing Chicdiff is [available here](https://doi.org/10.1093/bioinformatics/btz450): Cairns / Orchard / Malysheva & Spivakov. Chicdiff: a computational pipeline for detecting differential chromosomal interactions in Capture Hi-C data. Bioinformatics. 2019. AOP: btz450.


This repository contains the following folders:

- The Chicdiff R package     
- The ChicdiffData R data package with example Promoter Capture Hi-C datasets

Please refer to the Chicdiff R package vignette for more information.

*Installation instructions*

1. Make sure that you have R version >= 3.4.0. 

2. Install the R packages. An easy way to do this is by using functionality in devtools - run the following R code:

```{r}
install.packages("devtools")
library(devtools)
install_github("RegulatoryGenomicsGroup/chicdiff", subdir="Chicdiff")
```

Optionally, install the ChicdiffData package at the same time:

```{r}
install_github("RegulatoryGenomicsGroup/chicdiff", subdir="ChicdiffData", force=T)
```

This strategy downloads the repository multiple times. To avoid this, you can clone the github repo to a local folder and then install both packages using ``devtools::install(<path-to-package>)``.

If you encounter any problems, please [post an issue](https://github.com/RegulatoryGenomicsGroup/chicdiff/issues) or email the developers. In the email, include output from the R command ``sessionInfo()``, along with any error messages encountered.


*Contact information*

Chicdiff was developed by:

- Jonathan Cairns 
- Will Orchard
- Valeriya Malysheva
- Mikhail Spivakov ([mikhail.spivakov@lms.mrc.ac.uk](mailto:mikhail.spivakov@lms.mrc.ac.uk))

We started developing Chicdiff at the Regulatory Genomics Group, Babraham Institute, Cambridge UK. In July 2018, the group moved to MRC London Institute of Medical Sciences in London, where it is known as [Functional Gene Control](http://www.lms.mrc.ac.uk/groups/functional-gene-control) group.
