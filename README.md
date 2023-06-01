# saseR: Scalable Aberrant Splicing and Expression Retrieval

This repository contains the code to reproduce the results of the saseR manuscript. By running the worklfow.txt file, all datasets are downloaded from the respective sources, and the different workflows are ran to reproduce the figures of the manuscript. Note that some methods, can take a long time to run. Adaptations could be needed to run certain workflows on a local computer, e.g. the number of cores used to run STAR. Note that whenever 'RUV' is used as filename, this corresponds to the workflows or output of our saseR method.

## Overview of the repository

This repository contains the following main folders:

$Aberrant-expression$: This folder stores the main workflows to reproduce the results of the aberrant expression benchmarks.

$Aberrant-splicing$: This folder stores the main workflows to reproduce the results of the aberrant splicing benchmarks.

$Differential-usage$: This folder stores the main workflows to reproduce the results of the differential usage benchmarks.

$Data$: This folder stores the main R workflows to wrangle the datasets used. Additional data wrangling is done in the workflow.txt script.

$Figures$: This folder stores the scripts to reproduce the Figures published in the manuscript.

$R-Functions$: This folder stores the functions used by saseR to perform aberrant expression or splicing analyses.

Additionally, a renv.lock file is included for reproducibility of the results, and three .yml files to install different conda environments used for tha analyses. saseRPaper.yml is the main conda environment used throughout the manuscript, rsem.1.3.0 was used to perform RSEM simulations, which required different depencies to avoid errors, and Fraser.2.0.yml was used to perform the analyses with FRASER 2.0, which uses the Intron Jaccard Index. When FRASER 2.0 was benchmarked, also the renv.lock file in the Aberrant-splicing file was used, rather than the renv.lock file from the main directory.


