# Background

This is the git repo for the work presented in Burchill & Jaeger _How reliable are standard reading time analyses? 
Hierarchical bootstrap reveals substantial power over-optimism and scale-dependent Type I error inflation_ (https://osf.io/uymfp/).

We originally developed a set of libraries and scripts to run the simulation studies presented in the paper over a remote compute cluster. This code is highly specific to the structure of the remote compute cluster that we used back then. We include it in this repo primarily for anyone who is interested in retracing the exact decisions we made for Studies 3 and 4 in the main text, as well as Auxiliary Studies 2a-d in the supplementary information (see folder `legacy_code`). For the resubmission of the paper, we additionally developed scripts that can be run---with sufficient patience---on a modern laptop. We used these new scripts to conduct additional studies (Studies 1 and 2 in the main text and Auxiliary Study 2e in the supplementary information). These scripts serve as a good starting point for anyone interested in extending our studies, or in conducting parametric or non-parametric simulation studies on their own data (see folder `start_here`).

## The directory: `start_here`

This directory is designed to be more easily replicable and is where interested parties should start if they wish to run our code themselves.

### Subdirectories

* `runs_and_scripts`: contains three such scripts, each corresponding to one of the studies we presented in the paper. Each of these scripts will run for many days. For example, using 20 cores of a 2013 MacPro with 64 GB RAM, the script for Study 2 took about 12 days to complete. 

* `og_data`: the input data used by the scripts is stored here.

* `functions`: the scripts call functions stored in the `functions` folder, and use two libraries developed as part of this project (`burchill/zplyr`, which provides convenience functions, and `burchill/cs`, which provides functionality for distributing work across cores, including remote cores). 

* `models`, `bb_files`, `mm_files`: the scripts generate interim output into these folders. The `models` folder will contain the parametric models, `bb_files` will contain the sampling indices for the BATAs, and `mm_files` will contain the model summaries for the fully fleshed out BATAs. 

* `collated_files`: this is where the "final" results files compiled from all the separate mm_files will be saved.  The git repo has included example outputs of the collated files but not of the original `mm_files` and `bb_files` (since they contained many dozens GB of data).

## The directory: `legacy_code`

The result files and original code from the remote compute cluster (and R code used to generate figures for the paper) can be found in the subdirectory `legacy_code`. We offer these files for those curious about the specifics of our previous code--getting the code here working on one's own compute cluster is almost definitely not worth it, although the code could be used as inspiration. Given that the code in `start_here` has been radically simplified, we encourage those interested to start there to understand the basics more easily.

### Subdirectories

* `legacy_data`: contains the results files and a few of the original data files necessary for certain plots.

* `plots_and_figures`: contains the R code used to generate the figures in the paper. 

* `previous_functions_and_scripts`: contains the original remote cluster R code and script files. Many of the files here may never have directly contributed to the results in the manuscript, and represent previous iterations and unpublished / draft work, etc. In order to help direct attention to the scripts that directly start the data generation of the data used in the manuscript, the relevant files have been retroactively renamed `Study_<n>_..._script.R` corresponding to the final study naming conventions in the manuscript.


