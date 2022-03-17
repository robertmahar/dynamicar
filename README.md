# Decision-theoretic adaptive SMART design

This repository contains the scripts used in developing a decision-theoretic adaptive SMART design.

# Contributors

- Robert Mahar, Melbourne School of Population and Global Health

# Overview

A simulation study. Investigates the operating characteristics of a sequential multiple assignment randomised trial (SMART) design that used Bayesian decision-theory and adaptive randomization. Note that the functions and analyses scripts contained herein were developed with this specific analysis in mind, are not intended to be used 'as is' for other analyses, and may require substantial modification and testing if used for other purposes. 

# Requirements

- Installation of R
- Installation of Stan and Rstan
- R libraries as specified in `DESCRIPTION`, `runSim.R`, `runStats.R`
- Access to the a High Performance Cluster (for efficient simulation) 
- If using a cluster, a user-created `.Renviron` containing the environment variable `PATH_TO_PKG_ON_CLUSTER` which defines the file path to the package on the users cluster: e.g. PATH_TO_PKG_ON_CLUSTER = "/home/[yourdirectory]/dynamicar"


# Workflow 

## Submit the simulation script on your favourite cluster using .

Our analysis runs 64 instances of `runSimMeerkat.R` on the HPC cluster. Each job is indexed from 1 to 64, with each job corresponding to a subset of the main simulation table (i.e. the `grid` object in `runSimMeerkat.R`) comprising 3000-ish simulation scenarios that are each run for 10 monte carlo samples. This may take a while. The output of each 64 instances is saved via `runSimMeerkat.R` with the corresponding index. The data then needs to be transferred from the HPC to the local machine for analysis. 

```{r}

qsub -t 1-64 runSimMeerkat.pbs

```

Alternatively, smaller simulations may be feasible on your local machine. If so, after modifying the simulation inputs, simply source `runSim.R` after loading the package.  

## Plot the results.

This script loads the 64 saved simulation results into memory and does plots the expected utilities in level plots, saving the outputs as a two `.pdf` files. 

```{r}

source("runStatsMeerkat.R"")

```

# End of document. 
