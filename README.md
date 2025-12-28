# pELAstic pelagic limnocorral amplicon sequencing

![Idealized Limnocorrals in Lake 378](figures/limnocorral.png)


This repository contains all the code used to generate the figures in the paper. The `docker_runner.sh` script pulls a docker image from my docker hub, which handles all the dependencies. The container opens into a R terminal, but you can change this behaviour in the `sh` script. The `data` directory contains the preprocessed phyloseq objects, the `output` directory contains output from several scripts, the `R` directory contains the scripts used to generate the figures, which are saved in to the `figures` directory. The R scripts `08_ancombc_stability.R` and `10_random_forests.R` are meant to be ran on a server. The ancombc script may take overnight to run and the Random forests script may take an hour to run. 

