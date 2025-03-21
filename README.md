# nf-host-rep
This repository contains a Nextflow pipeline for running the repertoire model from [Braga et al.](https://doi.org/10.1093/sysbio/syaa019).
This is part of my bachelors thesis project where I am developing faster inference MCMC algorithms for the host repertoire model in universal probabilistic programming language [TreePPL](https://treeppl.org/).

The pipeline does the following
- Automatically generates a dataset compatible with the model
- Run a reference RevBayes implementation of the host repertoire model
- Runs one or several TreePPL implementations of the host repertoire model
- Generates a `quarto` report with simple statistics of the output.


### Running the pipeline (in progress)
To run the pipeline with `podman` simply clone the repository and run
```
nextflow run . -profile local,podman,mcmc_dists --runname <RUN NAME> --niter <iterations>
```
which will run the models specified in `run_specifications/mcmc_dists_models.csv` with `<iterations>` iterations.
The raw results of the simulations will be written to `simulations/<RUN NAME>`.

To generate a `quarto` report of the pipeline output with MCMC traces and ESS plots, run
```
python python_helpers/generate_report.py <RUN NAME>
```

### Docker images
The pipeline is run using docker images of the compiler pipeline `miking`-`miking-dppl`-`treeppl`.
These, along with docker images for the `python` and `R` dependencies, are specified in the repository [ErikDanielsson/treeppl-dockerfiles](https://github.com/ErikDanielsson/treeppl-dockerfiles).