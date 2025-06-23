# nf-host-rep
This repository contains a Nextflow pipeline for running the repertoire model from [Braga et al.](https://doi.org/10.1093/sysbio/syaa019).
This is part of my Bachelor's thesis project where I am exploring the how the universal probabilistic programming language [TreePPL](https://treeppl.org/) compares to the reference implementation of the model in [RevBayes](https://revbayes.github.io/).

The pipeline does the following
- Automatically generates a dataset compatible with the model
- Run a reference RevBayes implementation of the host repertoire model
- Runs one or several TreePPL implementations of the host repertoire model
- Generates a `quarto` report with simple statistics of the output.


### Running the pipeline
To run the pipeline with `podman` or `docker` simply run 
```bash
nextflow run ErikDanielsson/nf-host-rep \
    -profile local,<podman|docker>,MA \
    --runname <RUN NAME> \
    --models run_specifications/final-tppl-comp-models.csv \
    --empirical_trees empirical_data/tppl_comparision \
    --params_config sim_params/final.csv \
    --niter <MCMC-ITERATIONS> \
    --sampling_period <MCMC-SAMPLING-PERIOD> \
    --nruns <NUMBER OF RUNS> \
```
which will run the models specified in `run_specifications/final-tppl-comp-models.csv` with `<MCMC-ITERATIONS>` iterations sampling the chain every `<MCMC-ITERATIONS>` iterations.
Each model will be run `<NUMBER OF RUNS>` times to allow for multiple chain analyses.
The models are run on _Nymphalini_ trees specified in `empirical_data/tppl_comparision` directory.
Host-parasite interactions are simulated using the parameters specified in `sim_params/final.csv`.
The raw results of the simulations will be written to `simulations/<RUN NAME>`.

### Quarto report
To get a quick overview of the results of the analyses, you can generate a `quarto` report of the pipeline output with MCMC traces and ESS plots.
To do this run
```
python python_helpers/generate_report.py <RUN NAME>
```
**Installing Quarto requirements**:
In a conda environment with `pip` installed, run
`pip install --upgrade -r python_helpers/requirements.txt`.
This will install the required python packages for generating the report.
The report is rendered with `quarto` and requires the `quarto` CLI to be installed.
You furthermore need to ensure there is a Jupyter kernel detectable by the `quarto` command. 
This can be done by running
```
python -m ipykernel install --user --name <env-name> --display-name "Python (<env-name>)"
```
after activating the conda environment `<env-name>`.

### Docker images
The pipeline is run using docker images of the compiler pipeline `miking`-`miking-dppl`-`treeppl`.
These, along with docker images for the `python` and `R` dependencies, are specified in the repository [ErikDanielsson/treeppl-dockerfiles](https://github.com/ErikDanielsson/treeppl-dockerfiles).