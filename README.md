# decouplerRBench snakemake workflow
`decoupleRBench snakemake workflow` evaluates the performance of biological activity 
inference methods using perturbation experiments. This workflow incorporates the regulon-enrichment package enricher. It builds on `decoupleR`, and 
more specifically the `decouple` wrapper function. As such, it requires the 
decoupleR package to be installed and it is recommended that the user is familiar 
with the [basics of decoupleR](https://saezlab.github.io/decoupleR/articles/decoupleR.html#basics-1).
The benchmark pipeline requires an input tibble with user-specified settings,
benchmark data in the form of a count table and a corresponding metadata table.

## Install
To install necessary environments and run workflow please run:

### Prepare decoupleR environment

#### Create conda environment
```
cd envs
conda env create --file=decouper_env.yaml
conda activate decoupler_env
```

#### Conda install important packages
```
conda install -c bioconda bioconductor-decoupler
conda install -c conda-forge r-devtools
R
```

#### Activate R and use devtools to download decoupleR packages
```r
library(devtools)
install_github("JEstabrook/decoupleR", ref='enrich_env') 
install_github('saezlab/decoupleRBench')
```

#### Install other necessary packages
```r
install.packages(‘doMC’)
install.packages(‘doRNG’)
install.packages(‘doParallel’)
quit()
```

### Prepare regulon-enrichment environment
```
conda create -n enrich_env
conda activate enrich_env 
conda install -c estabroj89 regulon-enrichment
```


## Evaluation
For a given `decoupleR` method, activities are inferred for each regulator and 
experiment. To evaluate their performance, all experiments are concatenated 
together to generate a response vector (whether a regulator is perturbed or not)
and a predictor vector (the regulator activities). Then, using different 
thresholds we can calculate AUROC and AUPRC for each method. Given that the true 
positive classes are limited by the regulators covered in the perturbation 
experiments, we use a downsampling strategy, where for each permutation an 
equal number of  negative classes are randomly sampled.
