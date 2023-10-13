import numpy as np
import pandas as pd 
import json
import os
from itertools import compress

COMPONENTS = ['local_enrichment']
NOISE_PERC = [0.2, 0.4, 0.6, 0.8, 1.0] 
GAUSS_SD = [1, 5, 10, 25, 50]
N_REP = [1] 

NETWORKS = ['pathway_commons_primary', 'pathway_commons_secondary', 'aracne_cesc', 'aracne_stad', 'aracne_read', 'aracne_paad', 'aracne_lusc', 'aracne_lihc', 'aracne_laml', 'aracne_kirp', 'aracne_kirc', 'aracne_gbm', 'aracne_blca']

with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())
for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


rule all:
    input:
        # Main pipeline
        expand(["decoupler_workflow/results/{network}/{component}/decoupler_subset_results.rds", 
            "decoupler_workflow/out_files/{network}/{component}_results.tsv", 
            "decoupler_workflow/out_files/{network}/{component}_summary_results.tsv", 
            "decoupler_workflow/results/{network}/{component}/decoupler_priori_weights.tsv",
            ], network=NETWORKS, component=COMPONENTS),
        # Gaussian noise
        expand(["decoupler_workflow/rna_gaussian_noise/expr_rna_noise_sd{gauss_sd}_rep{n_rep}.rds",
            # "decoupler_workflow/rna_gaussian_noise/rna_gaussian_noise_sd{gauss_sd}_rep{n_rep}.rds", 
            "decoupler_workflow/rna_gaussian_noise/{component}/out_files/rna_gaussian_noise_sd{gauss_sd}_rep{n_rep}_results.tsv",
            "decoupler_workflow/rna_gaussian_noise/{component}/out_files/rna_gaussian_noise_sd{gauss_sd}_rep{n_rep}_summary_results.tsv",
            ], component=COMPONENTS, gauss_sd=GAUSS_SD, n_rep=N_REP),
        # Random edges
        expand(["decoupler_workflow/random_edges/netw_random_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/random_edges/{component}/random_edges_frac{noise_perc}_rep{n_rep}.rds", 
            "decoupler_workflow/random_edges/{component}/out_files/random_edges_frac{noise_perc}_rep{n_rep}_results.tsv",
            "decoupler_workflow/random_edges/{component}/out_files/random_edges_frac{noise_perc}_rep{n_rep}_summary_results.tsv",
            ], component=COMPONENTS, noise_perc=NOISE_PERC, n_rep=N_REP),
        # Shuffle edges
        expand(["decoupler_workflow/shuffle_edges/netw_shuffle_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/shuffle_edges/{component}/shuffle_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/shuffle_edges/{component}/out_files/shuffle_edges_frac{noise_perc}_rep{n_rep}_results.tsv",
            "decoupler_workflow/shuffle_edges/{component}/out_files/shuffle_edges_frac{noise_perc}_rep{n_rep}_summary_results.tsv",
            ], component=COMPONENTS, noise_perc=NOISE_PERC, n_rep=N_REP),
        # Prune edges
        expand(["decoupler_workflow/prune_edges/netw_prune_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/prune_edges/{component}/prune_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/prune_edges/{component}/out_files/prune_edges_frac{noise_perc}_rep{n_rep}_results.tsv",
            "decoupler_workflow/prune_edges/{component}/out_files/prune_edges_frac{noise_perc}_rep{n_rep}_summary_results.tsv",
            ], component=COMPONENTS, noise_perc=NOISE_PERC, n_rep=N_REP),

include: "rules/decoupler_test_data.smk"