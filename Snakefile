import numpy as np
import pandas as pd 
import json
import os
from itertools import compress

COMPONENTS = ['local_enrichment', 'quant_nes', 'total_enrichment']

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
            ], network=NETWORKS, component=COMPONENTS)

include: "rules/decoupler_test_data.smk"