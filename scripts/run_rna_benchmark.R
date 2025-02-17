library(decoupleR)
library(decoupleRBench)
library(dplyr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
print(args)
expr_fname = args[1]
meta_fname = args[2]
network = args[3]
results_out = args[4]
enrich_slot = args[5]
enr_weights = args[6]

# Define data, metadata and network path
seed <- 1
nproc = 4

netw_fname = sprintf("input_data/%s_decoupler.rds", network)

stats_list = list(c('aucell','wmean','wsum','ulm','viper','gsva','ora','fgsea','udt','mdt','enricher'))
opts_list <- list(list(
  udt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, min_n = 20, seed=seed),
  mdt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, trees = 10, min_n = 20,
             nproc = nproc, seed=seed),
  aucell = list(nproc=nproc, seed=seed),
  wmean = list(times=100, sparse=TRUE, randomize_type = "rows", seed=seed),
  wsum = list(times=100, sparse=TRUE, randomize_type = "rows", seed=seed),
  ulm = list(sparse=FALSE, center=FALSE, na.rm=FALSE),
  viper = list(verbose=FALSE, pleiotropy=T, eset.filter=F),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list(n_up=300, n_bottom=300, n_background=20000, with_ties = TRUE),
  fgsea = list(times=100, nproc=nproc, seed=seed),
  enricher = list(scaler_type="robust", minsize=5, enr_type=enrich_slot, enr_weights=enr_weights)
))

# Design
design <- tibble(
  set_name = 'pathway_commons', # name of the set resource
  bench_name = "dbd", # name of the benchmark data
  stats_list = stats_list,
  opts_list = opts_list,
  bexpr_loc = expr_fname, # benchmark data location
  bmeta_loc = meta_fname, # metadata location
  source_loc = netw_fname, # set source location
  source_col = "tf", # source name of the gene set source
  target_col = "target", # target name of the set source
  filter_col = "confidence", # column by which we wish to filter
  filter_crit = list(c('A','B','C')) # criteria by which we wish to filter
)

# Run benchmark
result <- run_benchmark(
  .design = design, # provide input tibble
  .minsize = 5, # filter gene sets with size < 10
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = FALSE, # silently run the pipeline
  .downsample_pr = TRUE, # downsample TNs for precision-recall curve
  .downsample_roc = TRUE, # downsample TNs for ROC
  .downsample_times = 100, # downsampling iterations
  .url_bool = FALSE # whether to load from url
)

print(result)
# Save result
saveRDS(result,results_out) 
