library(dplyr)
library(tibble)
library(reticulate)

args = commandArgs(trailingOnly=TRUE)
print(args)
expr_fname = args[1]
netw_fname = args[2]
netw_noise_out = args[3]
noise_perc = args[4]
n_rep = args[5]
netw_noise_out_tsv = args[6]

# Randomly sample edges from expression and randomly replace with regulons in network
expr <- readRDS(expr_fname)
netw <- readRDS(netw_fname)
edges <- rownames(expr)

# Randomly sample edges from network and regulons to replace
set.seed(42)
n_sample <- ceiling(nrow(netw) * as.numeric(noise_perc))
rand_edges <- sample(edges, n_sample, replace = TRUE)
rand_regs <- sample(1:nrow(netw), n_sample)

# Replace random edges
netw$target[rand_regs] <- rand_edges

# Save 
saveRDS(netw, netw_noise_out)

# Save for Priori
netw_priori <- netw %>% 
	filter(confidence %in% c("A", "B", "C")) %>%
	dplyr::select("tf", "target", "mor", "likelihood")
colnames(netw_priori) <- c("Regulator", "Target", "MoA", "Likelihood")

write.table(netw_priori, file = netw_noise_out_tsv, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")