library(dplyr)
library(tibble)
library(reticulate)

args = commandArgs(trailingOnly=TRUE)
print(args)
netw_fname = args[1]
netw_noise_out = args[2]
noise_perc = args[3]
n_rep = args[4]
netw_noise_out_tsv = args[5]

# Randomly sample edges from expression and randomly replace with regulons in network
netw <- readRDS(netw_fname)

# Randomly sample edges from network and regulons to replace
set.seed(42)
n_sample <- ceiling(nrow(netw) * as.numeric(noise_perc))
random_regs <- sample(1:nrow(netw), n_sample)

# Remove random regulons
pruned_netw <- netw[-random_regs,]

# Save 
saveRDS(pruned_netw, netw_noise_out)

# Save for Priori
netw_priori <- netw %>% 
	filter(confidence %in% c("A", "B", "C")) %>%
	dplyr::select("tf", "target", "mor", "likelihood")
colnames(netw_priori) <- c("Regulator", "Target", "MoA", "Likelihood")

write.table(netw_priori, file = netw_noise_out_tsv, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")