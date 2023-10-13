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
edges <- netw$target

# Randomly sample edges from network and regulons to replace
# Make sure number of samples is even
set.seed(42)
n_sample <- nrow(netw) * as.numeric(noise_perc)
if (n_sample%%2 != 0) {
	n_sample <- n_sample - 1
}
random_regs <- sample(1:nrow(netw), n_sample)

# Select half of total to shuffle
regs_1 <- random_regs[1:(n_sample / 2)]
regs_2 <- random_regs[((n_sample / 2) + 1): n_sample]

edges_1 <- netw$target[regs_1]
edges_2 <- netw$target[regs_2]

# Shuffle edges
netw$target[regs_1] <- edges_2
netw$target[regs_2] <- edges_1

# Save 
saveRDS(netw, netw_noise_out)

# Save for Priori
netw_priori <- netw %>% 
	filter(confidence %in% c("A", "B", "C")) %>%
	dplyr::select("tf", "target", "mor", "likelihood")
colnames(netw_priori) <- c("Regulator", "Target", "MoA", "Likelihood")

write.table(netw_priori, file = netw_noise_out_tsv, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")