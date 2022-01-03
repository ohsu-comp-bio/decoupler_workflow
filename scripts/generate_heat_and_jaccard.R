library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(ggrepel)
library(patchwork)
library(ggplotify)


args = commandArgs(trailingOnly=TRUE)

results_rds = args[1]
out_fig = args[2]


get_corr_mat <- function(df, filt = 'pathway_commons'){
  df <- df %>%
    dplyr::filter(set_name == filt) %>%
    dplyr::select(set_name, filter_crit, activity) %>%
    tidyr::unnest(activity) %>%
    select(set_name, filter_crit, id, source, statistic, score) %>%
    pivot_wider(names_from=statistic, values_from=score) %>%
    select(-set_name, -filter_crit, -id, -source) %>% drop_na()
  corr_matrix <- matrix(0, ncol(df), ncol(df))
  colnames(corr_matrix) <- colnames(df)
  rownames(corr_matrix) <- colnames(df)
  for (name_a in colnames(df)) {
    for (name_b in colnames(df)) {
      corr_matrix[name_a,name_b] <- cor(x=abs(df[[name_a]]), y=abs(df[[name_b]]), method = 'spearman')
    }
  }
  corr_matrix
}

jacc_idx <- function(a,b){
  n_inter <- length(intersect(a,b))
  n_union <- length(union(a,b))
  n_inter / n_union
}

get_jacc_mat <- function(df, filt = 'pathway_commons'){
  # Jaccard
  stat_acts <- df %>%
    dplyr::filter(set_name == filt) %>%
    select(activity) %>%
    unnest(cols=c(activity)) %>%
    select(statistic, source, id, score)

  n_top <- length(stat_acts$source %>% unique())
  n_top <- ceiling(n_top * 0.05)
  stat_acts <- stat_acts %>%
    group_by(statistic, id) %>%
    arrange(desc(abs(score))) %>%
    slice_head(n=n_top) %>%
    group_by(statistic, id) %>%
    select(-score) %>%
    nest(data=c(source)) %>%
    pivot_wider(id_cols = id, names_from = statistic, values_from = data) %>%
    column_to_rownames('id')

  jacc_matrix <- matrix(0, ncol(stat_acts), ncol(stat_acts))
  colnames(jacc_matrix) <- colnames(stat_acts)
  rownames(jacc_matrix) <- colnames(stat_acts)

  for (name_a in colnames(stat_acts)) {
    for (name_b in colnames(stat_acts)) {
      jacs <- map_dbl(rownames(stat_acts), function(sample){
        sources_a <- stat_acts[sample,name_a][[1]]$source
        sources_b <- stat_acts[sample,name_b][[1]]$source
        jacc_idx(sources_a, sources_b)
      })
      jacc_matrix[name_a,name_b] <- round(mean(jacs), 4)
    }
  }
  jacc_matrix
}

get_mat_plot <- function(mat, main = 'Spearman correlation', palette = 'Greens'){
  # Get order of clustered methods
  cor_heat <- pheatmap(mat, cluster_rows = T,
                       cluster_cols = T, silent=T)
  idxs <- cor_heat$tree_row$order
  mat <- mat[idxs,idxs]

  # Plot
  celldim <- 10
  mat_heat <- pheatmap(mat, color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name =palette)))(100),
                       display_numbers=F, cluster_rows = F, number_color='black', border_color=NA,
                       cluster_cols = T, na_col=NA, cellwidth = celldim, cellheight = celldim,
                       legend_breaks = c(0, 0.25, 0.50, 0.75, 1.0), legend=T,
                       show_rownames = T, show_colnames = T, main=main,
                       width = 4, height=4, silent=T)
  return(as.ggplot(mat_heat))
}


# Read
rna_results <- readRDS(results_rds)
rna_result <- rna_results@bench_res

# Generate matrices
rna_corr_mat <- get_corr_mat(rna_result, 'pathway_commons')
rna_jacc_mat <- get_jacc_mat(rna_result, 'pathway_commons')
print(paste0('Median rna corr: ', median(rna_corr_mat[upper.tri(rna_corr_mat)])))
print(paste0('Median rna jacc: ', median(rna_jacc_mat[upper.tri(rna_jacc_mat)])))

# Generate plots
rna_corr_plot <- get_mat_plot(rna_corr_mat, main='Spearman correlation', palette='Greens')
rna_jac_plot <- get_mat_plot(rna_jacc_mat, main='Jaccard index', palette='Reds')

# Merge together and save
pdf(file = out_fig, 
    width = 8, # The width of the plot in inches
    height = 9) # The height of the plot in inches
rna_corr_plot + rna_jac_plot +
  plot_annotation(tag_levels = 'A')
dev.off()

