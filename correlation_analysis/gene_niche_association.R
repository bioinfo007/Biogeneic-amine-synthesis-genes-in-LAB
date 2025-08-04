# -----------------------------------------------------------------------------
# Usage:
#   Rscript gene_niche_association.R
#
# Description:
#   This script tests for associations between gene presence/absence and
#   isolation source using Fisher’s exact test or Chi-squared test, applies
#   Benjamini–Hochberg FDR correction, and plots the top significant genes.
#
# Required input files:
#   1. gene_presence_absence_matrix.csv
#        - Must contain a "Genome" column and gene presence/absence columns (0/1).
#   2. metadata.csv
#        - Must contain "Biosample_ID" and "Isolation_source" columns.
#
# Notes:
#   - Provide the proper **absolute or relative paths** to the input files
#     in the script where they are read:
#       gene_matrix <- read.csv("<path>/gene_presence_absence_matrix.csv", stringsAsFactors = FALSE)
#       metadata    <- read.csv("<path>/metadata.csv", stringsAsFactors = FALSE)
#   - The output files will be saved in the working directory:
#       * gene_niche_association_results.csv
#       * significant_niche_specific_genes.csv
#       * top_genes_barplot.png
#
# Dependencies:
#   install.packages(c("dplyr", "ggplot2", "tibble"))
# -----------------------------------------------------------------------------
# Required libraries -----------------------------------------------------
library(dplyr)
library(ggplot2)
library(tibble)

# 1. Load Data -----------------------------------------------------------
gene_matrix <- read.csv("gene_presence_absence_matrix.csv", stringsAsFactors = FALSE)
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# 2. Filter and match samples --------------------------------------------

# Keep only those genomes present in metadata
gene_matrix_filt <- gene_matrix %>% filter(Genome %in% metadata$Biosample_ID)

# Match metadata to gene matrix order
metadata_filt <- metadata %>%
  filter(Biosample_ID %in% gene_matrix_filt$Genome) %>%
  slice(match(gene_matrix_filt$Genome, Biosample_ID))

# Sanity check
stopifnot(all(gene_matrix_filt$Genome == metadata_filt$Biosample_ID))

cat(sprintf("Analyzing %d genomes from metadata.\n", nrow(gene_matrix_filt)))

# 3. Prepare data for testing --------------------------------------------

gene_data <- gene_matrix_filt %>% select(-Genome)
isolation <- factor(metadata_filt$Isolation_source)

# 4. Filter genes by prevalence ------------------------------------------

n_samples <- nrow(gene_data)
presence_counts <- colSums(gene_data)

min_presence <- 5
max_presence <- 0.95 * n_samples

keep_genes <- (presence_counts >= min_presence) & (presence_counts <= max_presence)
gene_data_filtered <- gene_data[, keep_genes]

cat(sprintf("Kept %d genes out of %d after filtering by prevalence.\n", sum(keep_genes), ncol(gene_data)))

# 5. Association testing -------------------------------------------------

results <- data.frame(Gene = colnames(gene_data_filtered), p_value = NA_real_)

for(i in seq_along(results$Gene)){
  gene <- results$Gene[i]
  tab <- table(gene_data_filtered[[gene]], isolation)
  
  if(all(dim(tab) == 2)){
    test <- fisher.test(tab)
    results$p_value[i] <- test$p.value
  } else {
    chi <- suppressWarnings(chisq.test(tab))
    if(any(chi$expected < 5)){
      test <- chisq.test(tab, simulate.p.value = TRUE, B = 10000)
      results$p_value[i] <- test$p.value
    } else {
      results$p_value[i] <- chi$p.value
    }
  }
}

results$adj_p_value <- p.adjust(results$p_value, method = "BH")

write.csv(results, "gene_niche_association_results.csv", row.names = FALSE)

sig_genes <- results %>% filter(adj_p_value < 0.05) %>% arrange(adj_p_value)
write.csv(sig_genes, "significant_niche_specific_genes.csv", row.names = FALSE)

cat(sprintf("Found %d significant genes at FDR < 0.05\n", nrow(sig_genes)))
if(nrow(sig_genes) == 0) stop("No significant genes found. Adjust thresholds or check data.")

# 6. Prepare barplot data -----------------------------------------------

top_genes <- head(sig_genes$Gene, 10)

barplot_data <- data.frame(
  Gene = rep(top_genes, each = length(levels(isolation))),
  Isolation = rep(levels(isolation), times = length(top_genes)),
  Count = NA_integer_
)

for(g in top_genes){
  for(n in levels(isolation)){
    barplot_data$Count[barplot_data$Gene == g & barplot_data$Isolation == n] <-
      sum(gene_data_filtered[[g]][isolation == n])
  }
}

# Function for significance stars
get_stars <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

# Stars dataframe
stars_df <- sig_genes %>%
  filter(Gene %in% top_genes) %>%
  mutate(Stars = sapply(adj_p_value, get_stars))

# Fix Count not found error
barplot_data <- as_tibble(barplot_data)

# Max bar heights for placing stars
max_counts <- barplot_data %>%
  group_by(Gene) %>%
  summarize(max_count = max(Count, na.rm = TRUE), .groups = "drop")

# Join stars and positions
stars_df <- left_join(stars_df, max_counts, by = "Gene")

# 7. Plot barplot with stars ---------------------------------------------

p <- ggplot(barplot_data, aes(x = Gene, y = Count, fill = Isolation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(title = "Presence Counts of Top Significant Genes by Isolation Source", y = "Count of Presence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = stars_df,
            aes(x = Gene, y = max_count + 5, label = Stars),
            inherit.aes = FALSE,
            size = 6,
            vjust = 0)

print(p)

# 8. Save plot -----------------------------------------------------------

ggsave("top_genes_barplot.png", p, width = 10, height = 6)

cat("Barplot with significance stars saved as 'top_genes_barplot.png'.\n")
