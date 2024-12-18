library(limma)
library(edgeR)

normal_counts <- read.csv("train_normal_counts.csv")
meta_data <- read.csv("train_meta_data.csv")
labels <- factor(meta_data$class)
# Create a design matrix
design <- model.matrix(~0 + labels)
colnames(design) <- levels(labels)
fit <- lmFit(normal_counts, design)
contrast.matrix <- makeContrasts(
  Fibrosis_vs_Normal = `Fibrosis` - Normal,
  levels = design
)

# Apply contrasts to the fit
fit2 <- contrasts.fit(fit, contrast.matrix)

# Empirical Bayes moderation to get p-values
fit2 <- eBayes(fit2)
# Get the top DEGs for the Fibrosis vs Normal comparison
top_genes_fib_vs_norm <- topTable(fit2, coef = "Fibrosis_vs_Normal", adjust.method = "BH", number = Inf)

# View the top DEGs
head(top_genes_fib_vs_norm)
write.csv(top_genes_fib_vs_norm, "DEGs_Fibrosis_from_Normal.csv")
filtered_genes_fib_vs_norm <- top_genes_fib_vs_norm[1:300,]
genes_fib_vs_norm_names <- rownames(filtered_genes_fib_vs_norm)
common_genes <- intersect(rownames(normal_counts), genes_fib_vs_norm_names)
selected_normal_counts <- normal_counts[common_genes, ]
write.csv(selected_normal_counts, "subset_data1.csv")
write.csv(meta_data, "meta_data1.csv")


meta_data <- subset(meta_data, Simplified_class != 'Normal')
fibrosis_normal_counts <- normal_counts[as.integer(rownames(meta_data))]
labels <- factor(meta_data$Simplified_class)
# Create a design matrix
design <- model.matrix(~0 + labels)
colnames(design) <- levels(labels)
fit <- lmFit(fibrosis_normal_counts, design)
contrast.matrix <- makeContrasts(
  AdvancedFibrosis_vs_Fibrosis = `Advanced_fibrosis` - Non_advanced_Fibrosis,
  levels = design
)

# Apply contrasts to the fit
fit2 <- contrasts.fit(fit, contrast.matrix)

# Empirical Bayes moderation to get p-values
fit2 <- eBayes(fit2)
# Get the top DEGs for the Advanced Fibrosis vs Fibrosis comparison
top_genes_adv_vs_fib <- topTable(fit2, coef = "AdvancedFibrosis_vs_Fibrosis", adjust.method = "BH", number = Inf)
write.csv(top_genes_adv_vs_fib, "DEGs_Advanced_from_NonAdvanced.csv")
filtered_genes_adv_vs_fib <- top_genes_adv_vs_fib[1:300,]
genes_adv_vs_fib_names <- rownames(filtered_genes_adv_vs_fib)
common_genes <- intersect(rownames(fibrosis_normal_counts), genes_adv_vs_fib_names)
selected_normal_counts <- fibrosis_normal_counts[common_genes, ]
write.csv(selected_normal_counts, "subset_data2.csv")
write.csv(meta_data, "meta_data2.csv")