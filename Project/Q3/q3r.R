library(limma)
library(edgeR)

normal_counts <- read.csv("train_normal_counts.csv")
meta_data <- read.csv("train_meta_data.csv")
meta_data$Simplified_class[meta_data$Simplified_class != 'Normal'] <- 'Fibrosis'
labels <- factor(meta_data$Simplified_class)

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
write.csv(top_genes_fib_vs_norm, "DEGs_Fibrosis_from_Normal.csv")
filtered_genes_fib_vs_norm <- top_genes_fib_vs_norm[1:500,]
genes_fib_vs_norm_names <- rownames(filtered_genes_fib_vs_norm)
common_genes <- intersect(rownames(normal_counts), genes_fib_vs_norm_names)
selected_normal_counts <- normal_counts[common_genes, ]

write.csv(selected_normal_counts, "subset_data1.csv")
write.csv(meta_data, "meta_data1.csv")


