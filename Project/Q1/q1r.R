library(limma)
library(edgeR)

normal_counts <- read.csv("train_normal_counts.csv")
meta_data <- read.csv("train_meta_data.csv")
labels <- factor(meta_data$Simplified_class)

# Create a design matrix
design <- model.matrix(~0 + labels)
colnames(design) <- levels(labels)
fit <- lmFit(normal_counts, design)
contrast.matrix <- makeContrasts(
    AdvancedFibrosis_vs_Normal = `Advanced_fibrosis` - Normal,
    Fibrosis_vs_Normal = Non_advanced_Fibrosis - Normal,
    AdvancedFibrosis_vs_Fibrosis = `Advanced_fibrosis` - Non_advanced_Fibrosis,
    levels = design
)

# Apply contrasts to the fit
fit2 <- contrasts.fit(fit, contrast.matrix)

# Empirical Bayes moderation to get p-values
fit2 <- eBayes(fit2)

# Get the top DEGs for the Advanced Fibrosis vs Normal comparison
top_genes_adv_vs_norm <- topTable(fit2, coef = "AdvancedFibrosis_vs_Normal", adjust.method = "BH", number = Inf)

# Get the top DEGs for the Fibrosis vs Normal comparison
top_genes_fib_vs_norm <- topTable(fit2, coef = "Fibrosis_vs_Normal", adjust.method = "BH", number = Inf)

# Get the top DEGs for the Advanced Fibrosis vs Fibrosis comparison
top_genes_adv_vs_fib <- topTable(fit2, coef = "AdvancedFibrosis_vs_Fibrosis", adjust.method = "BH", number = Inf)

filtered_genes_adv_vs_norm <- top_genes_adv_vs_norm[1:200,]
filtered_genes_fib_vs_norm <- top_genes_fib_vs_norm[1:200,]
filtered_genes_adv_vs_fib <- top_genes_adv_vs_fib[1:200,]

genes_adv_vs_norm_names <- rownames(filtered_genes_adv_vs_norm)
genes_fib_vs_norm_names <- rownames(filtered_genes_fib_vs_norm)
genes_adv_vs_fib_names <- rownames(filtered_genes_adv_vs_fib)

combined_gene_names <- unique(c(genes_adv_vs_norm_names,
                                genes_fib_vs_norm_names,
                                genes_adv_vs_fib_names))

common_genes <- intersect(rownames(normal_counts), combined_gene_names)
selected_normal_counts <- normal_counts[common_genes, ]

write.csv(selected_normal_counts, "subset_data.csv")
