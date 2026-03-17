# ==============================================================================
# 5. High-Dimensional Case Study (PKD Metabolomics)
# ==============================================================================
# Demonstrating scaling to 932 variables where Welch's F is unavailable.
library(dplyr)
library(F2Adonis2)

set.seed(123)

cat("\n--- Loading and Preprocessing PKD Data ---\n")

# Read raw data
pkd_raw <- read.csv("PKD_Urine_Metabolomics_LCMS_Data_110221_2groups.csv",
                    check.names = FALSE, stringsAsFactors = FALSE)

# Extract Metadata and Features
group_labels <- as.character(pkd_raw[1, -1])
group_labels <- factor(group_labels)
sample_names <- colnames(pkd_raw)[-1]
metabolite_names <- pkd_raw[-1, 1]

# Format Matrix (Samples as rows, Metabolites as columns)
Y_data <- as.matrix(pkd_raw[-1, -1])
mode(Y_data) <- "numeric"
Y_data <- t(Y_data)
rownames(Y_data) <- sample_names
colnames(Y_data) <- metabolite_names

# Missing Value Handling (>50% threshold)
missing_pct <- colMeans(is.na(Y_data))
Y_data <- Y_data[, missing_pct <= 0.5]

# Imputation (method='min', replaced with 1/5 of min positive value)
for(i in 1:ncol(Y_data)) {
  col_data <- Y_data[, i]
  if(any(is.na(col_data))) {
    min_positive <- min(col_data[col_data > 0], na.rm = TRUE)
    Y_data[is.na(col_data), i] <- min_positive / 5
  }
}

# Normalization (GroupPQN using Pool samples)
pool_idx <- grepl("Pool", rownames(Y_data), ignore.case = TRUE)
if(any(pool_idx)) {
  pool_data <- Y_data[pool_idx, , drop = FALSE]
  sample_data <- Y_data[!pool_idx, , drop = FALSE]
  sample_groups <- group_labels[!pool_idx]

  ref_spectrum <- apply(pool_data, 2, median, na.rm = TRUE)

  normalize_pqn <- function(sample_row, reference) {
    quotients <- sample_row / reference
    scaling_factor <- median(quotients[is.finite(quotients) & quotients > 0], na.rm = TRUE)
    return(sample_row / scaling_factor)
  }

  Y_normalized <- t(apply(sample_data, 1, normalize_pqn, reference = ref_spectrum))
  Y_data <- Y_normalized
  group_labels <- factor(sample_groups)
}

# Log10 Transformation and Auto-scaling
Y_data <- log10(Y_data)
AutoNorm <- function(x){(x - mean(x))/sd(x, na.rm=T)}
Y_pkd_final <- apply(Y_data, 2, AutoNorm)

cat(sprintf("Final Processed Data: %d samples x %d metabolites\n",
            nrow(Y_pkd_final), ncol(Y_pkd_final)))

# ==============================================================================
# 2. PERFORM ROBUST F2 ANALYSIS
# ==============================================================================
cat("\n--- Running Robust F2 PERMANOVA ---\n")

# Prepare Metadata Frame
pkd_meta <- data.frame(Group = group_labels)

# Run F2 with 999 permutations and Bias-Adjusted Bootstrapping
# Using Euclidean distance as data is log-transformed and scaled
pkd_f2_results <- F2_adonis2(
  formula = Y_pkd_final ~ Group,
  data = pkd_meta,
  method = "euclidean",
  permutations = 999,
  bootstrap = TRUE,
  bias.adjust = TRUE
)

print(pkd_f2_results)
