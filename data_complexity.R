###############################################
# 0. LIBRARIES
###############################################

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  OpenML,
  data.table,
  dplyr,
  RSpectra,
  ggplot2,
  tidyr,
  forcats,
  dendextend,
  RColorBrewer,
  scales
)

###############################################
# 1. IDS
###############################################

stable_ids <- c(44956, 44957, 44958, 44959, 44963,
                44964, 45012, 44971, 44977)

###############################################
# 2. COMPLEXITY FUNCTION
###############################################

compute_complexity_safe <- function(id) {
  
  cat("Processing:", id, "\n")
  
  tryCatch({
    
    ds <- getOMLDataSet(data.id = id)
    ds_name <- ds$desc$name
    data <- as.data.frame(ds$data)
    
    target <- ds$target.features
    if (is.null(target)) target <- names(data)[ncol(data)]
    
    y <- data[[target]]
    X <- data[, setdiff(names(data), target), drop = FALSE]
    X_num <- X[, sapply(X, is.numeric), drop = FALSE]
    
    n <- nrow(data)
    p <- ncol(X_num)
    
    n_p_ratio <- n / (p + 1)
    missing_ratio <- mean(is.na(data))
    target_var <- var(y, na.rm = TRUE)
    
    if (p > 20) {
      cols <- sample(p, 20)
      X_small <- X_num[, cols, drop = FALSE]
    } else {
      X_small <- X_num
    }
    
    if (ncol(X_small) > 1) {
      cor_mat <- cor(X_small, use = "pairwise.complete.obs")
      feature_redundancy <- mean(abs(cor_mat[upper.tri(cor_mat)]), na.rm = TRUE)
    } else {
      feature_redundancy <- NA
    }
    
    intrinsic_dim <- 1
    
    if (ncol(X_num) > 2) {
      X_scaled <- scale(X_num)
      k <- min(10, ncol(X_scaled) - 1)
      s <- RSpectra::svds(X_scaled, k = k)
      var_explained <- cumsum(s$d^2) / sum(s$d^2)
      intrinsic_dim <- which(var_explained >= 0.95)[1]
      if (is.na(intrinsic_dim)) intrinsic_dim <- k
    }
    
    if (n > 200) {
      idx <- sample(n, 200)
      local_var <- var(y[idx], na.rm = TRUE)
    } else {
      local_var <- var(y, na.rm = TRUE)
    }
    
    data.frame(
      dataset_name = ds_name,
      dataset_id = id,
      n = n,
      p = p,
      n_p_ratio = n_p_ratio,
      missing_ratio = missing_ratio,
      target_variance = target_var,
      feature_redundancy = feature_redundancy,
      intrinsic_dim_95 = intrinsic_dim,
      local_target_variance = local_var
    )
    
  }, error = function(e) {
    
    data.frame(
      dataset_name = paste0("ID_", id),
      dataset_id = id,
      n = NA, p = NA, n_p_ratio = NA,
      missing_ratio = NA,
      target_variance = NA,
      feature_redundancy = NA,
      intrinsic_dim_95 = NA,
      local_target_variance = NA
    )
  })
}

###############################################
# 3. RUN PIPELINE
###############################################

results <- lapply(stable_ids, compute_complexity_safe)
final_df <- bind_rows(results)

###############################################
# 4. PREPARE MATRIX
###############################################

df_mat <- final_df

rownames(df_mat) <- make.unique(df_mat$dataset_name)
df_mat$dataset_name <- NULL
df_mat$dataset_id   <- NULL

colnames(df_mat) <- c(
  "Sample Size",
  "Features",
  "n/p Ratio",
  "Missing Rate",
  "Target Variance",
  "Feature Redundancy",
  "Intrinsic Dim.",
  "Local Variance"
)

df_mat <- df_mat[, apply(df_mat, 2, function(x) sd(x, na.rm = TRUE) > 0)]

###############################################
# 5. Z-SCORE NORMALIZATION
###############################################

df_scaled <- scale(df_mat)

df_scaled[is.na(df_scaled)] <- 0
df_scaled[is.infinite(df_scaled)] <- 0




row_clust <- hclust(dist(df_scaled), method = "ward.D2")
col_clust <- hclust(dist(t(df_scaled)), method = "ward.D2")

df_scaled <- df_scaled[row_clust$order, col_clust$order]



plot_df <- as.data.frame(df_scaled)
plot_df$Dataset <- rownames(plot_df)

plot_long <- pivot_longer(
  plot_df,
  cols = -Dataset,
  names_to = "Metric",
  values_to = "Zscore"
)

plot_long$Dataset <- factor(
  plot_long$Dataset,
  levels = rev(unique(plot_long$Dataset))
)

plot_long$Metric <- factor(
  plot_long$Metric,
  levels = unique(plot_long$Metric)
)



ggplot(plot_long, aes(Metric, Dataset, fill = Zscore)) +
  
  geom_tile(color = "white", linewidth = 0.6) +
  
  geom_text(
    aes(label = round(Zscore, 2)),
    size = 3.2,
    family = "serif",
    color = "black"
  ) +
  
  scale_fill_gradient2(
    low = "#2C3E50",
    mid = "#F7F7F7",
    high = "#7F0000",
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish,
    name = "Z-score"
  ) +
  
  labs(
    title = "",
    subtitle = "",
    x = NULL,
    y = NULL,
    caption = ""
  ) +
  
  theme_minimal(base_family = "serif", base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5,
      margin = ggplot2::margin(b = 8)
    ),
    
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5,
      color = "grey30",
      margin = ggplot2::margin(b = 14)
    ),
    
    plot.margin = ggplot2::margin(15, 20, 15, 15)
  )