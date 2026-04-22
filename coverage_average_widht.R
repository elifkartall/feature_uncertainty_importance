###############################################
# 0. LIBRARIES 
###############################################

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  OpenML,
  data.table,
  dplyr,
  RSpectra,
  tidyr,
  flextable,
  officer
)

###############################################
# 1. IDS & 2. SAFE COMPLEXITY FUNCTION (Aynı kalıyor)
###############################################

stable_ids <- c(44956, 44957, 44958, 44959, 44963,
                44964, 45012, 44971, 44977)

compute_complexity_safe <- function(id) {
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
      X_small <- X_num[, sample(p, 20), drop = FALSE]
    } else {
      X_small <- X_num
    }
    
    feature_redundancy <- if(ncol(X_small) > 1) {
      cor_mat <- cor(X_small, use = "pairwise.complete.obs")
      mean(abs(cor_mat[upper.tri(cor_mat)]), na.rm = TRUE)
    } else { NA }
    
    intrinsic_dim <- 1
    if (ncol(X_num) > 2) {
      X_scaled <- scale(X_num)
      k <- min(10, ncol(X_scaled) - 1)
      s <- RSpectra::svds(X_scaled, k = k)
      var_explained <- cumsum(s$d^2) / sum(s$d^2)
      intrinsic_dim <- which(var_explained >= 0.95)[1]
      if (is.na(intrinsic_dim)) intrinsic_dim <- k
    }
    
    local_var <- if(n > 200) var(y[sample(n, 200)], na.rm = TRUE) else var(y, na.rm = TRUE)
    
    data.frame(
      Dataset = ds_name,
      n = n,
      p = p,
      n_p_ratio = n_p_ratio,
      missing_ratio = missing_ratio,
      target_var = target_var,
      redundancy = feature_redundancy,
      int_dim = intrinsic_dim,
      local_var = local_var
    )
  }, error = function(e) data.frame(Dataset = paste0("ID_", id), n=NA, p=NA, n_p_ratio=NA, 
                                    missing_ratio=NA, target_var=NA, redundancy=NA, 
                                    int_dim=NA, local_var=NA))
}

###############################################
# 3. RUN PIPELINE & 4. PREPARE DATA
###############################################

results <- lapply(stable_ids, compute_complexity_safe)
final_df <- bind_rows(results)


apa_df <- final_df %>%
  rename(
    `Sample Size (n)` = n,
    `Features (p)` = p,
    `n/p Ratio` = n_p_ratio,
    `Missing Rate` = missing_ratio,
    `Target Var.` = target_var,
    `Redundancy` = redundancy,
    `Int. Dim.` = int_dim,
    `Local Var.` = local_var
  )

###############################################
# 5. CREATE APA TABLE 
###############################################

apa_table <- flextable(apa_df) %>%
  # Ondalık sayıları düzenle (Örn: Virgülden sonra 2 basamak)
  colformat_double(digits = 2) %>%
  colformat_int(j = c("Sample Size (n)", "Features (p)", "Int. Dim.")) %>%
  
  # APA Teması uygula (Sadece üst, alt ve başlık altı çizgileri)
  theme_apa() %>%
  
  # Genişliği otomatik ayarla
  autofit() %>%
  
  # Yazı tipini Times New Roman yap (Akademik standart)
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  
  # Tablo altı notu ekle
  add_footer_lines("")

# Tabloyu görüntüle
print(apa_table)

