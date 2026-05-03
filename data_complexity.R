###############################################
# 0. KÜTÜPHANELER
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
# 1. VERİ SETİ ID'LERİ
###############################################

stable_ids <- c(44956, 44957, 44958, 44959, 44963,
                44964, 45012, 44971, 44977)

###############################################
# 2. GÜVENLİ KARMAŞIKLIK FONKSİYONU
###############################################

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
    
    target_var <- if (is.numeric(y)) var(y, na.rm = TRUE) else NA
    
    # Tekrarlılık
    feature_redundancy <- if (ncol(X_num) > 1) {
      
      if (ncol(X_num) > 20) {
        X_small <- X_num[, sample(ncol(X_num), 20), drop = FALSE]
      } else {
        X_small <- X_num
      }
      
      cor_mat <- suppressWarnings(
        cor(X_small, use = "pairwise.complete.obs")
      )
      
      mean(abs(cor_mat[upper.tri(cor_mat)]), na.rm = TRUE)
      
    } else {
      NA
    }
    
    # İçsel Boyut
    intrinsic_dim <- 1
    
    if (ncol(X_num) > 2) {
      
      X_scaled <- scale(X_num)
      X_scaled[is.na(X_scaled)] <- 0
      
      k <- min(10, ncol(X_scaled) - 1)
      s <- RSpectra::svds(X_scaled, k = k)
      
      var_explained <- cumsum(s$d^2) / sum(s$d^2)
      intrinsic_dim <- which(var_explained >= 0.95)[1]
      
      if (is.na(intrinsic_dim)) intrinsic_dim <- k
    }
    
    # Yerel Varyans
    local_var <- if (is.numeric(y)) {
      if (n > 200) var(sample(y, 200), na.rm = TRUE)
      else var(y, na.rm = TRUE)
    } else {
      NA
    }
    
    data.frame(
      "Veri Seti" = ds_name,
      "Örnek Sayısı" = n,
      "Özellik Sayısı" = p,
      "Örnek / Özellik Oranı" = n_p_ratio,
      "Eksik Veri Oranı" = missing_ratio,
      "Hedef Varyansı" = target_var,
      "Tekrarlılık" = feature_redundancy,
      "İçsel Boyut" = intrinsic_dim,
      "Yerel Varyans" = local_var,
      check.names = FALSE
    )
    
  }, error = function(e) {
    
    data.frame(
      "Veri Seti" = paste0("ID ", id),
      "Örnek Sayısı" = NA,
      "Özellik Sayısı" = NA,
      "Örnek / Özellik Oranı" = NA,
      "Eksik Veri Oranı" = NA,
      "Hedef Varyansı" = NA,
      "Tekrarlılık" = NA,
      "İçsel Boyut" = NA,
      "Yerel Varyans" = NA,
      check.names = FALSE
    )
  })
}

###############################################
# 3. ANALİZİ ÇALIŞTIR
###############################################

results <- lapply(stable_ids, compute_complexity_safe)
final_df <- bind_rows(results)

###############################################
# 4. LOG DÖNÜŞÜMÜ
###############################################

final_df <- final_df %>%
  mutate(
    `Örnek Sayısı` = log1p(`Örnek Sayısı`),
    `Özellik Sayısı` = log1p(`Özellik Sayısı`),
    `Hedef Varyansı` = log1p(`Hedef Varyansı`),
    `Yerel Varyans` = log1p(`Yerel Varyans`)
  )

###############################################
# 5. MIN-MAX NORMALİZASYONU
###############################################

normalize_minmax <- function(x) {
  
  if (all(is.na(x))) return(x)
  
  rng <- range(x, na.rm = TRUE)
  
  if (rng[1] == rng[2]) return(rep(0, length(x)))
  
  (x - rng[1]) / (rng[2] - rng[1])
}

norm_df <- final_df %>%
  mutate(
    across(
      .cols = where(is.numeric),
      .fns  = normalize_minmax
    )
  )

###############################################
# 6. MAKALEYE HAZIR APA TABLO
###############################################

complexity_table <- flextable(norm_df) %>%
  
  colformat_double(digits = 3) %>%
  
  border_remove() %>%
  
  hline_top(
    part = "header",
    border = fp_border(width = 2)
  ) %>%
  
  hline_bottom(
    part = "header",
    border = fp_border(width = 1)
  ) %>%
  
  hline_bottom(
    part = "body",
    border = fp_border(width = 2)
  ) %>%
  
  autofit() %>%
  
  font(fontname = "Times New Roman", part = "all") %>%
  
  fontsize(size = 10, part = "all") %>%
  
  bold(part = "header") %>%
  
  # HİZALAMA: Önce her şeyi ortalıyoruz
  align(align = "center", part = "all") %>%
  
  # HİZALAMA: Sadece "Veri Seti" sütununu (1. sütun) sola yaslıyoruz
  align(j = 1, align = "left", part = "all") %>%
  
  valign(valign = "center", part = "all") %>%
  
  add_footer_lines(" ")

###############################################
# 7. TABLOYU GÖSTER
###############################################

print(complexity_table)
