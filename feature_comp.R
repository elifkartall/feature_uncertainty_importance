###############################################
# 1. KURULUM VE GEREKLİ KÜTÜPHANELER
###############################################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  DALEX, predictset, ggplot2, randomForest, xgboost, 
  OpenML, dplyr, tidyr, scales, flextable, quantreg, quantregForest, officer
)

xgb_predict_final <- function(model, newdata) {
  if (!is.matrix(newdata)) {
    cnames <- colnames(newdata)
    newdata <- as.matrix(newdata)
    colnames(newdata) <- cnames
  }
  return(predict(model, newdata))
}

###############################################
# 2. MODEL ÜRETİCİ FONKSİYONU
###############################################
get_model_objects <- function(x_train, y_train, model_type = "lm") {
  x_mat <- as.matrix(x_train)
  
  base_model <- switch(model_type,
                       "lm" = lm(y ~ ., data = data.frame(y = y_train, x_train)),
                       "rf" = randomForest(x_train, y_train, ntree = 500),
                       "xgboost" = xgboost(data = x_mat, label = y_train, nrounds = 100, verbose = 0,
                                           params = list(objective = "reg:squarederror"))
  )
  
  m_lo <- make_model(
    train_fun = function(x, y) {
      if(model_type == "xgboost") xgboost(data = as.matrix(x), label = y, nrounds = 100, verbose = 0, params = list(objective = "reg:quantile", alpha = 0.05))
      else if(model_type == "rf") quantregForest::quantregForest(x, y, ntree = 500)
      else quantreg::rq(y ~ ., data = data.frame(y = y, x), tau = 0.05)
    },
    predict_fun = function(obj, x_new) {
      if(model_type == "xgboost") predict(obj, as.matrix(x_new))
      else if(model_type == "rf") as.numeric(predict(obj, x_new, what = 0.05))
      else as.numeric(predict(obj, as.data.frame(x_new)))
    },
    type = "regression"
  )
  
  m_hi <- make_model(
    train_fun = function(x, y) {
      if(model_type == "xgboost") xgboost(data = as.matrix(x), label = y, nrounds = 100, verbose = 0, params = list(objective = "reg:quantile", alpha = 0.95))
      else if(model_type == "rf") quantregForest::quantregForest(x, y, ntree = 500)
      else quantreg::rq(y ~ ., data = data.frame(y = y, x), tau = 0.95)
    },
    predict_fun = function(obj, x_new) {
      if(model_type == "xgboost") predict(obj, as.matrix(x_new))
      else if(model_type == "rf") as.numeric(predict(obj, x_new, what = 0.95))
      else as.numeric(predict(obj, as.data.frame(x_new)))
    },
    type = "regression"
  )
  
  list(base = base_model, m_lo = m_lo, m_hi = m_hi)
}

###############################################
# 3. ANA HESAPLAMA DÖNGÜSÜ
###############################################
stable_ids <- c(44956, 44957, 44958, 44959, 44963,
                44964, 45012, 44971, 44977)
results_storage <- list()

cat("Analiz Başlıyor...\n")

for(id in stable_ids) {
  cat("\n--- ID:", id, "İşleniyor ---\n")
  tryCatch({
    oml <- getOMLDataSet(data.id = id)
    df_clean <- oml$data %>% na.omit() %>%
      mutate(across(where(is.character), as.factor)) %>%
      mutate(across(where(is.factor), ~as.numeric(as.factor(.x))))
    
    target <- if (length(oml$target.features) == 0) names(df_clean)[ncol(df_clean)] else oml$target.features
    y <- df_clean[[target]]; x_df <- df_clean[, names(df_clean) != target, drop = FALSE]
    
    set.seed(123); idx <- sample(1:nrow(x_df), 0.7 * nrow(x_df))
    x_train <- x_df[idx, , drop=FALSE]; y_train <- y[idx]
    x_test  <- x_df[-idx, , drop=FALSE]; y_test  <- y[-idx]
    
    model_results <- list()
    for(m in c("lm", "rf", "xgboost")) {
      cat("  Model Çalışıyor:", m, "\n")
      objs <- get_model_objects(x_train, y_train, m)
      
      # R-Kare (Accuracy)
      preds_main <- if(m == "xgboost") xgb_predict_final(objs$base, x_test) else predict(objs$base, x_test)
      r2_val <- 1 - (sum((y_test - preds_main)^2) / sum((y_test - mean(y_test))^2))
      
      # Prediction PFI
      exp_pred <- explain(objs$base, data = if(m=="xgboost") as.matrix(x_test) else x_test, y = y_test, 
                          predict_function = if(m=="xgboost") xgb_predict_final else NULL, verbose = FALSE)
      pfi_pred <- model_parts(exp_pred, B = 50)
      
      # CQR, Coverage ve Width
      cqr_res <- conformal_cqr(as.matrix(x_train), y_train, objs$m_lo, objs$m_hi, as.matrix(x_test), alpha = 0.1)
      widths <- interval_width(cqr_res)
      cov_val <- predictset::coverage(cqr_res, y_test)
      width_val <- mean(widths)
      
      # Uncertainty PFI
      mod_unc <- if(m=="xgboost") xgboost(data=as.matrix(x_test), label=widths, nrounds=100, verbose=0)
      else if(m=="rf") randomForest(x=x_test, y=widths, ntree=500)
      else lm(widths ~ ., data=x_test)
      
      exp_unc <- explain(mod_unc, data = if(m=="xgboost") as.matrix(x_test) else x_test, y = widths, 
                         predict_function = if(m=="xgboost") xgb_predict_final else NULL, verbose = FALSE)
      pfi_unc <- model_parts(exp_unc, B = 50)
      
      model_results[[m]] <- list(
        pfi_pred = pfi_pred, 
        pfi_unc = pfi_unc,
        coverage = cov_val,
        mean_width = width_val,
        accuracy_r2 = r2_val
      )
    }
    results_storage[[as.character(id)]] <- list(dataset_name = oml$desc$name, results = model_results)
  }, error = function(e) cat("HATA (ID:", id, "):", e$message, "\n"))
}

###############################################
# 4. ANALİZ VE NORMALİZASYON SÜREÇLERİ
###############################################
alignment_df <- data.frame()
performance_list <- list()

for(id_key in names(results_storage)) {
  entry <- results_storage[[id_key]]
  for(m_name in names(entry$results)) {
    res <- entry$results[[m_name]]
    
    # Isı Haritası Verisi
    df_p <- res$pfi_pred %>% filter(!variable %in% c("_baseline_", "_full_model_")) %>%
      group_by(variable) %>% summarise(v_pred = mean(dropout_loss))
    df_u <- res$pfi_unc %>% filter(!variable %in% c("_baseline_", "_full_model_")) %>%
      group_by(variable) %>% summarise(v_unc = mean(dropout_loss))
    merged <- merge(df_p, df_u, by = "variable")
    
    v1 <- merged$v_pred; v2 <- merged$v_unc
    if(sd(v1) < 1e-12) v1 <- v1 + rnorm(length(v1), 0, 1e-12)
    if(sd(v2) < 1e-12) v2 <- v2 + rnorm(length(v2), 0, 1e-12)
    
    alignment_df <- rbind(alignment_df, data.frame(
      Dataset = entry$dataset_name, Model = toupper(m_name),
      Spearman = cor(v1, v2, method = "spearman"),
      Kendall = cor(v1, v2, method = "kendall"),
      Overlap = length(intersect(merged %>% arrange(desc(v_pred)) %>% slice(1:min(5,nrow(merged))) %>% pull(variable),
                                 merged %>% arrange(desc(v_unc)) %>% slice(1:min(5,nrow(merged))) %>% pull(variable))) / min(5,nrow(merged))
    ))
    
    # Performans Tablosu Verisi
    performance_list[[length(performance_list) + 1]] <- data.frame(
      Dataset = entry$dataset_name,
      Model = toupper(m_name),
      R2_Accuracy = res$accuracy_r2,
      Coverage = res$coverage,
      Avg_Width = res$mean_width
    )
  }
}

# Geniş Formata Çevirme
perf_wide <- do.call(rbind, performance_list) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(R2_Accuracy, Coverage, Avg_Width),
    names_glue = "{.value}_{Model}"
  )

# Width Normalizasyon Fonksiyonu
normalize_minmax <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# Width Dönüşümü: Log + Min-Max
perf_wide <- perf_wide %>%
  mutate(across(starts_with("Avg_Width"), ~ log1p(.))) %>%
  mutate(across(starts_with("Avg_Width"), normalize_minmax))

# Genel Ortalama Satırı
summary_wide <- perf_wide %>%
  summarise(Dataset = "GENEL ORTALAMA", across(where(is.numeric), ~ mean(., na.rm = TRUE)))

final_table_df <- bind_rows(perf_wide, summary_wide)

###############################################
# 5. GÖRSELLEŞTİRME: ISI HARİTASI 
###############################################
plot_df <- alignment_df %>% 
  pivot_longer(cols = c("Spearman", "Kendall", "Overlap"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = case_when(
    Metric == "Overlap" ~ "İlk 5 Örtüşme", 
    Metric == "Spearman" ~ "Spearman",
    Metric == "Kendall" ~ "Kendall",
    TRUE ~ Metric
  )) %>%
  mutate(Metric = factor(Metric, levels = c("Kendall", "Spearman", "İlk 5 Örtüşme")))

heatmap_plot <- ggplot(plot_df, aes(x = Model, y = Metric, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Value)), fontface = "bold", size = 3.5) +
  facet_wrap(~Dataset, ncol = 3) +
  scale_fill_gradient2(low = "#2C3E50", mid = "#F7F7F7", high = "#7F0000", midpoint = 0, limits = c(-1, 1)) +
  theme_bw() + 
  labs(x = "", y = "", title = "") +
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold"), 
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.position = "none"
  )

print(heatmap_plot)
