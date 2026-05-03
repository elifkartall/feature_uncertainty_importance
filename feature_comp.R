###############################################
# 1. KURULUM VE GEREKLİ KÜTÜPHANELER
###############################################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(DALEX, predictset, ggplot2, randomForest, xgboost, OpenML, 
               dplyr, tidyr, patchwork, quantreg, quantregForest, flextable, officer)

# XGBoost Tahmin Yardımcısı
xgb_predict_final <- function(model, newdata) {
  mat <- as.matrix(sapply(newdata, as.numeric))
  return(as.numeric(predict(model, mat)))
}

# CQR Eğitim Şablonları
get_cqr_templates <- function(model_type) {
  list(
    m_lo = make_model(
      train_fun = function(x, y) {
        if(model_type %in% c("rf", "xgboost")) {
          quantregForest::quantregForest(as.data.frame(x), y, ntree = 200)
        } else {
          df <- as.data.frame(x); df$y <- y
          quantreg::rq(y ~ ., data = df, tau = 0.05)
        }
      },
      predict_fun = function(obj, x_new) {
        if(model_type %in% c("rf", "xgboost")) as.numeric(predict(obj, as.data.frame(x_new), what = 0.05))
        else as.numeric(predict(obj, as.data.frame(x_new)))
      },
      type = "regression"
    ),
    m_hi = make_model(
      train_fun = function(x, y) {
        if(model_type %in% c("rf", "xgboost")) {
          quantregForest::quantregForest(as.data.frame(x), y, ntree = 200)
        } else {
          df <- as.data.frame(x); df$y <- y
          quantreg::rq(y ~ ., data = df, tau = 0.95)
        }
      },
      predict_fun = function(obj, x_new) {
        if(model_type %in% c("rf", "xgboost")) as.numeric(predict(obj, as.data.frame(x_new), what = 0.95))
        else as.numeric(predict(obj, as.data.frame(x_new)))
      },
      type = "regression"
    )
  )
}

###############################################
# 2. ANA ANALİZ DÖNGÜSÜ
###############################################
stable_ids <- c(44956, 44957, 44958, 44959, 44963, 44964, 45012, 44971, 44977)
results_storage <- list()
performance_list <- list()

message(">>> ANALİZ BAŞLADI...")

for(id in stable_ids) {
  message(sprintf("\n--- İşleniyor ID: %s ---", id))
  
  oml <- getOMLDataSet(data.id = id)
  df_clean <- oml$data %>% 
    na.omit() %>% 
    mutate(across(where(is.character), as.factor)) %>% 
    mutate(across(where(is.factor), ~as.numeric(as.factor(.x))))
  
  colnames(df_clean) <- make.names(colnames(df_clean), unique = TRUE)
  target_raw <- if (length(oml$target.features) == 0) names(df_clean)[ncol(df_clean)] else oml$target.features
  target <- make.names(target_raw, unique = TRUE)
  
  y <- df_clean[[target]]
  x_df <- as.data.frame(df_clean[, names(df_clean) != target, drop = FALSE]) 
  
  set.seed(123)
  idx <- sample(1:nrow(x_df), 0.7 * nrow(x_df))
  x_train <- x_df[idx, , drop=FALSE]; y_train <- y[idx]
  x_test  <- x_df[-idx, , drop=FALSE]; y_test  <- y[-idx]
  
  model_results <- list()
  
  for(m in c("lm", "rf", "xgboost")) {
    message(paste("    > Model:", m))
    
    # Base Model Eğitimi
    base_mod <- if(m=="lm") {
      lm(y~., data=cbind(y=y_train, as.data.frame(x_train)))
    } else if(m=="rf") {
      randomForest(x=as.data.frame(x_train), y=y_train, ntree=200)
    } else {
      xgboost(data=as.matrix(sapply(x_train, as.numeric)), label=y_train, nrounds=100, verbose=0, params=list(objective="reg:squarederror"))
    }
    
    # CQR ve Width Hesaplama
    tm_cqr <- get_cqr_templates(m)
    cqr_res <- conformal_cqr(x = x_train, y = y_train, 
                             model_lower = tm_cqr$m_lo, model_upper = tm_cqr$m_hi, 
                             x_new = x_test, alpha = 0.1, cal_fraction = 0.5)
    widths <- interval_width(cqr_res)
    
    # PERFORMANS METRİKLERİNİ KAYDET
    y_pred_base <- if(m=="xgboost") xgb_predict_final(base_mod, x_test) else predict(base_mod, x_test)
    ss_res <- sum((y_test - y_pred_base)^2); ss_tot <- sum((y_test - mean(y_test))^2)
    r2_val <- 1 - (ss_res / ss_tot)
    
    performance_list[[paste0(id, "_", m)]] <- data.frame(
      Dataset = oml$desc$name, Model = toupper(m), R2_Accuracy = r2_val,
      Coverage = mean(y_test >= cqr_res$lower & y_test <= cqr_res$upper),
      Avg_Width = mean(widths)
    )
    
    # PFI Hesaplamaları
    exp_p <- explain(base_mod, data = as.data.frame(x_test), y = y_test, 
                     predict_function = if(m=="xgboost") xgb_predict_final else NULL, 
                     label = paste(toupper(m), "Pred"), verbose = FALSE)
    pfi_p <- model_parts(exp_p, B = 25)
    
    mod_u <- if(m=="xgboost") {
      xgboost(data=as.matrix(sapply(x_test, as.numeric)), label=widths, nrounds=100, verbose=0, params=list(objective="reg:squarederror"))
    } else if(m=="rf") {
      randomForest(x=as.data.frame(x_test), y=widths, ntree=200)
    } else {
      lm(widths~., data=cbind(widths=widths, as.data.frame(x_test)))
    }
    
    exp_u <- explain(mod_u, data = as.data.frame(x_test), y = widths, 
                     predict_function = if(m=="xgboost") xgb_predict_final else NULL, 
                     label = paste(toupper(m), "Unc"), verbose = FALSE)
    pfi_u <- model_parts(exp_u, B = 25)
    
    model_results[[m]] <- list(pfi_pred=pfi_p, pfi_unc=pfi_u)
  }
  results_storage[[as.character(id)]] <- list(dataset_name=oml$desc$name, results=model_results)
}

###############################################
# 3. PFI PLOT 
###############################################
example_id <- names(results_storage)[4]
res <- results_storage[[example_id]]$results

plot_custom_pfi <- function(pfi_obj, title_text, bar_color) {
  df_pfi <- as.data.frame(pfi_obj) %>%
    filter(!variable %in% c("_full_model_", "_baseline_")) %>%
    group_by(variable) %>% summarise(dropout_loss = mean(dropout_loss)) %>%
    arrange(desc(dropout_loss)) %>% slice_head(n = 5)
  
  ggplot(df_pfi, aes(x = reorder(variable, dropout_loss), y = dropout_loss)) +
    geom_col(fill = bar_color, width = 0.7) +
    coord_flip() + ggtitle(title_text) + theme_bw() +
    labs(x = "", y = "RMSE loss") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          panel.grid.major.y = element_blank(),
          axis.text = element_text(size = 12, color = "black"))
}

p1 <- plot_custom_pfi(res$lm$pfi_pred, "LM Tahmin", "#7F0000")
p2 <- plot_custom_pfi(res$rf$pfi_pred, "RF Tahmin", "#7F0000")
p3 <- plot_custom_pfi(res$xgboost$pfi_pred, "XGB Tahmin", "#7F0000")
p4 <- plot_custom_pfi(res$lm$pfi_unc, "LM Belirsizlik", "#808080")
p5 <- plot_custom_pfi(res$rf$pfi_unc, "RF Belirsizlik", "#808080")
p6 <- plot_custom_pfi(res$xgboost$pfi_unc, "XGB Belirsizlik", "#808080")

final_pfi_plot <- (p1 | p2 | p3) / (p4 | p5 | p6)
print(final_pfi_plot)

###############################################
# 4. ISI HARİTASI 
###############################################
plot_data_list <- list()
for(id_str in names(results_storage)) {
  for(m in names(results_storage[[id_str]]$results)) {
    p_pfi <- results_storage[[id_str]]$results[[m]]$pfi_pred
    u_pfi <- results_storage[[id_str]]$results[[m]]$pfi_unc
    
    d_p <- as.data.frame(p_pfi) %>% filter(!variable %in% c("_full_model_", "_baseline_")) %>%
      group_by(variable) %>% summarise(v_p = mean(dropout_loss))
    d_u <- as.data.frame(u_pfi) %>% filter(!variable %in% c("_full_model_", "_baseline_")) %>%
      group_by(variable) %>% summarise(v_u = mean(dropout_loss))
    
    merged <- inner_join(d_p, d_u, by = "variable")
    top_p <- merged %>% arrange(desc(v_p)) %>% slice_head(n = 5) %>% pull(variable)
    top_u <- merged %>% arrange(desc(v_u)) %>% slice_head(n = 5) %>% pull(variable)
    
    plot_data_list[[paste0(id_str, m)]] <- data.frame(
      Dataset = results_storage[[id_str]]$dataset_name, Model = toupper(m),
      Spearman = cor(merged$v_p, merged$v_u, method = "spearman"),
      Kendall = cor(merged$v_p, merged$v_u, method = "kendall"),
      Overlap = length(intersect(top_p, top_u)) / 5
    )
  }
}

plot_df <- do.call(rbind, plot_data_list) %>%
  pivot_longer(cols = c("Spearman", "Kendall", "Overlap"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(case_when(Metric=="Overlap"~"İlk 5 Örtüşme", TRUE~Metric), 
                         levels = c("Kendall", "Spearman", "İlk 5 Örtüşme")))

heatmap_plot <- ggplot(plot_df, aes(x = Model, y = Metric, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Value)), fontface = "bold", size = 8) +
  facet_wrap(~Dataset, ncol = 3) +
  scale_fill_steps2(low="#2C3E50", mid="#F7F7F7", high="#7F0000", midpoint=0, limits=c(-1,1), n.breaks=10) +
  theme_bw() + theme(strip.text = element_text(face="bold", size=16), 
                     axis.text.y = element_text(face="bold", size=18),
                     axis.text.x = element_text(face="bold", size=16), legend.position="none")

print(heatmap_plot)
