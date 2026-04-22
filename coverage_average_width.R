pacman::p_load(tidyr, dplyr, flextable, officer)

###############################################
# 1. VERİYİ GENİŞ FORMATA ÇEVİR
###############################################

perf_wide <- perf_df %>%
  pivot_wider(
    names_from = Model,
    values_from = c(Coverage, Avg_Width),
    names_glue = "{.value}_{Model}"
  )

###############################################
# 2. NORMALİZASYON FONKSİYONU
###############################################

normalize_minmax <- function(x) {
  rng <- range(x, na.rm = TRUE)
  
  if (rng[1] == rng[2]) {
    return(rep(0, length(x)))
  }
  
  (x - rng[1]) / (rng[2] - rng[1])
}

###############################################
# 3. WIDTH DEĞERLERİNİ DÖNÜŞTÜR
###############################################
# Sadece Avg_Width sütunlarına uygulanır
# Önce log1p(): büyük değerleri sıkıştırır
# Sonra Min-Max: 0-1 aralığına getirir

perf_wide <- perf_wide %>%
  mutate(
    across(
      starts_with("Avg_Width"),
      ~ log1p(.)
    )
  ) %>%
  mutate(
    across(
      starts_with("Avg_Width"),
      normalize_minmax
    )
  )

###############################################
# 4. GENEL ORTALAMA SATIRI
###############################################

summary_wide <- perf_wide %>%
  summarise(
    Dataset = "Ortalama",
    across(where(is.numeric), ~ mean(., na.rm = TRUE))
  )

final_df <- bind_rows(perf_wide, summary_wide)

###############################################
# 5. TABLO OLUŞTUR
###############################################

summary_table <- flextable(final_df) %>%
  
  set_header_labels(
    Dataset = "Veriseti",
    Coverage_LM = "LM",
    Coverage_RF = "RF",
    Coverage_XGBOOST = "XGBOOST",
    Avg_Width_LM = "LM",
    Avg_Width_RF = "RF",
    Avg_Width_XGBOOST = "XGBOOST"
  ) %>%
  
  add_header_row(
    values = c("", "Kapsama (Coverage)", "Ortalama Genişlik (Width)"),
    colwidths = c(1, 3, 3)
  ) %>%
  
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
  
  colformat_double(digits = 3) %>%
  
  bold(i = nrow(final_df)) %>%
  bold(part = "header") %>%
  
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  
  autofit() %>%
  
  set_caption("")


summary_table
