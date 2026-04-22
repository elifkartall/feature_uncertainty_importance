###############################################
# R2 + COVERAGE + WIDTH TABLOSU
###############################################

library(dplyr)
library(tidyr)
library(flextable)
library(officer)

perf_df <- do.call(rbind, performance_list)

final_df <- perf_df %>%
  select(Dataset, Model, R2_Accuracy, Coverage, Avg_Width) %>%
  pivot_wider(
    names_from = Model,
    values_from = c(R2_Accuracy, Coverage, Avg_Width)
  )

###############################################
# WIDTH NORMALIZE
###############################################

normalize_minmax <- function(x){
  rng <- range(x, na.rm = TRUE)
  if(rng[1] == rng[2]) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

final_df <- final_df %>%
  mutate(across(starts_with("Avg_Width"), ~ log1p(.))) %>%
  mutate(across(starts_with("Avg_Width"), normalize_minmax))

###############################################
# ORTALAMA SATIRI
###############################################

avg_row <- final_df %>%
  summarise(
    Dataset = "Ortalama",
    across(where(is.numeric), mean, na.rm = TRUE)
  )

final_df <- bind_rows(final_df, avg_row)

###############################################
# TABLO
###############################################

ft <- flextable(final_df)

ft <- add_header_row(
  ft,
  values = c("", "R-Kare", "Coverage", "Norm Width"),
  colwidths = c(1,3,3,3)
)

ft <- set_header_labels(
  ft,
  Dataset = "Veriseti",
  
  R2_Accuracy_LM = "LM",
  R2_Accuracy_RF = "RF",
  R2_Accuracy_XGBOOST = "XGB",
  
  Coverage_LM = "LM",
  Coverage_RF = "RF",
  Coverage_XGBOOST = "XGB",
  
  Avg_Width_LM = "LM",
  Avg_Width_RF = "RF",
  Avg_Width_XGBOOST = "XGB"
)

ft <- border_remove(ft)
ft <- hline_top(ft, part = "header", border = fp_border(width = 2))
ft <- hline_bottom(ft, part = "header", border = fp_border(width = 1))
ft <- hline_bottom(ft, part = "body", border = fp_border(width = 2))

ft <- bold(ft, part = "header")
ft <- bold(ft, i = nrow(final_df))

ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = 1, align = "left", part = "all")

ft <- colformat_double(ft, digits = 3)
ft <- autofit(ft)

ft
