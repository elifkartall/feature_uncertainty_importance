perf_df_wide <- do.call(rbind, performance_list) %>%
  pivot_wider(names_from = Model, values_from = c(R2_Accuracy, Coverage, Avg_Width), names_glue = "{.value}_{Model}")

normalize_minmax <- function(x){ (x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T)) }
perf_df_wide <- perf_df_wide %>% mutate(across(starts_with("Avg_Width"), ~ normalize_minmax(log1p(.))))

final_tab_data <- bind_rows(perf_df_wide, perf_df_wide %>% summarise(Dataset="GENEL ORTALAMA", across(where(is.numeric), mean)))

ft <- flextable(final_tab_data) %>%
  add_header_row(values = c("", "R² (Accuracy)", "Coverage", "Normalized Width"), colwidths = c(1, 3, 3, 3)) %>%
  set_header_labels(Dataset = "Veriseti", R2_Accuracy_LM="LM", R2_Accuracy_RF="RF", R2_Accuracy_XGBOOST="XGB",
                    Coverage_LM="LM", Coverage_RF="RF", Coverage_XGBOOST="XGB",
                    Avg_Width_LM="LM", Avg_Width_RF="RF", Avg_Width_XGBOOST="XGB") %>%
  border_remove() %>% hline_top(part="header", border=fp_border(width=2)) %>% 
  hline_bottom(part="header", border=fp_border(width=1.5)) %>% hline_bottom(part="body", border=fp_border(width=2)) %>%
  bold(part="header") %>% bold(i=nrow(final_tab_data)) %>% align(align="center", part="all") %>% autofit()

print(ft)
