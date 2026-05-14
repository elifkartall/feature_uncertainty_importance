# Explaining Feature Contributions to Predictive Performance and Predictive Uncertainty

> Although machine learning models can achieve high predictive accuracy, the factors influencing why a model makes a prediction and how reliable that prediction is may not always depend on the same variables. This study investigates the extent to which feature importance structures explaining predictive performance overlap with those explaining model uncertainty.
>
> Within the scope of the study, Linear Regression, Random Forest, and XGBoost models were evaluated on 9 regression datasets with different structural characteristics. Conformalized Quantile Regression was used for uncertainty modeling, and Permutation Feature Importance analyses were conducted for both predictive performance and prediction interval width.
>
> The obtained feature importance rankings were compared using Spearman, Kendall, and Top-5 overlap metrics. The findings indicate that, especially in complex models, prediction and uncertainty explanation structures may rely on different feature patterns.

---

# Technologies Used

<p align="left">
  <img src="https://img.shields.io/badge/R_Language-276DC3?style=for-the-badge&logo=r&logoColor=white" />
  <br>
  <img src="https://img.shields.io/badge/dplyr-1a162d?style=for-the-badge&logo=tidyverse&logoColor=white" />
  <img src="https://img.shields.io/badge/tidyr-1a162d?style=for-the-badge&logo=tidyverse&logoColor=white" />
  <img src="https://img.shields.io/badge/ggplot2-1a162d?style=for-the-badge&logo=tidyverse&logoColor=white" />
  <img src="https://img.shields.io/badge/DALEX-008080?style=for-the-badge" />
  <img src="https://img.shields.io/badge/xgboost-EA4335?style=for-the-badge" />
  <img src="https://img.shields.io/badge/randomForest-228B22?style=for-the-badge" />
</p>

---

# Project Structure

```bash
feature-uncertainty-importance/
│
├── data/                                 # Dataset loading process is included within the scripts.
│
├── scripts/
│   ├── feature_comp.R                    # Prediction vs uncertainty feature importance comparisons
│   ├── coverage_average_width.R          # Coverage and interval width analyses
│   ├── data_complexity.R                 # Data complexity metrics
│   └── korelasyon_analizi.R              
│
└── README.md
```

# Findings

## 1. Prediction vs Uncertainty Relationship Across Models

<img width="1440" height="736" alt="Gemini_Generated_Image_e4cn9ze4cn9ze4cn" src="https://github.com/user-attachments/assets/260d992c-b8f9-4d33-97f0-81433f05e5c9" />

The results demonstrate the extent to which variables important for predictive performance overlap with variables important for uncertainty estimation.

- Linear Regression models generally exhibit high Spearman and Kendall correlations.
- In Random Forest models, correlation values vary depending on the dataset.
- In XGBoost models, negative correlations are observed in some datasets.

These findings suggest that, especially in complex models, prediction and uncertainty explanation structures may rely on different feature patterns.

---

## 2. Relationship Between Data Complexity and Explanations

<img width="1345" height="784" alt="Gemini_Generated_Image_26om0i26om0i26om" src="https://github.com/user-attachments/assets/cd11c779-d9ab-4f65-84a5-fa9150f06cd2" />


This heatmap illustrates the relationship between data complexity measures and the alignment of prediction and uncertainty feature importance structures.

According to the results:

- As the number of features increases, the agreement between prediction and uncertainty importance rankings decreases.
- Correlation values tend to decrease at higher intrinsic dimensionality levels.
- Particularly in XGBoost models, feature importance structures diverge more clearly as data complexity increases.

These findings indicate that, within complex data structures, models may construct prediction and uncertainty mechanisms using different information subspaces.

---

# Conclusion

This study demonstrates that predictive performance explanations and uncertainty explanation structures in machine learning models do not always rely on the same variables.

The findings highlight that not only model accuracy but also model uncertainty should be interpreted and explained. This approach may contribute to developing more reliable and transparent artificial intelligence systems, particularly in high-risk decision-support applications.

---

# Contact
**Elif Kartal**  
📧 ds.elifkartal@gmail.com
