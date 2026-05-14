# Explaining Feature Contributions to Predictive Performance and Predictive Uncertainty

> **Özet:**  
> Makine öğrenmesi modelleri yüksek doğruluk sağlayabilse de, modelin neden bu tahmini yaptığı ve ne kadar güvenilir olduğu her zaman aynı değişkenlerden etkilenmeyebilir. Bu çalışmada, tahmin performansını açıklayan feature importance yapıları ile model belirsizliğini açıklayan feature importance yapılarının ne ölçüde örtüştüğü incelenmiştir.
>
> Çalışma kapsamında Linear Regression, Random Forest ve XGBoost modelleri; farklı yapısal özelliklere sahip 9 regresyon veri seti üzerinde değerlendirilmiştir. Belirsizlik modellemesi için Conformalized Quantile Regression yaklaşımı kullanılmış, hem tahmin performansı hem de prediction interval width üzerinden Permutation Feature Importance analizleri gerçekleştirilmiştir.
>
> Elde edilen feature importance sıralamaları Spearman, Kendall ve Top-5 overlap metrikleri üzerinden karşılaştırılmıştır. Bulgular, özellikle karmaşık modellerde tahmin ve belirsizlik açıklamaları yapılarının farklı değişken örüntülerine dayanabildiğini göstermektedir.

---

# Kullanılan Teknolojiler

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
# Proje Yapısı

```bash
feature-uncertainty-importance/
│
├── data/                                 # Kullanılan veri setlerinün yüklenmesi aşaması scriptte yer almaktadır.
│
├── scripts/
│   ├── feature_comp.R                    # Tahmin ve belirsizlik değişken önemleri karşılaştırmaları
│   ├── coverage_average_width.R          # Coverage ve interval width analizleri
│   ├── data_complexity.R                 # Veri karmaşıklığı metrikleri
│   └── korelasyon_analizi.R              
│
└── README.md
```


# Bulgular

## 1. Modeller Bazında Prediction vs Uncertainty İlişkisi

<img width="856" height="499" alt="model_veri_metrik" src="https://github.com/user-attachments/assets/4f07245d-4df1-41e2-a57f-a74596cb2935" />


Sonuçlar tahmin performansı için önemli olan değişkenler ile belirsizlik için önemli olan değişkenlerin ne ölçüde örtüştüğünü göstermektedir.

- Linear Regression modellerinde Spearman ve Kendall korelasyonlarının çoğunlukla yüksek olduğu görülmektedir.
- Random Forest modellerinde korelasyon değerleri veri setine göre değişkenlik göstermektedir.
- XGBoost modellerinde ise bazı veri setlerinde negatif korelasyonlar gözlenmiştir.

Bu bulgular, özellikle karmaşık modellerde prediction ve uncertainty explanations yapılarının farklı feature örüntülerine dayanabileceğini göstermektedir.

---

## 2. Veri Karmaşıklığı ve Açıklama İlişkisi

<img width="856" height="499" alt="karmaşıklık_korelasyon" src="https://github.com/user-attachments/assets/c7d7dfc9-7016-48b3-a5d9-c9e49b84b52d" />


Bu heatmap, veri karmaşıklığı ölçütleri ile tahmin ve belirsizlik önemleri ilişkisi arasındaki bağlantıyı göstermektedir.

Sonuçlara göre:

- Feature sayısı arttıkça prediction ve uncertainty importance sıralamaları arasındaki uyum azalmaktadır.
- Yüksek içsel boyut seviyelerinde korelasyonların düştüğü görülmektedir.
- Özellikle XGBoost modellerinde veri karmaşıklığı arttıkça değişken önemleri yapıları belirgin şekilde ayrışmaktadır.

Bu durum, karmaşık veri yapılarında modelin tahmin ve belirsizlik süreçlerini farklı bilgi alt uzayları üzerinden oluşturabileceğini göstermektedir.

---

# Sonuç

Bu çalışma, makine öğrenmesi modellerinde tahmin performansı ile  belirsizlik açıklamaları yapılarının her zaman aynı değişkenlere dayanmadığını göstermektedir.

Elde edilen sonuçlar, yalnızca model doğruluğunun değil, model belirsizliğinin de açıklanmasının gerekli olduğunu ortaya koymaktadır. Bu yaklaşım özellikle yüksek riskli karar destek sistemlerinde daha güvenilir ve şeffaf yapay zeka uygulamaları geliştirilmesine katkı sağlayabilir.

---


---

# İletişim
**Elif Kartal**  
📧 ds.elifkartal@gmail.com

