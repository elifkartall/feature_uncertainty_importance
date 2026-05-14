# Explaining Feature Contributions to Predictive Performance and Predictive Uncertainty

> **Özet:**  
> Makine öğrenmesi modelleri yüksek doğruluk sağlayabilse de, modelin neden bu tahmini yaptığı ve ne kadar güvenilir olduğu her zaman aynı değişkenlerden etkilenmeyebilir. Bu çalışmada, tahmin performansını açıklayan feature importance yapıları ile model belirsizliğini açıklayan feature importance yapılarının ne ölçüde örtüştüğü incelenmiştir.
>
> Çalışma kapsamında Linear Regression, Random Forest ve XGBoost modelleri; farklı yapısal özelliklere sahip 9 regresyon veri seti üzerinde değerlendirilmiştir. Belirsizlik modellemesi için Conformalized Quantile Regression (CQR) yaklaşımı kullanılmış, hem tahmin performansı hem de prediction interval width üzerinden Permutation Feature Importance (PFI) analizleri gerçekleştirilmiştir.
>
> Elde edilen feature importance sıralamaları Spearman, Kendall ve Top-5 overlap metrikleri üzerinden karşılaştırılmıştır. Bulgular, özellikle karmaşık modellerde prediction ve uncertainty explanations yapılarının farklı değişken örüntülerine dayanabildiğini göstermektedir.

---

# 📂 Proje Yapısı

```bash
feature-uncertainty-importance/
│
├── data/                                 # Kullanılan veri setlerinün yüklenmesi aşaması scriptte yer almaktadır.
│
├── scripts/
│   ├── feature_comp.R                    # Prediction vs uncertainty feature importance karşılaştırmaları
│   ├── coverage_average_width.R          # Coverage ve interval width analizleri
│   ├── data_complexity.R                 # Veri karmaşıklığı metrikleri
│   └── korelasyon_analizi.R              # Spearman / Kendall analizleri
│
└── README.md
```

---

# 🧠 Çalışma Akışı

Bu çalışma dört temel aşamadan oluşmaktadır:

1. **Veri Karmaşıklığı Analizi**  
   Veri setlerinin örneklem büyüklüğü, feature sayısı, redundancy, intrinsic dimension ve local variance gibi karmaşıklık ölçümleri hesaplandı.

2. **Belirsizlik Modellemesi**  
   Conformalized Quantile Regression (CQR) yaklaşımı ile prediction interval'ları oluşturuldu.

3. **PFI Tabanlı Açıklanabilirlik Analizi**  
   - Prediction performance için feature importance
   - Prediction interval width için feature importance
   
   ayrı ayrı hesaplandı.

4. **Karşılaştırmalı Analiz**  
   Elde edilen importance sıralamaları:
   - Spearman correlation
   - Kendall correlation
   - Top-5 overlap
   
   metrikleri üzerinden karşılaştırıldı.

---

# 🛠️ Kullanılan Teknolojiler

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

# Bulgular

## 1. Modeller Bazında Prediction vs Uncertainty İlişkisi

<p align="center">
  <img src="images/model_heatmap.png" width="950"/>
</p>

Bu heatmap, prediction performance için önemli olan değişkenler ile uncertainty üretiminde önemli olan değişkenlerin ne ölçüde örtüştüğünü göstermektedir.

- Linear Regression modellerinde Spearman ve Kendall korelasyonlarının çoğunlukla yüksek olduğu görülmektedir.
- Random Forest modellerinde korelasyon değerleri veri setine göre değişkenlik göstermektedir.
- XGBoost modellerinde ise bazı veri setlerinde negatif korelasyonlar gözlenmiştir.

Bu bulgular, özellikle karmaşık modellerde prediction ve uncertainty explanations yapılarının farklı feature örüntülerine dayanabileceğini göstermektedir.

---

## 2. Veri Karmaşıklığı ve Açıklama İlişkisi

<p align="center">
  <img src="images/complexity_heatmap.png" width="950"/>
</p>

Bu heatmap, veri karmaşıklığı ölçütleri ile prediction–uncertainty importance ilişkisi arasındaki bağlantıyı göstermektedir.

Sonuçlara göre:

- Feature sayısı arttıkça prediction ve uncertainty importance sıralamaları arasındaki uyum azalmaktadır.
- Yüksek intrinsic dimension seviyelerinde korelasyonların düştüğü görülmektedir.
- Özellikle XGBoost modellerinde veri karmaşıklığı arttıkça feature importance yapıları belirgin şekilde ayrışmaktadır.

Bu durum, karmaşık veri yapılarında modelin prediction ve uncertainty süreçlerini farklı bilgi alt uzayları üzerinden oluşturabileceğini göstermektedir.

---

# 📌 Sonuç

Bu çalışma, makine öğrenmesi modellerinde prediction performance açıklamaları ile uncertainty explanations yapılarının her zaman aynı feature'lara dayanmadığını göstermektedir.

Özellikle yüksek performanslı ve karmaşık modellerde:

- Prediction explanations
- Uncertainty explanations

farklı değişken örüntülerine dayanabilmektedir.

Elde edilen sonuçlar, yalnızca model doğruluğunun değil, model belirsizliğinin de açıklanmasının gerekli olduğunu ortaya koymaktadır. Bu yaklaşım özellikle yüksek riskli karar destek sistemlerinde daha güvenilir ve şeffaf yapay zeka uygulamaları geliştirilmesine katkı sağlayabilir.

---



---

# İletişim
**Elif Kartal**  
📧 ds.elifkartal@gmail.com

