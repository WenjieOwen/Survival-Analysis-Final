---
title: "descriptive_and_exploratory_analysis"
author: "Yixin Zheng"
date: '2025-11-25'
output:
  pdf_document: default
  html_document: default
---



# Data Preparation
## Data Wrangling (wenjie)

``` r
# Import data & data cleaning 
cirrhosis <- read.csv("data/cirrhosis.csv") |> clean_names()

# Status convert to event indicator (1 = death, 0 = censored)
cirrhosis <- cirrhosis |>
  mutate(
    event = case_when(
      status == "D" ~ 1,
      status %in% c("C", "CL") ~ 0,
      TRUE ~ NA_real_
    ),
    sex = factor(sex),
    drug = factor(drug),
    ascites = factor(ascites),
    hepatomegaly = factor(hepatomegaly),
    spiders = factor(spiders),
    edema = factor(edema, levels = c("N", "S", "Y"), ordered = TRUE),
    stage = factor(stage) # Yixin: "Stage: histologic stage of disease (1, 2, 3, or 4)"
  )

cirrhosis <- cirrhosis |>
  mutate(age_years = age / 365.25)

skim(cirrhosis)
```


Table: Data summary

|                         |          |
|:------------------------|:---------|
|Name                     |cirrhosis |
|Number of rows           |418       |
|Number of columns        |22        |
|_______________________  |          |
|Column type frequency:   |          |
|character                |1         |
|factor                   |7         |
|numeric                  |14        |
|________________________ |          |
|Group variables          |None      |


**Variable type: character**

|skim_variable | n_missing| complete_rate| min| max| empty| n_unique| whitespace|
|:-------------|---------:|-------------:|---:|---:|-----:|--------:|----------:|
|status        |         0|             1|   1|   2|     0|        3|          0|


**Variable type: factor**

|skim_variable | n_missing| complete_rate|ordered | n_unique|top_counts                   |
|:-------------|---------:|-------------:|:-------|--------:|:----------------------------|
|drug          |       106|          0.75|FALSE   |        2|D-p: 158, Pla: 154           |
|sex           |         0|          1.00|FALSE   |        2|F: 374, M: 44                |
|ascites       |       106|          0.75|FALSE   |        2|N: 288, Y: 24                |
|hepatomegaly  |       106|          0.75|FALSE   |        2|Y: 160, N: 152               |
|spiders       |       106|          0.75|FALSE   |        2|N: 222, Y: 90                |
|edema         |         0|          1.00|TRUE    |        3|N: 354, S: 44, Y: 20         |
|stage         |         6|          0.99|FALSE   |        4|3: 155, 4: 144, 2: 92, 1: 21 |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate|     mean|      sd|      p0|      p25|      p50|      p75|     p100|hist  |
|:-------------|---------:|-------------:|--------:|-------:|-------:|--------:|--------:|--------:|--------:|:-----|
|id            |         0|          1.00|   209.50|  120.81|    1.00|   105.25|   209.50|   313.75|   418.00|▇▇▇▇▇ |
|n_days        |         0|          1.00|  1917.78| 1104.67|   41.00|  1092.75|  1730.00|  2613.50|  4795.00|▅▇▆▃▂ |
|age           |         0|          1.00| 18533.35| 3815.85| 9598.00| 15644.50| 18628.00| 21272.50| 28650.00|▂▆▇▅▁ |
|bilirubin     |         0|          1.00|     3.22|    4.41|    0.30|     0.80|     1.40|     3.40|    28.00|▇▁▁▁▁ |
|cholesterol   |       134|          0.68|   369.51|  231.94|  120.00|   249.50|   309.50|   400.00|  1775.00|▇▁▁▁▁ |
|albumin       |         0|          1.00|     3.50|    0.42|    1.96|     3.24|     3.53|     3.77|     4.64|▁▂▇▇▁ |
|copper        |       108|          0.74|    97.65|   85.61|    4.00|    41.25|    73.00|   123.00|   588.00|▇▂▁▁▁ |
|alk_phos      |       106|          0.75|  1982.66| 2140.39|  289.00|   871.50|  1259.00|  1980.00| 13862.40|▇▁▁▁▁ |
|sgot          |       106|          0.75|   122.56|   56.70|   26.35|    80.60|   114.70|   151.90|   457.25|▇▇▁▁▁ |
|tryglicerides |       136|          0.67|   124.70|   65.15|   33.00|    84.25|   108.00|   151.00|   598.00|▇▂▁▁▁ |
|platelets     |        11|          0.97|   257.02|   98.33|   62.00|   188.50|   251.00|   318.00|   721.00|▅▇▃▁▁ |
|prothrombin   |         2|          1.00|    10.73|    1.02|    9.00|    10.00|    10.60|    11.10|    18.00|▇▅▁▁▁ |
|event         |         0|          1.00|     0.39|    0.49|    0.00|     0.00|     0.00|     1.00|     1.00|▇▁▁▁▅ |
|age_years     |         0|          1.00|    50.74|   10.45|   26.28|    42.83|    51.00|    58.24|    78.44|▂▆▇▅▁ |

## Data Cleaning (yixin)

``` r
# Yixin: restrict to randomized trial patients (non-missing drug)
cirr_trial <- cirrhosis |> filter(!is.na(drug))

# Yixin: store counts for text/reporting
n_total <- nrow(cirrhosis)
n_trial <- nrow(cirr_trial)
n_nonrand <- n_total - n_trial

# Yixin: baseline variables that will go into the tables/tests
baseline_vars <- c("age_years", "sex", "drug", "bilirubin", "albumin", 
                   "prothrombin", "ascites", "hepatomegaly", "spiders", "edema")

# Yixin: complete-case dataset for baseline comparisons
cirr_complete <- cirr_trial |> drop_na(all_of(baseline_vars))
```

# Descriptive and Exploratory Analysis
We assessed baseline balance between the Placebo and D-penicillamine groups using summaries of key demographic and clinical variables.
## Basic EDA (wenjie)

``` r
cirrhosis |>
  count(status) |>
  ggplot(aes(x = status, y = n, fill = status)) +
  geom_col() +
  labs(title = "Distribution of Survival Status")
```

![](descriptive_and_exploratory_analysis_files/figure-latex/wenjie_basic_EDA-1.pdf)<!-- --> 

``` r
ggplot(cirrhosis, aes(x = status, y = age, fill = status)) +
  geom_boxplot() +
  labs(title = "Age Distribution by Survival Status", y = "Age (in days)", x = "Status") +
  theme_minimal()
```

![](descriptive_and_exploratory_analysis_files/figure-latex/wenjie_basic_EDA-2.pdf)<!-- --> 

``` r
cirrhosis |>
  pivot_longer(cols = c(albumin, bilirubin), names_to = "marker", values_to = "value") |>
  ggplot(aes(x = status, y = value, fill = status)) +
  geom_boxplot() +
  facet_wrap(~marker, scales = "free_y") +
  labs(title = "Key Lab Markers by Survival Status") +
  theme_minimal()
```

![](descriptive_and_exploratory_analysis_files/figure-latex/wenjie_basic_EDA-3.pdf)<!-- --> 


## Continuous Baseline Variables (yixin)
For continuous baseline variables, 
we compared the two treatment groups using both:
Two-sample t-tests (for differences in mean values), 
and Wilcoxon rank-sum tests (nonparametric robustness checks).
The corresponding hypotheses were:
$$
\begin{aligned}
H_{0}:~& \mu_{\text{Placebo}} = \mu_{\text{Dpen}} \\
H_{1}:~& \mu_{\text{Placebo}} \neq \mu_{\text{Dpen}}
\end{aligned}
$$


``` r
# Yixin: baseline tables + hypothesis tests by treatment group
# randomized patients only (cirr_trial), using cirr_complete for continuous vars.

# Yixin: we choose clinically relevant continuous baseline variables with minimal missingness.
# Including variables like copper or cholesterol would greatly reduce complete-case sample size ,
# so we follow the Mayo PBC literature (Dickson et al.) and summarize only key prognostic markers.
cont_vars <- c("age_years", "bilirubin", "albumin", "prothrombin")

cont_summary <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x   <- sd(x, na.rm = TRUE)
  med_x  <- median(x, na.rm = TRUE)
  q1_x   <- quantile(x, 0.25, na.rm = TRUE)
  q3_x   <- quantile(x, 0.75, na.rm = TRUE)
  
  list(
    mean = mean_x,
    sd   = sd_x,
    med  = med_x,
    q1   = q1_x,
    q3   = q3_x
  )
  }

baseline_cont <- lapply(cont_vars, function(v) {
  x_placebo <- cirr_complete |> # Yixin: extract values by treatment group from complete-case dataset
    filter(drug == "Placebo") |>
    pull(!!sym(v))
  
  x_dpen <- cirr_complete |>
    filter(drug == "D-penicillamine") |>
    pull(!!sym(v))
  
  s0 <- cont_summary(x_placebo)  # Yixin: Placebo summary
  s1 <- cont_summary(x_dpen)     # Yixin: D-penicillamine summary
  
# Yixin: two-sample t-test and Wilcoxon rank-sum test (D-penicillamine vs Placebo)
  t_res <- t.test(x_dpen, x_placebo, var.equal = FALSE)
  w_res <- wilcox.test(x_dpen, x_placebo, exact = FALSE)
  
  tibble(
    Variable = v,
    Placebo_mean_SD = sprintf("%.2f (%.2f)", s0$mean, s0$sd),
    Dpen_mean_SD = sprintf("%.2f (%.2f)", s1$mean, s1$sd),
    Placebo_median_IQR = sprintf("%.2f (%.2f, %.2f)", s0$med, s0$q1, s0$q3),
    Dpen_median_IQR = sprintf("%.2f (%.2f, %.2f)", s1$med, s1$q1, s1$q3),
    t_test_p_value = t_res$p.value,
    wilcoxon_rank_sum_p_value = w_res$p.value)
  }) |>
  bind_rows()

# Yixin: replace variable codes with full names for the table
baseline_cont <- baseline_cont |>
  mutate(
    Variable = dplyr::recode(
      Variable,
      age_years   = "Age (years)",
      bilirubin   = "Serum bilirubin (mg/dL)",
      albumin     = "Albumin (g/dL)",
      prothrombin = "Prothrombin time (seconds)"
      )
    )

# Yixin: column labels for the table header
colnames(baseline_cont) <- c(
  "Baseline variable",
  "Placebo mean (SD)",
  "D-penicillamine mean (SD)",
  "Placebo median (IQR)",
  "D-penicillamine median (IQR)",
  "t-test p-value",
  "Wilcoxon p-value"
  )

# Yixin: final baseline table for continuous variables
knitr::kable(
  baseline_cont,
  digits = 3,
  caption = "Baseline continuous variables by treatment group (randomized complete cases)"
  )
```



Table: Baseline continuous variables by treatment group (randomized complete cases)

|Baseline variable          |Placebo mean (SD) |D-penicillamine mean (SD) |Placebo median (IQR) |D-penicillamine median (IQR) | t-test p-value| Wilcoxon p-value|
|:--------------------------|:-----------------|:-------------------------|:--------------------|:----------------------------|--------------:|----------------:|
|Age (years)                |48.58 (9.96)      |51.42 (11.01)             |48.11 (41.43, 55.80) |51.93 (42.98, 58.90)         |          0.018|            0.020|
|Serum bilirubin (mg/dL)    |3.65 (5.28)       |2.87 (3.63)               |1.30 (0.72, 3.60)    |1.40 (0.80, 3.20)            |          0.133|            0.842|
|Albumin (g/dL)             |3.52 (0.40)       |3.52 (0.44)               |3.54 (3.34, 3.78)    |3.56 (3.21, 3.83)            |          0.874|            0.951|
|Prothrombin time (seconds) |10.80 (1.14)      |10.65 (0.85)              |10.60 (10.00, 11.40) |10.60 (10.03, 11.00)         |          0.199|            0.588|

From table above, we see that, for these 4 key continuous baseline variables, only age has a statistically significant difference between treatment groups (t-test $p \approx 0.018$; Wilcoxon $p \approx 0.020$), 
meaning that patients in the D-penicillamine group were slightly older at baseline. 
For bilirubin, albumin, and prothrombin time, we did not reject the null hypothesis as the p-values obatined were large. Overall, continuous variables are pretty balanced between treatment groups beside age. 

## Categorical Baseline Variables (yixin)
For categorical baseline variables, 
we evaluated group differences using:
Pearson’s chi-squared test, when all expected cell counts were adequate
Fisher’s exact test, when sparse cells were present.
the corresponding hypotheses were:
$$
\begin{aligned}
H_{0}:~& P(X = k \mid \text{Placebo}) = P(X = k \mid \text{Dpen}) 
\quad \text{for all } k, \\[6pt]
H_{1}:~& \exists\, k \text{ such that } 
P(X = k \mid \text{Placebo}) \neq P(X = k \mid \text{Dpen})
\end{aligned}
$$


``` r
# Yixin: use the randomized dataset, drop NAs
cat_vars <- c("sex", "ascites", "hepatomegaly", "spiders", "edema", "stage")

baseline_cat <- lapply(cat_vars, function(v) {
  dat_v <- cirr_trial |>
    select(drug, !!sym(v)) |>
    drop_na()  # Yixin: remove missing for drugs only
  
  tab  <- table(dat_v[[v]], dat_v$drug)
  prop <- prop.table(tab, margin = 2)  # column % within each treatment group
  
  # Yixin: choose Pearson chi-squared if all expected counts are reasonable, otherwise choose Fisher
  if (any(tab < 5)) {
    test <- fisher.test(tab)
    test_name <- "Fisher's exact"
  } else {
    test <- chisq.test(tab, correct = FALSE)
    test_name <- "Chi-squared"
  }
  
  df <- as.data.frame(tab)
  df$Percent <- round(100 * as.vector(prop), 1)
  df$Label   <- paste0(df$Freq, " (", df$Percent, "%)")
  
  wide <- tidyr::pivot_wider(
    df,
    id_cols  = Var1,
    names_from = Var2,
    values_from = Label
  )
  
  names(wide)[1] <- "Level"
  wide$Variable <- v
  wide$Test     <- test_name
  wide$Pvalue   <- test$p.value
  wide
}) |>
  bind_rows()

# Yixin: replace Variable codes with full names
baseline_cat <- baseline_cat |>
  mutate(
    Variable = dplyr::recode(
      Variable,
      sex = "Sex",
      ascites = "Ascites",
      hepatomegaly = "Hepatomegaly",
      spiders = "Spider angiomas",
      edema = "Edema",
      stage = "Histologic stage"
    )
  )

# Yixin: reorder columns for readability and give clearer names
baseline_cat <- baseline_cat |>
  select(Variable, Level, Placebo, `D-penicillamine`, Test, Pvalue)

colnames(baseline_cat) <- c(
  "Baseline variable",
  "Level",
  "Placebo n (%)",
  "D-penicillamine n (%)",
  "Test type",
  "p-value"
)

knitr::kable(
  baseline_cat,
  digits = 3,
  caption = "Baseline categorical variables by treatment group (randomized patients)"
)
```



Table: Baseline categorical variables by treatment group (randomized patients)

|Baseline variable |Level |Placebo n (%) |D-penicillamine n (%) |Test type      | p-value|
|:-----------------|:-----|:-------------|:---------------------|:--------------|-------:|
|Sex               |F     |139 (90.3%)   |137 (86.7%)           |Chi-squared    |   0.326|
|Sex               |M     |15 (9.7%)     |21 (13.3%)            |Chi-squared    |   0.326|
|Ascites           |N     |144 (93.5%)   |144 (91.1%)           |Chi-squared    |   0.433|
|Ascites           |Y     |10 (6.5%)     |14 (8.9%)             |Chi-squared    |   0.433|
|Hepatomegaly      |N     |67 (43.5%)    |85 (53.8%)            |Chi-squared    |   0.069|
|Hepatomegaly      |Y     |87 (56.5%)    |73 (46.2%)            |Chi-squared    |   0.069|
|Spider angiomas   |N     |109 (70.8%)   |113 (71.5%)           |Chi-squared    |   0.885|
|Spider angiomas   |Y     |45 (29.2%)    |45 (28.5%)            |Chi-squared    |   0.885|
|Edema             |N     |131 (85.1%)   |132 (83.5%)           |Chi-squared    |   0.877|
|Edema             |S     |13 (8.4%)     |16 (10.1%)            |Chi-squared    |   0.877|
|Edema             |Y     |10 (6.5%)     |10 (6.3%)             |Chi-squared    |   0.877|
|Histologic stage  |1     |4 (2.6%)      |12 (7.6%)             |Fisher's exact |   0.205|
|Histologic stage  |2     |32 (20.8%)    |35 (22.2%)            |Fisher's exact |   0.205|
|Histologic stage  |3     |64 (41.6%)    |56 (35.4%)            |Fisher's exact |   0.205|
|Histologic stage  |4     |54 (35.1%)    |55 (34.8%)            |Fisher's exact |   0.205|
From the table above, we see that, for all categorical variables, the Chi-squared or Fisher’s exact tests produced non-significant p-values. Therefore, we fail to reject the null hypothesis, indicating there's no evidence of imbalance in the distribution of any categorical baseline variables. The two randomized groups is therefore comparable with respect to these characteristics.


