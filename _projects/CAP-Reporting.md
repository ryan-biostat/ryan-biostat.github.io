---
layout: page
title: CAP Lab Instrument QC Report
description: statistical reporting via RMarkdown for measuring instrument performance
img: /assets/img/lab_report.png
importance: 3
category: Biostatistics
---

#### **About**

The Advanced Genomics Lab operates with two complementary missions: a clinical arm focused on patient diagnostics and a research arm dedicated to developing and improving genomic technologies. As a Biostatistician, my primary role is within the research arm, where I build statistical methods for analyzing genomic data and interpreting sequencing results. Secondarily, I get the opportunity to collaborate with the clinical laboratory whenever statistical expertise is needed to support diagnostic workflows.

Our lab relies on several high-performance instruments for sequencing, and maintaining CAP accreditation requires routine evaluation of these machines using manufacturer-defined quality control (QC) metrics. This project of developing a statistical reporting framework for QC assessment of our Scanner and Fluidics systems originated as an open-ended request from the clinical team. It became an awesome opportunity to apply statistical methodology in a real laboratory setting and strengthen my skills in report generation and scientific communication.



## **Purpose**

Clinical laboratories operate under strict regulatory standards to ensure the accuracy, reliability, and safety of diagnostic testing. One of the highest benchmarks in the United States is [CAP accreditation](https://www.cap.org/laboratory-improvement/accreditation/laboratory-accreditation-program), issued by the College of American Pathologists. CAP-accredited labs must demonstrate that their instruments, workflows, and analytical processes consistently meet rigorous quality standards. A key component of this accreditation is the ongoing evaluation of instrument performance using manufacturer-defined quality control (QC) metrics.

In our lab, instrument systems directly support sequencing, so the clinical arm depends on equivalent performance across machines. Subtle deviations can signal emerging issues that may affect sample integrity or raise CAP compliance concerns. This project addresses that need by establishing a routine, statistically grounded reporting system that analyzes QC data on a recurring schedule. Statistical methods provide the framework to distinguish true changes from natural variability and quantify uncertainty, while structured report generation translates these findings into summaries for clinical leadership and potential auditors.

## **Methods**

Generate quarterly statistical reports for our CAP lab to evaluate QC metrics.

1. #### **Understand the Problem**

   To begin this project, I met with our lab director and clinical staff to develop a clear understanding of the instruments involved and how their quality-control metrics are measured. These discussions helped establish context for my reporting. Our laboratory relies on two Scanners and six Fluidics machines, each contributing critical steps in sample processing.  These conversations formed the foundation for the statistical approach by showing how metrics were currently being monitored and what gaps a more formal reporting system needed to address.

2. #### **Interpret the Data**

   The next step was to ingest the available data and evaluate its structure, completeness, and reliability. I began by importing the raw QC files from each Scanner and Fluidics machine, checking for inconsistencies such as missing fields or unexpected value ranges. During this process I worked directly with the manufacturer to clarify how certain values were calculated and what constituted normal operating variation. At the same time, I collaborated with the clinical team to establish a quarterly plan  for data retrieval. This created a reliable pipeline for ongoing reporting and positioned the lab for sustainable long-term monitoring.

3. #### **Generate a Report**

   I developed an automated reporting framework using [RMarkdown](https://rmarkdown.rstudio.com/) to integrate statistical analysis, visualization, and written interpretation into a single reproducible document. The report executes a series of statistical tests tailored to each QC metric, evaluates model assumptions analyses, and generating conclusions based on p-values. For the Fluidics systems, I additionally implemented a [mixed-effects model](https://www.statsmodels.org/stable/mixed_linear.html) to identify which factors (like machine, run date, or sample batch) exert the strongest influence on QC outcomes. This modeling approach allowed the report to move beyond simple comparisons and toward an understanding of underlying drivers in instrument behavior. The resulting RMarkdown report is robust and produces repeatable summaries that can be regenerated each quarter as new data is generated.

4. #### **Communicate Results**

   After completing the analyses, I met with the lab director and clinical staff to walk through the key findings. I presented the workflow from data acquisition through modeling, explained my reasoning behind each analytical step, and highlighted the most relevant outcomes for both each system. My goal was to ensure that the statistical results were understandable, actionable, and aligned with the lab’s operational needs and CAP expectations. These reports continue to be generated and aid in machine evaluation for our lab.

## **Code Sample**

{% highlight R linenos %}

title: "Scanner Comparison - May 2025 Scanner Comparison"
author: "Ryan Gallagher"
date: "`r Sys.Date()`"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(readxl)
library(kableExtra)
library(GGally)
library(htmltools)
library(gridExtra)


T1.file = "REDACTED"
T2.file = "REDACTED"

s3 = read.table(T1.file, header=TRUE, sep='\t') %>% mutate(Scanner = "Scanner 1", Batch="632")
s4 = read.table(T2.file, header=TRUE, sep='\t') %>% mutate(Scanner = "Scanner 2", Batch="632")

names(s3)[6] <- "SNPQC.GreaterThanOrEqual.15.00"
names(s4)[6] <- "SNPQC.GreaterThanOrEqual.15.00"

s = rbind(s3, s4)

t1 = s %>% filter(Scanner == "Scanner 1")
t2 = s %>% filter(Scanner == "Scanner 2")
```

### Analysis Description

The objective of this report is to determine whether there is a measurable difference between the scanners Scanner 1 & Scanner 2. Three metrics are provided: `MAPD`, `SNPQC`, `Waviness.SD` as well as identifiers for scanner. This report will summarize these metrics per scanner then runs statistical tests to determine if the groups are significantly different.

The data used in this report are:

-   `REDACTED` for May 15th, 2025

This batch consisted of 11 patient samples and 1 control (`CytoRef103`).

### Summary Statistics

The objective of this section is to describe the data and determine the distribution within each group to then allow us to compare the data across scanners.

For each metric within each group, we summarize the data and test whether they approximately follow a normal distribution using the Shapiro-Wilks Test. Our interpretation is that if the p-value is **less than 0.05 then we reject the null** **that our data is approximately normal**.

```{r all, fig.align = 'center', echo=FALSE}

library(dplyr)
library(knitr)
library(kableExtra)

summary_table.s = s %>%
  group_by(Scanner) %>%
  summarize(
    N = n(),
    MAPD_mean = mean(MAPD, na.rm = TRUE),
    MAPD_sd = sd(MAPD, na.rm = TRUE),
    MAPD_shapiro_p = if (N > 2) shapiro.test(MAPD)$p.value else NA,
    SNPQC_mean = mean(SNPQC, na.rm = TRUE),
    SNPQC_sd = sd(SNPQC, na.rm = TRUE),
    SNPQC_shapiro_p = if (N > 2) shapiro.test(SNPQC)$p.value else NA,
    Waviness_SD_mean = mean(`Waviness.SD`, na.rm = TRUE),
    Waviness_SD_sd = sd(`Waviness.SD`, na.rm = TRUE),
    Waviness_SD_shapiro_p = if (N > 2) shapiro.test(`Waviness.SD`)$p.value else NA
  ) %>%
  ungroup()

summary_table.s %>%
  kable(
    digits = 2,      
    col.names = c("Scanner", "N", 
                  "Mean", "SD", "Shapiro p-value", 
                  "Mean", "SD", "Shapiro p-value", 
                  "Mean", "SD", "Shapiro p-value"),
    caption = "Summary of Metrics by Scanner with Shapiro-Wilk Test Results"
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"), position = "center") %>%
  add_header_above(c(" " = 2, "MAPD" = 3, "SNPQC" = 3, "Waviness SD" = 3), bold = TRUE, align = "c") %>%
  column_spec(2, border_right = TRUE) %>% 
  column_spec(5, border_right = TRUE) %>%
  column_spec(8, border_right = TRUE) %>%
  column_spec(11, border_right = TRUE)

```

We find that all metrics report a Shapiro-Wilk p-value greater than 0.05. Thus, **we accept the null** that our data is approximately normal. This will aid in accepting assumptions necessary to compare our scanner groups - Scanner 1 vs. Scanner 2.

We'll run a Bartlett test for homogeneity of variances to determine whether the variance is roughly equal for the same metric across groups. Our interpretation is that if the Bartlett p-value is **less than 0.05 then we reject the null that the variances in each of the groups are the same**.

```{r bartlett, fig.align='center', echo=FALSE}
variables = c("MAPD", "SNPQC", "Waviness.SD")
bartlett_results = lapply(variables, function(var) {
  test = bartlett.test(s[[var]] ~ s$Scanner)
  data.frame(
    Variable = var,
    Statistic = round(test$statistic, 2),
    p_value = round(test$p.value, 4)
  )
})

bartlett_table = bind_rows(bartlett_results) %>% 
  rownames_to_column() %>%
  select(-rowname)

bartlett_table %>%
  kable(
    digits = 2,
    align = "c",
    col.names = c("Variable", "Bartlett Statistic", "p-value"),
    caption = "Bartlett's Test for Equality of Variances"
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

We find that our p-values are greater than 0.05, and **we can accept the null** hypothesis that the variances for each metric are the same across groups. This will allow us to accept assumptions of the statistical test in the next section.

### Two-Sample T-Test for Difference of Groups

We have satisfied the following assumptions for each metric in our data:

-   Data is normally distributed within each group

-   Observations are independent

-   Variance is homogeneous between groups

Thus, the Two-Sample T-Test will be employed to determine whether the the means of the two groups are equal. Our interpretation is that if the Two-Sample T-Test p-value is **less than 0.05 then we reject the null that there is no significant difference between the means of the two groups**.

```{r things, fig.align='center', echo=FALSE}

variables = c("MAPD", "SNPQC", "Waviness.SD")

t_test_results = lapply(variables, function(var) {
  test = t.test(t1[[var]], t2[[var]], alternative = "two.sided", paired = FALSE, var.equal = TRUE)
  data.frame(
    Variable = var,
    T_Statistic = round(test$statistic, 2),
    P_Value = round(test$p.value, 4)
  )
})

t_test_table = bind_rows(t_test_results) %>% 
  rownames_to_column() %>%
  select(-rowname)

t_test_table %>%
  kable(
    digits = 2,
    align = "c",
    col.names = c("Variable", "T Statistic", "P-Value"),
    caption = "Two-Sample T-Test Results"
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))

```

We find that our p-values are greater than 0.05, and **we can accept the null hypothesis that there is no significant difference between the means of the two groups**.

### Conclusion

Based on the metrics `MAPD`, `SNPQC`, and `Waviness SD`, **we are confidence that Scanner 1 and Scanner 2 demonstrate comparable performance.**

```{my notes, include=FALSE}
'''
### Discussion of `SNPQC` Distribution

It should be stated first that the Shapiro-Wilks test for normality is a rather powerful test for sample sizes such as ours and that we should be confident in the p-value / our conclusion for normality. Further, we see a distinct agreement in `SNPQC` distribution across both scanners, so the determination of congruence between scanners remains. However, the normality test is reporting a p-value close to 0.05 for both Scanners, and this section will investigate. The distribution for this variable is:

'''
```

{% endhighlight %}