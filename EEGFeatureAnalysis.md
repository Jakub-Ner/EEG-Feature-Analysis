---
title: "EEG Feature Analysis"
output:
  html_document: 
    toc: true
    fig_caption: true
    keep_md: true
    df_print: tibble
  html_notebook: default
  pdf_document: default
  word_document: default
---

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter.*

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.

## Introduction

To prepare EEG analysis we utilized (Utilities for Electroencephalographic (EEG) Analysis • eegUtils (craddm.github.io).

We analyze [eeg datasets](https://github.com/jordan-bird/eeg-feature-generation/tree/master) from ["A Study on Mental State Classification using EEG-based Brain-Machine Interface"](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8710576)

## Analysis preparation

We install, import dependencies, set const values and import datasets:


```r
remotes::install_github("craddm/eegUtils")
```

```
## Using GitHub PAT from the git credential store.
```

```
## Skipping install of 'eegUtils' from a github remote, the SHA1 (aeb0ec00) has not changed since last install.
##   Use `force = TRUE` to force installation
```

```r
library(ggplot2)
```

```
## Warning: package 'ggplot2' was built under R version 4.3.3
```

```r
library(eegUtils)
```

```
## 
## Attaching package: 'eegUtils'
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```r
library(readxl)
```

```
## Warning: package 'readxl' was built under R version 4.3.3
```

```r
library(tibble)
```

```
## Warning: package 'tibble' was built under R version 4.3.3
```

```r
options(digits=15) # show whole numbers

# CONSTS
SAMPLING_RATE <- 250


a_concentrating_1_init <- read.csv("./dataset//subjecta-concentrating-1.csv")
a_concentrating_2_init <- read.csv("./dataset//subjecta-concentrating-2.csv")
a_neutral_1_init <- read.csv("./dataset//subjecta-neutral-1.csv")
a_neutral_2_init <- read.csv("./dataset//subjecta-neutral-2.csv")
a_relaxed_1_init <- read.csv("./dataset//subjecta-relaxed-1.csv")
a_relaxed_2_init <- read.csv("./dataset//subjecta-relaxed-2.csv")

b_concentrating_1_init <- read.csv("./dataset//subjectb-concentrating-1.csv")
b_concentrating_2_init <- read.csv("./dataset//subjectb-concentrating-2.csv")
b_neutral_1_init <- read.csv("./dataset//subjectb-neutral-1.csv")
b_neutral_2_init <- read.csv("./dataset//subjectb-neutral-2.csv")
b_relaxed_1_init <- read.csv("./dataset//subjectb-relaxed-1.csv")
b_relaxed_2_init <- read.csv("./dataset//subjectb-relaxed-2.csv")

c_concentrating_1_init <- read.csv("./dataset//subjectc-concentrating-1.csv")
c_concentrating_2_init <- read.csv("./dataset//subjectc-concentrating-2.csv")
c_neutral_1_init <- read.csv("./dataset//subjectc-neutral-1.csv")
#c_neutral_2_init <- read.csv("./dataset//subjectc-neutral-2.csv")
c_relaxed_1_init <- read.csv("./dataset//subjectc-relaxed-1.csv")
c_relaxed_2_init <- read.csv("./dataset//subjectc-relaxed-2.csv")

d_concentrating_1_init <- read.csv("./dataset//subjectd-concentrating-1.csv")
# d_concentrating_2 has to few rows
d_neutral_1_init <- read.csv("./dataset//subjectd-neutral-1.csv")
d_neutral_2_init <- read.csv("./dataset//subjectd-neutral-2.csv")
d_relaxed_1_init <- read.csv("./dataset//subjectd-relaxed-1.csv")
d_relaxed_2_init <- read.csv("./dataset//subjectd-relaxed-2.csv")
```


```r
d_neutral_2_init
```

```
## # A tibble: 15,204 × 6
##     timestamps   TP9   AF7    AF8  TP10 Right.AUX
##          <dbl> <dbl> <dbl>  <dbl> <dbl>     <dbl>
##  1 1533058130. -1.46  6.84 12.7   -77.1     33.7 
##  2 1533058130.  8.30  6.35 -2.93  -81.1    -18.6 
##  3 1533058130. 21.0   7.81  0.488 -75.2     24.9 
##  4 1533058130. 20.5   6.84 16.1   -64.5     65.4 
##  5 1533058130. 11.7   9.77 10.7   -63.0     76.2 
##  6 1533058130.  6.35  9.28  3.91  -70.8     75.2 
##  7 1533058130. 13.2   8.79  0.977 -76.2     44.9 
##  8 1533058130. 28.8  13.7   2.93  -75.7     61.5 
##  9 1533058131. 33.2  17.1  15.1   -71.8     66.4 
## 10 1533058131. 22.9  20.5  18.1   -69.8      1.95
## # ℹ 15,194 more rows
```

## Pre-processing

Firstly we drop first and last 2 seconds of recordings, because it takes time for a subject to obtain a state (concentration, neutral, relaxation).


```r
ROWS_TO_DROP <-2 * SAMPLING_RATE
drop <- function(dataset) {
  dataset <- dataset[-c(1:ROWS_TO_DROP, (nrow(dataset) - ROWS_TO_DROP + 1):10000), ]
  return(dataset)
}

a_concentrating_1 <- drop(a_concentrating_1_init)
a_concentrating_2 <- drop(a_concentrating_2_init)
a_neutral_1 <- drop(a_neutral_1_init)
a_neutral_2 <- drop(a_neutral_2_init)
a_relaxed_1 <- drop(a_relaxed_1_init)
a_relaxed_2 <- drop(a_relaxed_2_init)

b_concentrating_1 <- drop(b_concentrating_1_init)
b_concentrating_2 <- drop(b_concentrating_2_init)
b_neutral_1 <- drop(b_neutral_1_init)
b_neutral_2 <- drop(b_neutral_2_init)
b_relaxed_1 <- drop(b_relaxed_1_init)
b_relaxed_2 <- b_relaxed_2_init[0:9998,]

c_concentrating_1 <- drop(c_concentrating_1_init)
c_concentrating_2 <- drop(c_concentrating_2_init)
c_neutral_1 <- drop(c_neutral_1_init)
#c_neutral_2 <- drop(c_neutral_2_init)
c_relaxed_1 <- drop(c_relaxed_1_init)
c_relaxed_2 <- drop(c_relaxed_2_init)

d_concentrating_1 <- drop(d_concentrating_1_init)
d_neutral_1 <- drop(d_neutral_1_init)
d_neutral_2 <- drop(d_neutral_2_init)
d_relaxed_1 <- drop(d_relaxed_1_init)
d_relaxed_2 <- drop(d_relaxed_2_init)

d_neutral_2
```

```
## # A tibble: 9,998 × 6
##     timestamps   TP9   AF7   AF8  TP10 Right.AUX
##          <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl>
##  1 1533058132.  39.6 23.4   44.9 -34.7     53.7 
##  2 1533058132.  26.9 23.9   31.7 -32.2     66.4 
##  3 1533058132.  18.6 19.5   15.1 -39.1      4.88
##  4 1533058132.  28.3 16.6   11.2 -42.0     13.7 
##  5 1533058132.  41.0 15.1   20.5 -39.6      5.37
##  6 1533058132.  41.0 14.2   24.4 -36.1     15.6 
##  7 1533058132.  28.3 15.6   14.6 -43.5     21.5 
##  8 1533058132.  15.1 12.7   15.6 -54.7     32.2 
##  9 1533058132.  17.6  8.79  16.1 -53.2     84.0 
## 10 1533058132.  33.7 11.7   17.6 -46.9    108.  
## # ℹ 9,988 more rows
```

In the experiment Right.AUX was used for calibration purposes. More specifically as NZ. That's why we rename the mentioned column


```r
rename_to_NZ <- function(dataset){
  #names(dataset)[names(dataset) == "Right.AUX"] <- "NZ"
  return(dataset)
} 


a_concentrating_1 <- rename_to_NZ(a_concentrating_1)
a_concentrating_2 <- rename_to_NZ(a_concentrating_2)
a_neutral_1 <- rename_to_NZ(a_neutral_1)
a_neutral_2 <- rename_to_NZ(a_neutral_2)
a_relaxed_1 <- rename_to_NZ(a_relaxed_1)
a_relaxed_2 <- rename_to_NZ(a_relaxed_2)

b_concentrating_1 <- rename_to_NZ(b_concentrating_1)
b_concentrating_2 <- rename_to_NZ(b_concentrating_2)
b_neutral_1 <- rename_to_NZ(b_neutral_1)
b_neutral_2 <- rename_to_NZ(b_neutral_2)
b_relaxed_1 <- rename_to_NZ(b_relaxed_1)
#b_relaxed_2 <- rename_to_NZ(b_relaxed_2)

c_concentrating_1 <- rename_to_NZ(c_concentrating_1)
c_concentrating_2 <- rename_to_NZ(c_concentrating_2)
c_neutral_1 <- rename_to_NZ(c_neutral_1)
#c_neutral_2 <- rename_to_NZ(c_neutral_2)
c_relaxed_1 <- rename_to_NZ(c_relaxed_1)
c_relaxed_2 <- rename_to_NZ(c_relaxed_2)

d_concentrating_1 <- rename_to_NZ(d_concentrating_1)
d_neutral_1 <- rename_to_NZ(d_neutral_1)
d_neutral_2 <- rename_to_NZ(d_neutral_2)
d_relaxed_1 <- rename_to_NZ(d_relaxed_1)
d_relaxed_2 <- rename_to_NZ(d_relaxed_2)

a_concentrating_1
```

```
## # A tibble: 9,998 × 6
##     timestamps   TP9   AF7     AF8   TP10 Right.AUX
##          <dbl> <dbl> <dbl>   <dbl>  <dbl>     <dbl>
##  1 1533222562. 11.7   41.0  -45.4   -4.88     -40.0
##  2 1533222562. 18.1   36.6    6.84  -7.81      13.7
##  3 1533222562. 35.2   35.6   89.4   -4.39     -30.8
##  4 1533222562. 40.5   36.1   76.7    3.91     -46.4
##  5 1533222562. 23.4   37.1  -40.5   -1.95      34.7
##  6 1533222562.  1.46  35.2  -42.0   -9.77      67.9
##  7 1533222562.  4.39  29.8  -50.8  -10.3       48.3
##  8 1533222562. 28.3   32.2   22.5   -3.42      60.5
##  9 1533222562. 41.0   37.6   33.2    1.95      94.2
## 10 1533222562. 20.0   35.6 -120.    -5.37      48.8
## # ℹ 9,988 more rows
```

To make analysis easier we wrap our dataframes with eeg_epochs


```r
wrap <- function(dataset) return(eeg_epochs(
  dataset[,-1], 
  srate = SAMPLING_RATE, 
  epochs=tibble(epoch=c(1)), 
  timings=data.frame(time=seq(1, nrow(dataset)), epoch=rep(1, nrow(dataset)))
  )
)

a_concentrating_1_obj <- wrap(a_concentrating_1)
a_concentrating_2_obj <- wrap(a_concentrating_2)
a_neutral_1_obj <- wrap(a_neutral_1)
a_neutral_2_obj <- wrap(a_neutral_2)
a_relaxed_1_obj <- wrap(a_relaxed_1)
a_relaxed_2_obj <- wrap(a_relaxed_2)

b_concentrating_1_obj <- wrap(b_concentrating_1)
b_concentrating_2_obj <- wrap(b_concentrating_2)
b_neutral_1_obj <- wrap(b_neutral_1)
b_neutral_2_obj <- wrap(b_neutral_2)
b_relaxed_1_obj <- wrap(b_relaxed_1)
b_relaxed_2_obj <- wrap(b_relaxed_2)

c_concentrating_1_obj <- wrap(c_concentrating_1)
c_concentrating_2_obj <- wrap(c_concentrating_2)
c_neutral_1_obj <- wrap(c_neutral_1)
#c_neutral_2_obj <- wrap(c_neutral_2)
c_relaxed_1_obj <- wrap(c_relaxed_1)
c_relaxed_2_obj <- wrap(c_relaxed_2)

d_concentrating_1_obj <- wrap(d_concentrating_1)
d_neutral_1_obj <- wrap(d_neutral_1)
d_neutral_2_obj <- wrap(d_neutral_2)
d_relaxed_1_obj <- wrap(d_relaxed_1)
d_relaxed_2_obj <- wrap(d_relaxed_2)

a_concentrating_1_obj
```

```
## Epoched EEG data
## 
## Number of channels	: 5 
## Number of epochs	: 1 
## Epoch limits		: 1 - 9998 seconds
## Electrode names		: TP9 AF7 AF8 TP10 Right.AUX 
## Sampling rate		: 250  Hz
## Reference		:
```

In order to remove artifacts we apply bandpass filter. High-pass filter at 1 Hz removes muscle activity, blinking. Low-pass filter at 100 Hz drops electrical interference and environmental artifacts. Also we applied notch filter to remove electrical interference (50 Hz in UK)


```r
apply_filter <- function(eeg_obj){
  eeg_obj <- eeg_filter(eeg_obj, method = "iir", low_freq = 50, high_freq = 49)
  return(
    eeg_filter(eeg_obj, method = "iir", low_freq = 1, high_freq = 100)
  )
}

a_concentrating_1_obj <- apply_filter(a_concentrating_1_obj)
a_concentrating_2_obj <- apply_filter(a_concentrating_2_obj)
a_neutral_1_obj <- apply_filter(a_neutral_1_obj)
a_neutral_2_obj <- apply_filter(a_neutral_2_obj)
a_relaxed_1_obj <- apply_filter(a_relaxed_1_obj)
a_relaxed_2_obj <- apply_filter(a_relaxed_2_obj)

b_concentrating_1_obj <- apply_filter(b_concentrating_1_obj)
b_concentrating_2_obj <- apply_filter(b_concentrating_2_obj)
b_neutral_1_obj <- apply_filter(b_neutral_1_obj)
b_neutral_2_obj <- apply_filter(b_neutral_2_obj)
b_relaxed_1_obj <- apply_filter(b_relaxed_1_obj)
b_relaxed_2_obj <- apply_filter(b_relaxed_2_obj)

c_concentrating_1_obj <- apply_filter(c_concentrating_1_obj)
c_concentrating_2_obj <- apply_filter(c_concentrating_2_obj)
c_neutral_1_obj <- apply_filter(c_neutral_1_obj)
#c_neutral_2_obj <- apply_filter(c_neutral_2_obj)
c_relaxed_1_obj <- apply_filter(c_relaxed_1_obj)
c_relaxed_2_obj <- apply_filter(c_relaxed_2_obj)

d_concentrating_1_obj <- apply_filter(d_concentrating_1_obj)
d_neutral_1_obj <- apply_filter(d_neutral_1_obj)
d_neutral_2_obj <- apply_filter(d_neutral_2_obj)
d_relaxed_1_obj <- apply_filter(d_relaxed_1_obj)
d_relaxed_2_obj <- apply_filter(d_relaxed_2_obj)
```

## Analysis

In order to analyse Brain-Computer-Interface data, we need to utilize simplifications, due to low spatial resolution and noise sensitivity. We leverage left-brain, right-brain dominance theory (Corballis, 2014; Joseph, 1988).

### Spatial analysis

Below we visualize placement of **TP9 (**temporal-parietal) **AF7 (**auricular-frontal**) AF8 TP10 NZ** electrodes on the subjects heads with averaged signals around the area. Blue areas suggest decreased activity, while red increased. Due to low spatial resolution (non-invasive and small amount of channels), analysis results are ambiguous and it's hard to draw conclusions.


```r
show_topplot <- function(eeg_obj, plot_title) {
  return (topoplot(electrode_locations(eeg_obj, drop=TRUE), chan_marker = "name") + ggtitle(plot_title))
}


show_topplot(a_concentrating_1_obj, "A - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
show_topplot(a_concentrating_2_obj, "A - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

```r
show_topplot(b_concentrating_1_obj, "B - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-7-3.png)<!-- -->

```r
show_topplot(b_concentrating_2_obj, "B - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-7-4.png)<!-- -->

```r
show_topplot(c_concentrating_1_obj, "C - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-7-5.png)<!-- -->

```r
show_topplot(c_concentrating_2_obj, "C - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-7-6.png)<!-- -->

```r
show_topplot(d_concentrating_1_obj, "D - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-7-7.png)<!-- -->


```r
show_topplot(a_neutral_1_obj, "A - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
show_topplot(a_neutral_2_obj, "A - neutral - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```r
show_topplot(b_neutral_1_obj, "B - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

```r
show_topplot(b_neutral_2_obj, "B - neutral - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-8-4.png)<!-- -->

```r
show_topplot(c_neutral_1_obj, "C - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-8-5.png)<!-- -->

```r
show_topplot(d_neutral_1_obj, "D - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-8-6.png)<!-- -->

```r
show_topplot(d_neutral_2_obj, "D - neutral - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-8-7.png)<!-- -->


```r
par(mfrow=c(4,2))

show_topplot(a_relaxed_1_obj, "A - relaxed - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
show_topplot(a_relaxed_2_obj, "A - relaxed - 2") 
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

```r
show_topplot(b_relaxed_1_obj, "B - relaxed - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-9-3.png)<!-- -->

```r
show_topplot(b_relaxed_2_obj, "B - relaxed - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-9-4.png)<!-- -->

```r
show_topplot(c_relaxed_1_obj, "C - relaxed - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-9-5.png)<!-- -->

```r
show_topplot(c_relaxed_2_obj, "C - relaxed - 2") 
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-9-6.png)<!-- -->

```r
show_topplot(d_relaxed_1_obj, "D - relaxed - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-9-7.png)<!-- -->

```r
show_topplot(d_relaxed_2_obj, "D - relaxed - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-9-8.png)<!-- -->

### Frequency Analysis (PSD)

To determine contribution of different brain waves, we calculate Power Spectral Density. To obtain less noisy results we use Welch's method instead of Fast Fourier Transform.


```r
show_psd <- function(eeg_obj, title) {
  return (plot_psd(eeg_obj, freq_range = c(1, 100), ) + ggtitle(title))
}
 
show_psd(a_concentrating_1_obj, "A - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
show_psd(a_concentrating_2_obj, "A - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-10-2.png)<!-- -->

```r
show_psd(b_concentrating_1_obj, "B - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-10-3.png)<!-- -->

```r
show_psd(b_concentrating_2_obj, "B - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-10-4.png)<!-- -->

```r
show_psd(c_concentrating_1_obj, "C - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-10-5.png)<!-- -->

```r
show_psd(c_concentrating_2_obj, "C - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-10-6.png)<!-- -->

```r
show_psd(d_concentrating_1_obj, "D - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-10-7.png)<!-- -->

```r
#show_psd(d_concentrating_2_obj, "D - concentraitng - 2")
```


```r
show_psd(a_neutral_1_obj, "A - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
show_psd(a_neutral_2_obj, "A - neutral - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

```r
show_psd(b_neutral_1_obj, "B - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-11-3.png)<!-- -->

```r
show_psd(b_neutral_2_obj, "B - neutral - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-11-4.png)<!-- -->

```r
show_psd(c_neutral_1_obj, "C - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-11-5.png)<!-- -->

```r
show_psd(d_neutral_1_obj, "D - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-11-6.png)<!-- -->

```r
show_psd(d_neutral_2_obj, "D - neutral - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-11-7.png)<!-- -->


```r
show_psd(a_relaxed_1_obj, "A - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
show_psd(a_relaxed_2_obj, "A - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

```r
show_psd(a_relaxed_1_obj, "B - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-12-3.png)<!-- -->

```r
show_psd(a_relaxed_2_obj, "B - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-12-4.png)<!-- -->

```r
show_psd(a_relaxed_1_obj, "C - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-12-5.png)<!-- -->

```r
show_psd(a_relaxed_2_obj, "C - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-12-6.png)<!-- -->

```r
show_psd(a_relaxed_1_obj, "D - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-12-7.png)<!-- -->

```r
show_psd(a_relaxed_2_obj, "D - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-12-8.png)<!-- -->

We noticed that Right.AUX recordings doesn't help in pattern recognition (values always oscillates around similar power for a given frequencies, despite the mental state). That's why for averaged frequency-time analysis we drop this channel.

Below we clearly see, that during concentration brain waves have higher amplitudes. The lowest brain activity is during relaxation.


```r
show_timecourse <- function(eeg_obj, title) {
  return (plot_timecourse(eeg_obj, electrode = c("AF7", "AF8", "TP10", "TP10")) + ggtitle(title))
}
 
show_timecourse(a_concentrating_1_obj, "A - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
show_timecourse(a_concentrating_2_obj, "A - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-13-2.png)<!-- -->

```r
show_timecourse(b_concentrating_1_obj, "B - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-13-3.png)<!-- -->

```r
show_timecourse(b_concentrating_2_obj, "B - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-13-4.png)<!-- -->

```r
show_timecourse(c_concentrating_1_obj, "C - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-13-5.png)<!-- -->

```r
show_timecourse(c_concentrating_2_obj, "C - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-13-6.png)<!-- -->

```r
show_timecourse(d_concentrating_1_obj, "D - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-13-7.png)<!-- -->


```r
show_timecourse(a_neutral_1_obj, "A - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
show_timecourse(a_neutral_2_obj, "A - neutral - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-14-2.png)<!-- -->

```r
show_timecourse(b_neutral_1_obj, "B - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-14-3.png)<!-- -->

```r
show_timecourse(b_neutral_2_obj, "B - neutral - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-14-4.png)<!-- -->

```r
show_timecourse(c_neutral_1_obj, "C - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-14-5.png)<!-- -->

```r
show_timecourse(d_neutral_1_obj, "D - neutral - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-14-6.png)<!-- -->

```r
show_timecourse(d_neutral_2_obj, "D - neutral - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-14-7.png)<!-- -->


```r
show_timecourse(a_relaxed_1_obj, "A - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
show_timecourse(a_relaxed_2_obj, "A - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-15-2.png)<!-- -->

```r
show_timecourse(a_relaxed_1_obj, "B - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-15-3.png)<!-- -->

```r
show_timecourse(a_relaxed_2_obj, "B - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-15-4.png)<!-- -->

```r
show_timecourse(a_relaxed_1_obj, "C - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-15-5.png)<!-- -->

```r
show_timecourse(a_relaxed_2_obj, "C - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-15-6.png)<!-- -->

```r
show_timecourse(a_relaxed_1_obj, "D - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-15-7.png)<!-- -->

```r
show_timecourse(a_relaxed_2_obj, "D - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-15-8.png)<!-- -->

### Feature extraction

In order to better identify patterns or characteristics in the EEG data that are associated with specific mental states, we apply difference between EEG data during different mental states and a baseline (averaged neutral) state to extract features that represent the differences in brain activity across these states.


```r
a_base = (a_neutral_1 + a_neutral_2) / 2
a_concentrating_1_obj <- wrap(a_concentrating_1 - a_base)
a_concentrating_2_obj <- wrap(a_concentrating_2 - a_base)
a_neutral_1_obj <- wrap(a_neutral_1 - a_base)
a_neutral_2_obj <- wrap(a_neutral_2 - a_base)
a_relaxed_1_obj <- wrap(a_relaxed_1 - a_base)
a_relaxed_2_obj <- wrap(a_relaxed_2 - a_base)

b_base = (b_neutral_1 + b_neutral_2) / 2
b_concentrating_1_obj <- wrap(b_concentrating_1 - b_base)
b_concentrating_2_obj <- wrap(b_concentrating_2 - b_base)
b_neutral_1_obj <- wrap(b_neutral_1 - b_base)
b_neutral_2_obj <- wrap(b_neutral_2 - b_base)
b_relaxed_1_obj <- wrap(b_relaxed_1 - b_base)
b_relaxed_2_obj <- wrap(b_relaxed_2 - b_base)

c_base = c_neutral_1
c_concentrating_1_obj <- wrap(c_concentrating_1 - c_base)
c_concentrating_2_obj <- wrap(c_concentrating_2 - c_base)
c_neutral_1_obj <- wrap(c_neutral_1 - c_base)
#c_neutral_2_obj <- wrap(c_neutral_2)
c_relaxed_1_obj <- wrap(c_relaxed_1 - c_base)
c_relaxed_2_obj <- wrap(c_relaxed_2 - c_base)

d_base = (d_neutral_1 + d_neutral_2) / 2
d_concentrating_1_obj <- wrap(d_concentrating_1 - d_base)
d_neutral_1_obj <- wrap(d_neutral_1 - d_base)
d_neutral_2_obj <- wrap(d_neutral_2 - d_base)
d_relaxed_1_obj <- wrap(d_relaxed_1 - d_base)
d_relaxed_2_obj <- wrap(d_relaxed_2 - d_base)
```

After applying the baseline results are more distinct for different states. For concentration power is distributed evenly, whereas for relaxation, the plot is left-skewed.


```r
show_psd <- function(eeg_obj, title) {
  return (plot_psd(eeg_obj, freq_range = c(1, 100), ) + ggtitle(title))
}
 
show_psd(a_concentrating_1_obj, "A - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
show_psd(a_concentrating_2_obj, "A - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

```r
show_psd(b_concentrating_1_obj, "B - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-17-3.png)<!-- -->

```r
show_psd(b_concentrating_2_obj, "B - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-17-4.png)<!-- -->

```r
show_psd(c_concentrating_1_obj, "C - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-17-5.png)<!-- -->

```r
show_psd(c_concentrating_2_obj, "C - concentraitng - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-17-6.png)<!-- -->

```r
show_psd(d_concentrating_1_obj, "D - concentraitng - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-17-7.png)<!-- -->

```r
#show_psd(d_concentrating_2_obj, "D - concentraitng - 2")
```


```r
show_psd(a_relaxed_1_obj, "A - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
show_psd(a_relaxed_2_obj, "A - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-18-2.png)<!-- -->

```r
show_psd(a_relaxed_1_obj, "B - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-18-3.png)<!-- -->

```r
show_psd(a_relaxed_2_obj, "B - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-18-4.png)<!-- -->

```r
show_psd(a_relaxed_1_obj, "C - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-18-5.png)<!-- -->

```r
show_psd(a_relaxed_2_obj, "C - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-18-6.png)<!-- -->

```r
show_psd(a_relaxed_1_obj, "D - relax - 1")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-18-7.png)<!-- -->

```r
show_psd(a_relaxed_2_obj, "D - relax - 2")
```

![](EEGFeatureAnalysis_files/figure-html/unnamed-chunk-18-8.png)<!-- -->

## Implications

Brain waves analysis enables effective interactions with a subject via BCI. In the article we showcase methods for removing artifacts from a signal (bandpass- and notch-filter). PSD was helpful during feature selection and the baseline is an example how to extract useful information from the signal. Leveraging this knowledge can lead to the creation of BCI systems that enable seamless and intuitive interactions between humans and software or computers, opening up new possibilities for communication, control, and interaction for various applications, including assistive technology, gaming, virtual reality, healthcare and education.
