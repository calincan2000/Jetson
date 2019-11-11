Lag correlation
================

## Data preparation

We begin by loading the necessary packages

``` r
library(tidyverse)
library(tidyquant)
library(timetk)

library(knitr)
library(kableExtra)
```

Then the loading the anomalies and the displacement csv

``` r
anomalies <- read_csv("anomalies.csv")

displacement <- read_csv("https://raw.githubusercontent.com/omdena/UNHCR/master/task7_data2heatmap_conversion/displacement_conflict_weekly.csv")
```

``` r
head(anomalies)
```

    ## # A tibble: 6 x 11
    ##   region date       observed  season trend remainder remainder_l1
    ##   <chr>  <date>        <dbl>   <dbl> <dbl>     <dbl>        <dbl>
    ## 1 Awdal  1982-01-01     41.4  1.16    64.3     -24.1        -17.4
    ## 2 Awdal  1982-01-08     41.9  0.976   64.1     -23.1        -17.4
    ## 3 Awdal  1982-01-15     44.7  0.506   63.9     -19.7        -17.4
    ## 4 Awdal  1982-01-22     48.1  0.0195  63.8     -15.7        -17.4
    ## 5 Awdal  1982-01-29     49.5 -0.681   63.6     -13.5        -17.4
    ## 6 Awdal  1982-02-05     51.5 -1.26    63.5     -10.7        -17.4
    ## # … with 4 more variables: remainder_l2 <dbl>, anomaly <chr>,
    ## #   recomposed_l1 <dbl>, recomposed_l2 <dbl>

``` r
head(displacement)
```

    ## # A tibble: 6 x 6
    ##   ts         region district   arrivals departures conflict
    ##   <date>     <chr>  <chr>         <dbl>      <dbl>    <dbl>
    ## 1 2016-01-03 Awdal  Baki              0          0        0
    ## 2 2016-01-03 Awdal  Borama            0          0        0
    ## 3 2016-01-03 Awdal  Lughaye           0          0        0
    ## 4 2016-01-03 Awdal  Zeylac            0          0        0
    ## 5 2016-01-03 Bakool Ceel Barde        0          0        0
    ## 6 2016-01-03 Bakool Rab Dhuure        0          0        1

We add two days to the date column of anomalies in order to make it
match with the displacement data set and we add weekly values in order
to have monthly values for the displacement data set. Then we join both
data sets together.

``` r
anomalies <- anomalies %>%
  mutate(date = date + 2)

displacement_monthly <- displacement %>%
  group_by(ts, region)%>%
  summarise(arrivals = sum(arrivals), 
            departures = sum(departures),
            conflict = sum(conflict))

anomalies_displacement <- inner_join(anomalies, displacement_monthly, by = c("region", "date"= "ts"))
```

We nest our data set.

``` r
anomalies_displacement_nest <- anomalies_displacement %>%
  group_by(region) %>%
  tidyr::nest(.)
```

## Lag correlation

### Between observed vegetation index values and departures

And look for lag correlation by district between the observed values in
the vegetation index and departures.

``` r
lag_correlation_observed <- anomalies_displacement_nest %>%
  mutate(correlation = map(.x = data, .y = region, ~ccf(x = .x[,"observed"], y = .x[,"departures"], lag = 25, main = .y)),
         correlation_acf = map(.x = correlation, ~as.tibble(.x$acf))) %>%
  unnest(correlation_acf) %>%
  select(-data, -correlation) %>%
  rename(ccf = V1) %>%
  group_by(region)%>%
  mutate(lag = seq(from = -25, to =, 25)) %>%
  ungroup()
```

![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-5.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-6.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-7.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-8.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-9.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-10.png)<!-- -->

    ## Warning: `as.tibble()` is deprecated, use `as_tibble()` (but mind the new semantics).
    ## This warning is displayed once per session.

![](Lag-correlation_files/figure-gfm/unnamed-chunk-7-11.png)<!-- -->

### Between remainder vegetation index values and departures

And look for lag correlation by district between the remainder values in
the vegetation index and departures.

``` r
lag_correlation_remainder <- anomalies_displacement_nest %>%
  mutate(correlation = map(.x = data, .y = region, ~ccf(x = .x[,"remainder"], y = .x[,"departures"], lag = 25, main = .y)),
         correlation_acf = map(.x = correlation, ~as.tibble(.x$acf))) %>%
  unnest(correlation_acf) %>%
  select(-data, -correlation) %>%
  rename(ccf = V1) %>%
  group_by(region)%>%
  mutate(lag = seq(from = -25, to =, 25)) %>%
  ungroup()
```

![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-7.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-8.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-9.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-10.png)<!-- -->![](Lag-correlation_files/figure-gfm/unnamed-chunk-8-11.png)<!-- -->

## Dynamic linear model

We import dynlm

``` r
library(dynlm)
library(lmtest)
library(sandwich)
library(broom)
```

We fit a simple dynamic linear model per district

``` r
dynlm_fit <- anomalies_displacement_nest %>%
  mutate(dynlm_fit = map(.x = data, ~dynlm(formula = departures ~ observed, data = .x)),
         coef_test = map(.x = dynlm_fit, ~lmtest::coeftest(.x, vcov. = vcovHAC)),
         coef_test_tidy = map(.x = coef_test, ~broom::tidy(.x))) %>%
  unnest(coef_test_tidy) %>%
  select(-data, -dynlm_fit, -coef_test) 

head(dynlm_fit)
```

    ## # A tibble: 6 x 6
    ## # Groups:   region [3]
    ##   region term        estimate std.error statistic p.value
    ##   <chr>  <chr>          <dbl>     <dbl>     <dbl>   <dbl>
    ## 1 Awdal  (Intercept)  -78.9       57.7    -1.37    0.173 
    ## 2 Awdal  observed       5.37       2.17    2.48    0.0141
    ## 3 Bakool (Intercept)  548.       692.      0.792   0.429 
    ## 4 Bakool observed      10.3       15.3     0.669   0.504 
    ## 5 Bari   (Intercept)  228.       137.      1.67    0.0973
    ## 6 Bari   observed       0.118      3.88    0.0304  0.976
