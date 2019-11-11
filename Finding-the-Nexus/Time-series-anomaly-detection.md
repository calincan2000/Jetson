Anomaly detection on vegetation indexes
================

We begin by loading the necessary packages

``` r
library(tidyverse)
library(anomalize)
library(tidyquant)
library(timetk)
library(sweep)
library(imputeTS)
```

Then we load the data directly from
github

``` r
displacement <- read_csv("https://raw.githubusercontent.com/omdena/UNHCR/master/task7_data2heatmap_conversion/displacement_conflict_weekly.csv")

vegetation_indexes <- read_csv("https://raw.githubusercontent.com/omdena/UNHCR/master/task7_data2heatmap_conversion/vegetation_indexes.csv")
```

We observe some negative values in the weeks column in the vegetation
index
    dataset.

``` r
unique(vegetation_indexes$week)
```

    ##  [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
    ## [18]  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34
    ## [35]  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51
    ## [52]  52 -28 -29 -30 -31 -32  -3  -4  -5  -6  -7  -8  -9 -10 -11 -12

So we transform it to have only positive values in that column.

``` r
vegetation_indexes <- vegetation_indexes %>%
  mutate(week = abs(week))
```

We check for data completeness, taking into account we only have 41
weeks for year 2019

``` r
n_regions <- n_distinct(vegetation_indexes$region)
n_weeks   <- n_distinct(vegetation_indexes$week)
n_years   <- n_distinct(vegetation_indexes$year) - 1
  
n_dimensions_complete <- n_regions*n_weeks*n_years + n_regions* 41

print(n_dimensions_complete)
```

    ## [1] 35370

Now we look at the number of lines of our data set.

``` r
nrow(vegetation_indexes)
```

    ## [1] 34443

So we have implicit missing data. We will turn implicit missing values
into explicit missing values.

``` r
vegetation_indexes <- vegetation_indexes %>%
  complete(region, year, week) %>%
  filter(!(year == 2019 & week > 41 ))

nrow(vegetation_indexes)
```

    ## [1] 35370

We have now the number of rows than what we should have. We then find
out the percentage of missing values.

``` r
vegetation_indexes <- vegetation_indexes %>%    
  mutate(is_na = ifelse(is.na(VHI), 1, 0))

round(sum(vegetation_indexes$is_na)/nrow(vegetation_indexes),3)
```

    ## [1] 0.026

We create a date column, so dates are easier to manipulate.

``` r
n_periods <- nrow(vegetation_indexes)/n_regions

vegetation_indexes <- vegetation_indexes %>%
  arrange(region, year, week) %>%
  group_by(region) %>%
  mutate(date = seq(as.Date("1982/1/1"), by = "week", length.out = n_periods),
         date = ymd(date))
```

Impute missing values

``` r
vegetation_indexes <- vegetation_indexes %>%
  group_by(region) %>%
  mutate(SMN = na_ma(SMN),
         SMT = na_ma(SMT),
         VCI = na_ma(VCI),
         TCI = na_ma(TCI),
         VHI = na_ma(VHI))
```

We group our data frame by region and then nest it, so it is easier to
manipulate. This way, all we do will be done on a per-region basis.

``` r
vegetation_indexes_nest <- vegetation_indexes %>%
  group_by(region) %>%
  nest()

vegetation_indexes_nest_ts <- vegetation_indexes_nest %>%
  mutate(data_ts = map(.x = data, 
                       .f = tk_ts,
                       select = c('SMN', 'SMT', 'VCI', 'TCI', 'VHI'),
                       start = 1982,
                       freq = 52))
```

We then perfom a seasonal decomposition, so that we can get the seasonal
and remainder components of our time series.

``` r
vegetation_indexes_anomalies <- vegetation_indexes_nest %>%
  mutate(ts_decompose = map(.x = data,
                            ~time_decompose(data = .x,
                                            target = VHI,
                                            method = "stl",
                                            trend = "5 years")))
```

So our results look like this,

``` r
head(vegetation_indexes_anomalies$ts_decompose[[1]], 10)
```

    ## # A time tibble: 10 x 5
    ## # Index: date
    ##    date       observed  season trend remainder
    ##    <date>        <dbl>   <dbl> <dbl>     <dbl>
    ##  1 1982-01-01     41.4  1.16    64.3    -24.1 
    ##  2 1982-01-08     41.9  0.976   64.1    -23.1 
    ##  3 1982-01-15     44.7  0.506   63.9    -19.7 
    ##  4 1982-01-22     48.1  0.0195  63.8    -15.7 
    ##  5 1982-01-29     49.5 -0.681   63.6    -13.5 
    ##  6 1982-02-05     51.5 -1.26    63.5    -10.7 
    ##  7 1982-02-12     52.9 -1.49    63.3     -8.95
    ##  8 1982-02-19     54.4 -1.33    63.2     -7.44
    ##  9 1982-02-26     55.8 -1.12    63.0     -6.06
    ## 10 1982-03-05     56.0 -0.330   62.9     -6.50

Then we detect anomalies using the remainder. In short, we will
establish lower and upper bounds for the remainder component of our
time-series. If the remainder component is out of the established
bounds, it will be flagged as an anomaly.

``` r
vegetation_indexes_anomalies <- vegetation_indexes_anomalies %>%
  mutate(anomalies = map(.x= ts_decompose,
                            ~anomalize(data = .x,
                                        target = remainder,
                                        method = "iqr",
                                        alpha = 0.2)))
```

So now we have a column that flags anomalies,

``` r
head(vegetation_indexes_anomalies$anomalies[[1]], 10)
```

    ## # A time tibble: 10 x 8
    ## # Index: date
    ##    date       observed  season trend remainder remainder_l1 remainder_l2
    ##    <date>        <dbl>   <dbl> <dbl>     <dbl>        <dbl>        <dbl>
    ##  1 1982-01-01     41.4  1.16    64.3    -24.1         -17.4         17.2
    ##  2 1982-01-08     41.9  0.976   64.1    -23.1         -17.4         17.2
    ##  3 1982-01-15     44.7  0.506   63.9    -19.7         -17.4         17.2
    ##  4 1982-01-22     48.1  0.0195  63.8    -15.7         -17.4         17.2
    ##  5 1982-01-29     49.5 -0.681   63.6    -13.5         -17.4         17.2
    ##  6 1982-02-05     51.5 -1.26    63.5    -10.7         -17.4         17.2
    ##  7 1982-02-12     52.9 -1.49    63.3     -8.95        -17.4         17.2
    ##  8 1982-02-19     54.4 -1.33    63.2     -7.44        -17.4         17.2
    ##  9 1982-02-26     55.8 -1.12    63.0     -6.06        -17.4         17.2
    ## 10 1982-03-05     56.0 -0.330   62.9     -6.50        -17.4         17.2
    ## # … with 1 more variable: anomaly <chr>

Now, we will recompose our time series. If the remainder component of
the observed value is flagged as an anomaly, then the observed value
will be flagged as an anomaly.

``` r
vegetation_indexes_anomalies <- vegetation_indexes_anomalies %>%
  mutate(ts_recompose = map(.x= anomalies,
                            ~time_recompose(data = .x)))
```

So we now have a data frame that looks like this,

``` r
head(vegetation_indexes_anomalies$ts_recompose[[1]], 10)
```

    ## # A time tibble: 10 x 10
    ## # Index: date
    ##    date       observed  season trend remainder remainder_l1 remainder_l2
    ##    <date>        <dbl>   <dbl> <dbl>     <dbl>        <dbl>        <dbl>
    ##  1 1982-01-01     41.4  1.16    64.3    -24.1         -17.4         17.2
    ##  2 1982-01-08     41.9  0.976   64.1    -23.1         -17.4         17.2
    ##  3 1982-01-15     44.7  0.506   63.9    -19.7         -17.4         17.2
    ##  4 1982-01-22     48.1  0.0195  63.8    -15.7         -17.4         17.2
    ##  5 1982-01-29     49.5 -0.681   63.6    -13.5         -17.4         17.2
    ##  6 1982-02-05     51.5 -1.26    63.5    -10.7         -17.4         17.2
    ##  7 1982-02-12     52.9 -1.49    63.3     -8.95        -17.4         17.2
    ##  8 1982-02-19     54.4 -1.33    63.2     -7.44        -17.4         17.2
    ##  9 1982-02-26     55.8 -1.12    63.0     -6.06        -17.4         17.2
    ## 10 1982-03-05     56.0 -0.330   62.9     -6.50        -17.4         17.2
    ## # … with 3 more variables: anomaly <chr>, recomposed_l1 <dbl>,
    ## #   recomposed_l2 <dbl>

We can then observe the results for the Bari region,

``` r
vegetation_indexes_anomalies$ts_recompose[[4]] %>%
  plot_anomalies(time_recomposed = TRUE, alpha_dots = 0.4, alpha_ribbon = 0.1, fill_ribbon = "grey80")
```

![](Time-series-anomaly-detection_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Now, we unnest all our data frames.

``` r
vegetation_indexes_anomalies_unnest <- vegetation_indexes_anomalies %>%
  unnest(ts_recompose) %>%
  select(region,  date, observed, season, trend, remainder, remainder_l1, 
         remainder_l2, anomaly, recomposed_l1, recomposed_l2)
```

Observe the structure

``` r
head(vegetation_indexes_anomalies_unnest)
```

    ## # A time tibble: 6 x 11
    ## # Index:  date
    ## # Groups: region [1]
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

And we count the number of anomalies

``` r
table(vegetation_indexes_anomalies_unnest$anomaly)
```

    ## 
    ##    No   Yes 
    ## 31417  3953
