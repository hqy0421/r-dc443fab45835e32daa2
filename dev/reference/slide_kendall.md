# Rolling Kendall's tau-b on a tsibble by time window

**\[experimental\]**

`slide_kendall()` computes Kendall's tau-b between two measured
variables over a time-based rolling window. For every observation the
window is the half-open time interval `(t - .width, t]`, where `t` is
the value of the index. Only records observed in *both* variables take
part in a window; missings are never imputed, and time points absent
from the data never contribute a sample.

Windows overlap, and records enter and leave them in time order. Instead
of recomputing the correlation for each window, the concordant and
discordant pair counts are carried forward: each record is inserted once
and removed once, and tie counts are updated alongside. This yields
`O(n log^2 n)` time and `O(n log n)` space per key for `n` records,
rather than the quadratic cost of recomputing each window.

## Usage

``` r
slide_kendall(.x, x, y, .width)

# S3 method for class 'tbl_ts'
slide_kendall(.x, x, y, .width)

# S3 method for class 'grouped_ts'
slide_kendall(.x, x, y, .width)
```

## Arguments

- .x:

  A tsibble with a `POSIXct` index.

- x, y:

  The two measured variables to correlate, as bare column names or
  character strings. They must be numeric.

- .width:

  A single positive number giving the window width in seconds. The
  window contains records with index in `(t - .width, t]`.

## Value

A tsibble with the same key, index and grouping as `.x`, plus the four
columns described in Details.

## Details

The result keeps the key, index and grouping of `.x` and appends four
columns:

- `n_valid`: the number of records in the window observed in both
  variables.

- `n_concordant`: the number of concordant pairs in the window.

- `n_discordant`: the number of discordant pairs in the window.

- `tau_b`: Kendall's tau-b, `(C - D) / sqrt((n0 - n1) * (n0 - n2))`,
  where `n0 = n_valid * (n_valid - 1) / 2`, and `n1` and `n2` count the
  pairs tied in `x` and in `y` respectively. It is `NA` when the
  denominator is zero, that is when no pair is untied on both sides.

All records sharing an index value receive the same window result; a
record whose own variables are missing does not enter any window but
still carries the value for its time point. Keys are processed
independently, so no state crosses key boundaries.

## See also

Other rolling tsibble:
[`slide_tsibble()`](https://tsibble.tidyverts.org/dev/reference/slide_tsibble.md),
[`stretch_tsibble()`](https://tsibble.tidyverts.org/dev/reference/stretch_tsibble.md),
[`tile_tsibble()`](https://tsibble.tidyverts.org/dev/reference/tile_tsibble.md)

## Examples

``` r
library(dplyr, warn.conflicts = FALSE)
tick <- tsibble(
  time = as.POSIXct("2020-01-01", tz = "UTC") + c(0, 30, 60, 90, 120),
  price = c(10, 11, 10.5, 12, 12.5),
  volume = c(100, 120, 90, 130, 125),
  index = time
)
tick %>%
  slide_kendall(price, volume, .width = 90)
#> # A tsibble: 5 x 7 [30s] <UTC>
#>   time                price volume n_valid n_concordant n_discordant  tau_b
#>   <dttm>              <dbl>  <dbl>   <int>        <dbl>        <dbl>  <dbl>
#> 1 2020-01-01 00:00:00  10      100       1            0            0 NA    
#> 2 2020-01-01 00:00:30  11      120       2            1            0  1    
#> 3 2020-01-01 00:01:00  10.5     90       3            2            1  0.333
#> 4 2020-01-01 00:01:30  12      130       3            3            0  1    
#> 5 2020-01-01 00:02:00  12.5    125       3            2            1  0.333

# Per-key, then average the daily correlation
library(lubridate)
tick %>%
  group_by_key() %>%
  slide_kendall(price, volume, .width = 90) %>%
  index_by(Date = ~ as_date(.)) %>%
  summarise(mean_tau_b = mean(tau_b, na.rm = TRUE))
#> # A tsibble: 1 x 2 [?]
#>   Date       mean_tau_b
#>   <date>          <dbl>
#> 1 2020-01-01      0.667
```
