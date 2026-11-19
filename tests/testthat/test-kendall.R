library(dplyr, warn.conflicts = FALSE)
library(lubridate)

# Reference implementation: recompute every window from scratch.
kendall_ref <- function(time, x, y, width) {
  n <- length(time)
  n_valid <- integer(n)
  n_conc <- double(n)
  n_disc <- double(n)
  tau <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    w <- which(time > time[i] - width & time <= time[i] & !is.na(x) & !is.na(y))
    n_valid[i] <- length(w)
    if (n_valid[i] >= 2L) {
      cmb <- utils::combn(w, 2L)
      sx <- sign(x[cmb[1L, ]] - x[cmb[2L, ]])
      sy <- sign(y[cmb[1L, ]] - y[cmb[2L, ]])
      n_conc[i] <- sum(sx * sy > 0)
      n_disc[i] <- sum(sx * sy < 0)
      n0 <- n_valid[i] * (n_valid[i] - 1L) / 2
      n1 <- sum(choose(table(x[w]), 2L))
      n2 <- sum(choose(table(y[w]), 2L))
      den <- (n0 - n1) * (n0 - n2)
      if (den > 0) tau[i] <- (n_conc[i] - n_disc[i]) / sqrt(den)
    }
  }
  list(n_valid = n_valid, n_concordant = n_conc, n_discordant = n_disc, tau_b = tau)
}

make_tbl <- function(time, x, y) {
  tsibble(time = time, x = x, y = y, index = time)
}

test_that("slide_kendall() matches a from-scratch reference on random data", {
  set.seed(2024)
  for (trial in 1:40) {
    n <- sample(1:60, 1L)
    tt <- as.POSIXct("2020-01-01", tz = "UTC") + cumsum(sample(1:4, n, TRUE))
    x <- sample(c(NA, 1:7), n, TRUE)
    y <- sample(c(NA, 1:7), n, TRUE)
    width <- sample(c(1, 2, 5, 20, 1e6), 1L)
    got <- slide_kendall(make_tbl(tt, x, y), x, y, .width = width)
    want <- kendall_ref(as.numeric(tt), x, y, width)
    expect_equal(got$n_valid, want$n_valid)
    expect_equal(got$n_concordant, want$n_concordant)
    expect_equal(got$n_discordant, want$n_discordant)
    expect_equal(got$tau_b, want$tau_b)
  }
})

test_that("slide_kendall() agrees with stats::cor(method = 'kendall')", {
  set.seed(7)
  x <- round(rnorm(40), 1)
  y <- round(rnorm(40), 1)
  tt <- as.POSIXct("2020-01-01", tz = "UTC") + 1:40
  res <- slide_kendall(make_tbl(tt, x, y), x, y, .width = 100)
  expect_equal(res$tau_b[40], cor(x, y, method = "kendall"))
})

test_that("the window is half-open, (t - width, t]", {
  tt <- as.POSIXct("2020-01-01", tz = "UTC") + c(0, 10, 20)
  tbl <- make_tbl(tt, c(1, 2, 3), c(1, 2, 3))
  expect_equal(slide_kendall(tbl, x, y, .width = 10)$n_valid, c(1L, 1L, 1L))
  expect_equal(slide_kendall(tbl, x, y, .width = 10.001)$n_valid, c(1L, 2L, 2L))
})

test_that("ties are handled exactly and a zero denominator gives NA", {
  tt <- as.POSIXct("2020-01-01", tz = "UTC") + 0:3
  all_tied <- make_tbl(tt, c(1, 1, 1, 1), c(2, 2, 2, 2))
  res <- slide_kendall(all_tied, x, y, .width = 10)
  # every record is in every window up to its own time point
  expect_equal(res$n_valid, 1:4)
  expect_equal(res$n_concordant, rep(0, 4))
  expect_equal(res$n_discordant, rep(0, 4))
  expect_true(all(is.na(res$tau_b)))

  # tied on one side only also zeroes the denominator
  x_tied <- make_tbl(tt, c(1, 1, 1, 1), c(1, 2, 3, 4))
  expect_true(all(is.na(slide_kendall(x_tied, x, y, .width = 10)$tau_b)))

  # a single untied pair is enough
  one_pair <- make_tbl(tt, c(1, 2, 2, 2), c(1, 2, 2, 2))
  expect_equal(slide_kendall(one_pair, x, y, .width = 10)$tau_b[2], 1)
})

test_that("missing values are neither imputed nor counted", {
  tt <- rep(as.POSIXct("2020-01-01", tz = "UTC") + 0:2, 2)
  tbl <- tsibble(
    time = tt, g = rep(c("a", "b"), each = 3),
    x = c(1, 2, 3, 4, 5, NA), y = c(2, 1, 4, 1, 6, 3),
    index = time, key = g
  )
  res <- slide_kendall(tbl, x, y, .width = 10)
  # the NA row is kept and takes the value of its own time point
  expect_equal(nrow(res), nrow(tbl))
  expect_false(is.na(res$tau_b[6]))
  expect_equal(res$tau_b[6], res$tau_b[5])
  # ... but it never enters a window
  expect_equal(res$n_valid, c(1L, 2L, 3L, 1L, 2L, 2L))
})

test_that("time gaps do not contribute samples", {
  tt <- as.POSIXct("2020-01-01", tz = "UTC") + c(0, 100, 200)
  tbl <- make_tbl(tt, c(1, 2, 3), c(1, 2, 3))
  expect_equal(slide_kendall(tbl, x, y, .width = 150)$n_valid, c(1L, 2L, 2L))
})

test_that("keys do not share state", {
  set.seed(99)
  grid <- expand.grid(sec = 1:40, k = c("A", "B", "C"), stringsAsFactors = FALSE)
  grid <- grid[order(grid$k, grid$sec), ]
  tt <- as.POSIXct("2020-01-01", tz = "UTC") + grid$sec
  x <- round(rnorm(nrow(grid)), 1)
  y <- round(rnorm(nrow(grid)), 1)
  tbl <- tsibble(time = tt, k = grid$k, x = x, y = y, index = time, key = k)
  res <- slide_kendall(tbl, x, y, .width = 6)

  for (kk in c("A", "B", "C")) {
    sub <- tbl[tbl$k == kk, ]
    want <- kendall_ref(as.numeric(sub$time), sub$x, sub$y, 6)
    got <- res[res$k == kk, ]
    expect_equal(got$n_valid, want$n_valid)
    expect_equal(got$tau_b, want$tau_b)
  }
})

test_that("key, index and grouping are preserved", {
  set.seed(5)
  grid <- expand.grid(sec = 1:10, k = c("a", "b"), g = c("x", "y"))
  grid <- grid[order(grid$k, grid$sec), ]
  tbl <- tsibble(
    time = as.POSIXct("2020-01-01", tz = "UTC") + grid$sec,
    k = grid$k, k2 = grid$g, x = runif(nrow(grid)), y = runif(nrow(grid)),
    index = time, key = c(k, k2)
  )
  res <- slide_kendall(tbl, x, y, .width = 5)
  expect_s3_class(res, "tbl_ts")
  expect_identical(key_vars(res), c("k", "k2"))
  expect_identical(index_var(res), "time")
  expect_identical(nrow(res), nrow(tbl))
  expect_identical(names(res), c(names(tbl), "n_valid", "n_concordant", "n_discordant", "tau_b"))

  grped <- group_by(tbl, k2)
  res_g <- slide_kendall(grped, x, y, .width = 5)
  expect_s3_class(res_g, "grouped_ts")
  expect_identical(group_vars(res_g), "k2")
  expect_equal(as_tibble(res_g)$tau_b, res$tau_b)
})

test_that("index_by() + summarise() gives a daily correlation report", {
  set.seed(1)
  grid <- expand.grid(sec = 0:86399, k = c("a", "b"), stringsAsFactors = FALSE)
  grid <- grid[grid$sec %% 3600 == 0, ]
  tbl <- tsibble(
    time = as.POSIXct("2020-01-01", tz = "UTC") + grid$sec + 86400 * rep(0:1, length.out = nrow(grid)),
    k = grid$k, x = round(rnorm(nrow(grid)), 2), y = round(rnorm(nrow(grid)), 2),
    index = time, key = k
  )
  res <- slide_kendall(tbl, x, y, .width = 3600 * 6)
  daily <- res %>%
    group_by_key() %>%
    index_by(Date = ~ as_date(.)) %>%
    summarise(
      mean_tau_b = mean(tau_b, na.rm = TRUE),
      n_valid_windows = sum(!is.na(tau_b))
    )
  expect_s3_class(daily, "tbl_ts")
  expect_identical(index_var(daily), "Date")
  expect_identical(key_vars(daily), "k")
  expect_equal(nrow(daily), 4L) # two keys x two days
  expect_equal(sort(unique(daily$k)), c("a", "b"))
  expect_true(all(daily$n_valid_windows > 0))
  expect_true(all(daily$mean_tau_b >= -1 & daily$mean_tau_b <= 1))
})

test_that("slide_kendall() handles empty and single-row inputs", {
  empty <- tsibble(
    time = as.POSIXct(character(0), tz = "UTC"),
    x = numeric(0), y = numeric(0), index = time
  )
  res <- slide_kendall(empty, x, y, .width = 10)
  expect_identical(nrow(res), 0L)
  expect_true(all(c("n_valid", "tau_b") %in% names(res)))

  one <- make_tbl(as.POSIXct("2020-01-01", tz = "UTC"), 1, 2)
  res1 <- slide_kendall(one, x, y, .width = 10)
  expect_equal(res1$n_valid, 1L)
  expect_true(is.na(res1$tau_b))
})

test_that("slide_kendall() accepts strings, difftime and unsorted rows", {
  tt <- as.POSIXct("2020-01-01", tz = "UTC") + c(5, 1, 3)
  tbl <- make_tbl(tt, c(3, 1, 2), c(3, 1, 2))
  by_name <- slide_kendall(tbl, "x", "y", .width = 10)
  by_bare <- slide_kendall(tbl, x, y, .width = as.difftime(10, units = "secs"))
  expect_equal(by_name$tau_b, by_bare$tau_b)

  # same data in a different row order gives the same answer at each time point
  shuffled <- tbl[c(2, 1, 3), ]
  res_shuffled <- slide_kendall(shuffled, x, y, .width = 10)
  res_sorted <- slide_kendall(tbl, x, y, .width = 10)
  ord <- order(res_shuffled$time)
  expect_equal(res_shuffled$tau_b[ord], res_sorted$tau_b[order(res_sorted$time)])
  expect_equal(res_shuffled$n_valid[ord], res_sorted$n_valid[order(res_sorted$time)])
})

test_that("slide_kendall() validates its inputs", {
  tbl <- make_tbl(as.POSIXct("2020-01-01", tz = "UTC") + 0:2, 1:3, 4:6)
  expect_error(slide_kendall(tbl, x, y, .width = -1), "positive number")
  expect_error(slide_kendall(tbl, x, y, .width = c(1, 2)), "positive number")
  expect_error(slide_kendall(tbl, x, y, .width = NA_real_), "positive number")
  expect_error(slide_kendall(tbl, x, x, .width = 5), "two different")
  expect_error(slide_kendall(tbl, time, y, .width = 5), "key or the index")

  chr_tbl <- tsibble(
    time = as.POSIXct("2020-01-01", tz = "UTC") + 0:2,
    x = letters[1:3], y = 1:3, index = time
  )
  expect_error(slide_kendall(chr_tbl, x, y, .width = 5), "numeric")

  date_tbl <- tsibble(date = as.Date("2020-01-01") + 0:2, x = 1:3, y = 4:6, index = date)
  expect_error(slide_kendall(date_tbl, x, y, .width = 5), "POSIXct")

  clash <- make_tbl(as.POSIXct("2020-01-01", tz = "UTC") + 0:2, 1:3, 4:6)
  clash$tau_b <- 1:3
  expect_error(slide_kendall(clash, x, y, .width = 5), "Can't overwrite")
})

test_that("slide_kendall() scales beyond quadratic recomputation", {
  skip_on_cran()
  set.seed(42)
  n <- 20000L
  tt <- as.POSIXct("2020-01-01", tz = "UTC") + sort(sample(1:(n * 4), n))
  x <- round(rnorm(n), 2)
  y <- round(x + rnorm(n), 2)
  tbl <- make_tbl(tt, x, y)
  elapsed <- system.time(res <- slide_kendall(tbl, x, y, .width = 600))[["elapsed"]]
  expect_lt(elapsed, 20)

  # spot-check a handful of windows against the reference
  want <- kendall_ref(as.numeric(tt), x, y, 600)
  idx <- sample(n, 5L)
  expect_equal(res$tau_b[idx], want$tau_b[idx])
  expect_equal(res$n_valid[idx], want$n_valid[idx])
  expect_equal(res$n_concordant[idx], want$n_concordant[idx])
})
