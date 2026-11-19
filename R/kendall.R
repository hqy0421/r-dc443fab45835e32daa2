#' Rolling Kendall's tau-b on a tsibble by time window
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' `slide_kendall()` computes Kendall's tau-b between two measured variables
#' over a time-based rolling window. For every observation the window is the
#' half-open time interval `(t - .width, t]`, where `t` is the value of the
#' index. Only records observed in *both* variables take part in a window;
#' missings are never imputed, and time points absent from the data never
#' contribute a sample.
#'
#' Windows overlap, and records enter and leave them in time order. Instead of
#' recomputing the correlation for each window, the concordant and discordant
#' pair counts are carried forward: each record is inserted once and removed
#' once, and tie counts are updated alongside. This yields `O(n log^2 n)` time
#' and `O(n log n)` space per key for `n` records, rather than the quadratic
#' cost of recomputing each window.
#'
#' @param .x A tsibble with a `POSIXct` index.
#' @param x,y The two measured variables to correlate, as bare column names or
#'   character strings. They must be numeric.
#' @param .width A single positive number giving the window width in seconds.
#'   The window contains records with index in `(t - .width, t]`.
#'
#' @details
#' The result keeps the key, index and grouping of `.x` and appends four
#' columns:
#' * `n_valid`: the number of records in the window observed in both variables.
#' * `n_concordant`: the number of concordant pairs in the window.
#' * `n_discordant`: the number of discordant pairs in the window.
#' * `tau_b`: Kendall's tau-b, `(C - D) / sqrt((n0 - n1) * (n0 - n2))`, where
#'   `n0 = n_valid * (n_valid - 1) / 2`, and `n1` and `n2` count the pairs tied
#'   in `x` and in `y` respectively. It is `NA` when the denominator is zero,
#'   that is when no pair is untied on both sides.
#'
#' All records sharing an index value receive the same window result; a record
#' whose own variables are missing does not enter any window but still carries
#' the value for its time point. Keys are processed independently, so no state
#' crosses key boundaries.
#'
#' @return A tsibble with the same key, index and grouping as `.x`, plus the
#'   four columns described in Details.
#' @family rolling tsibble
#' @export
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' tick <- tsibble(
#'   time = as.POSIXct("2020-01-01", tz = "UTC") + c(0, 30, 60, 90, 120),
#'   price = c(10, 11, 10.5, 12, 12.5),
#'   volume = c(100, 120, 90, 130, 125),
#'   index = time
#' )
#' tick %>% slide_kendall(price, volume, .width = 90)
#'
#' # Per-key, then average the daily correlation
#' library(lubridate)
#' tick %>%
#'   group_by_key() %>%
#'   slide_kendall(price, volume, .width = 90) %>%
#'   index_by(Date = ~ as_date(.)) %>%
#'   summarise(mean_tau_b = mean(tau_b, na.rm = TRUE))
slide_kendall <- function(.x, x, y, .width) {
  UseMethod("slide_kendall")
}

#' @export
slide_kendall.tbl_ts <- function(.x, x, y, .width) {
  x_var <- tidyselect::vars_pull(names(.x), {{ x }})
  y_var <- tidyselect::vars_pull(names(.x), {{ y }})
  idx_var <- index_var(.x)
  new_cols <- c("n_valid", "n_concordant", "n_discordant", "tau_b")

  if (identical(x_var, y_var)) {
    abort("`x` and `y` must be two different columns.")
  }
  if (any(c(x_var, y_var) %in% c(key_vars(.x), idx_var))) {
    abort("`x` and `y` cannot be the key or the index of the tsibble.")
  }
  if (!is.numeric(.x[[x_var]]) || !is.numeric(.x[[y_var]])) {
    abort("`x` and `y` must both be numeric columns.")
  }
  idx <- .x[[idx_var]]
  if (!inherits(idx, "POSIXct")) {
    abort("`slide_kendall()` requires a `POSIXct` index.")
  }
  if (inherits(.width, "difftime")) {
    .width <- as.numeric(.width, units = "secs")
  }
  if (!is.numeric(.width) || length(.width) != 1L || is.na(.width) || .width <= 0) {
    abort("`.width` must be a single positive number (in seconds).")
  }
  clashing <- intersect(new_cols, names(.x))
  if (length(clashing)) {
    abort(sprintf(
      "Can't overwrite existing column%s %s.",
      if (length(clashing) > 1L) "s" else "",
      paste0("`", clashing, "`", collapse = ", ")
    ))
  }

  n <- vec_size(.x)
  res <- list(
    n_valid = integer(n),
    n_concordant = double(n),
    n_discordant = double(n),
    tau_b = double(n)
  )
  if (n > 0L) {
    idx_num <- as.numeric(idx) # POSIXct -> seconds since epoch
    xv <- as.numeric(.x[[x_var]])
    yv <- as.numeric(.x[[y_var]])
    # key_rows() refers to rows grouped by key. The blocks are usually in time
    # order already, but subsetting can leave the `ordered` flag stale, so
    # check the index itself before handing a block to the sweep.
    for (rows in key_rows(.x)) {
      if (length(rows) < 1L) next
      idx_rows <- idx_num[rows]
      if (is.unsorted(idx_rows)) rows <- rows[order(idx_rows)]
      out <- .Call(C_slide_kendall, idx_num[rows], xv[rows], yv[rows], .width)
      res$n_valid[rows] <- out[[1L]]
      res$n_concordant[rows] <- out[[2L]]
      res$n_discordant[rows] <- out[[3L]]
      res$tau_b[rows] <- out[[4L]]
    }
  }

  tbl <- as_tibble(.x)
  for (nm in new_cols) tbl[[nm]] <- res[[nm]]
  update_meta(tbl, .x, ordered = is_ordered(.x), interval = interval(.x))
}

#' @export
slide_kendall.grouped_ts <- slide_kendall.tbl_ts
