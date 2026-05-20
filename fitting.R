# Optimization functions
library(optimx)

pregroup <- function(df) {
  df |>
    group_by(Species, Sex, Age) |>
    summarize(
      Weight = sum(Weight * N) / sum(N),
      Variance = poolvar(Variance, N, "pool"),
      N = sum(N),
      .groups = "drop"
    )
}

llhood <- function(obs, preds, pooled_variance) {
  sum(dnorm(obs, preds, sd = sqrt(pooled_variance), log = TRUE))
}

poolvar <- function(x, N, fct = NULL) {
  fct <- fct %||% factor(seq_along(x))
  pooled <- "pool" %in% fct && length(fct) == 1L
  if (pooled) {
    fct <- factor(rep(1, length(x)))
  }
  stopifnot({
    length(x) == length(N)
    length(x) == length(fct)
  })
  spx <- split(x, f = fct)
  spN <- split(N, f = fct)
  out <- numeric(length(levels(fct)))
  it <- 1L
  for (lvl in levels(fct)) {
    dent <- sum(spN[[lvl]] - 1)
    numt <- sum((spN[[lvl]] - 1) * spx[[lvl]])
    out[it] <- numt / dent
    it <- it + 1L
  }
  if (pooled) {
    unique(out)
  } else {
    out[as.integer(fct)]
  }
}

llobjective <- function(pars, input, obs, sigma, mfun) {
  preds <- do.call(mfun, c(list(input), as.list(pars)))
  ll <- llhood(obs, preds, sigma)
  -1 * ll # bobyqa minimizes values
}

opt_fit <- function(df, fct = NULL, mfun = logistic_mod) {
  if (is.null(fct)) {
    fct <- factor(seq_len(NROW(df)))
  } else {
    fct <- factor(df[[fct]])
  }
  this_fun <- mfun
  this_fun_name <- this_fun[["name"]]
  x <- df[["Age"]]
  y <- df[["Weight"]]
  yvar <- df[["Variance"]]
  nvar <- df[["N"]]
  pooledVar <- poolvar(yvar, nvar, fct)
  this_starts <- this_fun[["starts"]](x, y)
  print(this_starts)
  this_lower <- this_fun[["lower_lim"]]
  this_upper <- this_fun[["upper_lim"]]
  this_fit <- optimx(
    par = this_starts,
    fn = llobjective,
    lower = this_lower,
    upper = this_upper,
    method = "bobyqa",
    input = x, obs = y, sigma = pooledVar, mfun = this_fun[["predict"]]
  ) |> suppressWarnings()
  fit_out <- cbind(
    data.frame(model = this_fun_name, optimizer = rownames(this_fit)),
    as.data.frame(this_fit)
  )
  rownames(fit_out) <- NULL
  return(fit_out)
}

fit_all_groups <- function(df, grp_fct = "Sex", err_fct = NULL, mfun = list(logistic_mod)) {
  stopifnot({
    length(grp_fct) > 0
    length(err_fct) > 0 || is.null(err_fct)
  })
  if (length(grp_fct) > 1) {
    df$grp_fct <- factor(apply(sapply(grp_fct, \(x) df[[x]]), 1, paste0, collapse = " "))
  } else {
    df$grp_fct <- factor(df[[grp_fct]])
  }
  grp_fct <- paste0(grp_fct, collapse = "_")
  if (length(err_fct) > 1) {
    df$err_fct <- apply(sapply(err_fct, \(x) df[[x]]), 1, paste0, collapse = " ")
  } else {
    df$err_fct <- df[[err_fct]]
  }
  grp_split_df <- split(df, df$grp_fct)
  grp_split_df <- lapply(
    grp_split_df,
    function(x) {
      each_fun <- lapply(
        mfun,
        function(f) {
          opt_fit(x, fct = "err_fct", mfun = f)
        }
      )
      do.call(rbind, each_fun)
    }
  )
  out <- do.call(rbind, grp_split_df)
  out <- cbind(setNames(data.frame(rownames(out)), grp_fct), out)
  rownames(out) <- NULL
  return(out)
}
