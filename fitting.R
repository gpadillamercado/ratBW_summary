# Optimization functions
library(optimx)
# 1. Split the data that will have different model parameters
# 2. Make a factor variable for any sigma_splits (same model parameters, different sigmas)
# 2a. Save the levels of the sigma_splits and make a list of starts, lower, and upper.
# 3. Fit each split by main data_group, but for each of these:
# 3a. Predict a set of parameters for the whole data group.
# 3b. Split by sigma_split and set the respective sigma.
# 3c. Unsplit the data and estimate the log-likelihood.


fit_data <- function(dat, data_fct, sigma_fct, mod_fun) {
  # Split the data for each data_fct
  stopifnot(!is.null(dat[[data_fct]]))
  stopifnot(!is.null(dat[[sigma_fct]]))
  stopifnot(all(c("Weight", "Age") %in% names(dat)))
  big_split <- factor(dat[[data_fct]])
  big_lvls <- levels(big_split)
  dat_split <- split(dat, big_split)
  grp_fit <- list()
  it <- 1L
  for (this in big_lvls) {
    this_group <- dat_split[[this]]
    # Make a sigma_fct group column for this
    this_group[["sigma_grp"]] <- factor(paste0("sigma_", this_group[[sigma_fct]]))
    grp_sigmas <- levels(this_group$sigma_grp)
    this_sigmas <- list(
      starts = as.list(setNames(rep(1E-2, length(grp_sigmas)), grp_sigmas)),
      lower = as.list(setNames(rep(1E-5, length(grp_sigmas)), grp_sigmas)),
      upper = as.list(setNames(rep(1, length(grp_sigmas)), grp_sigmas))
    )
    grp_fit[[it]] <- fit_group(this_group, this_sigmas, mod_fun)
    it <- it + 1L
  }
  grp_fit <- setNames(grp_fit, big_lvls)
  return(grp_fit)
}

# Define negative log-likelihood function here
llobjective <- function(par, dat, mod_fun) {
  # Which parameters are sigmas? With values.
  sigma_pars <- par[grepl("sigma", names(par))]
  model_pars <- par[!grepl("sigma", names(par))]
  preds <- do.call(mod_fun, c(list(Age = dat$Age), as.list(model_pars)))
  obs <- dat$Weight
  res <- obs - preds

  # Set sigma based on the sigma_grp column
  sigma <- unname(unlist(sigma_pars[dat$sigma_grp]))
  nll <- -1 * sum(dnorm(res, 0, sd = sqrt(sigma), log = TRUE))
  return(nll)
}

fit_group <- function(dat, fct_lst = NULL, mfun = logistic_mod) {
  stopifnot(c("sigma_grp", "Weight", "Age") %in% names(dat))
  # Combine parameter starts & bounds with sigma starts & bounds
  this_start <- c(
    mfun$starts(dat$Age, dat$Weight),
    fct_lst$starts
  )
  this_lower <- c(
    mfun$lower_lim,
    fct_lst$lower
  )
  this_upper <- c(
    mfun$upper_lim,
    fct_lst$upper
  )

  fit <- try(
    optimx::opm(
      unlist(this_start),
      llobjective,
      lower = unlist(this_lower),
      upper = unlist(this_upper),
      method = "bobyqa",
      hessian = FALSE,
      dat = dat,
      mod_fun = mfun$predict
    ) |> as.data.frame()
  )
  if (inherits(fit, "try-error")) {
    return(NULL)
  }
  fit$method <- "bobyqa"
  rownames(fit) <- NULL
  return(fit)
}

get_coefs <- function(fit) {
  as.list(fit[c("A", "k", "tmid")])
}
