# code/ras_ss.R
#
# These read the run config (n_snps, n_disc, n_targ, skip1, skip2,
# min_window_size, max_window_size, prune_filter, R_true, R_emp, slope_thresh,
# davies_thresh, best_ws) from the global environment AT CALL TIME, so define
# those in the analysis Rmd before calling any of these.

# generate n fake patients whose genotypes are correlated based on the LD data
sim_genotypes <- function(n, R) {
  X <- rmvnorm(n, sigma = R) # draw n ppl from multivariate normal w/ covariance R, applys LD correlations to generate accurate simulated SNPs
  X <- round(X) # DNA isn't continuous, genotypic dosage have to be whole number, so nearest integer
  pmin(pmax(X, 0), 2) # 0, 1, 2 mutated alleles, clip so dosage doesn't go out of range
}

# LD pruning quality control w/ r^2 < 0.2
# greedy forward pruning, drop SNPs in high LD with a kept SNP
prune_ld <- function(R, thresh) {
  keep <- rep(TRUE, nrow(R)) # start by assuming every SNP survives pruning
  for (i in seq_len(nrow(R) - 1)) { # go through each SNP in order, i is the "anchor" SNP we are comparing everything after it against
    if (!keep[i]) next # if anchor SNP already pruned earlier, it can't be used to prune anything else
    for (j in (i + 1):nrow(R)) { # compare anchor SNP i against every SNP j that comes after it
      if (keep[j] && R[i, j]^2 > thresh) keep[j] <- FALSE # if SNP j is still in and its LD correlation squared w/ SNP i is above r^2 < 0.2 cutoff, it's redudant so discard
    }
  }
  keep # return which SNPs survived pruning
}

# runs linear regression for every SNP to measure correlation with phenotype y
# calculate beta-hat-j, marginal effect size, from training split
get_marginal_stats <- function(X, y) {
  # Vectorized OLS for massive speedup
  yc <- y - mean(y) # center phenotype (subtract mean) to prep for regression
  Xc <- sweep(X, 2, colMeans(X), "-") # center every SNP column as well by subtracting SNP's avg genotype dosage
  sxx <- colSums(Xc^2) # sum of squares for each SNP, demoninator of regression slope formula
  sxy <- as.numeric(crossprod(Xc, yc)) # all 300 SNP-vs-phenotype regressions in one matrix multiply instead of loop (covariance in slope formula)
  beta <- sxy / sxx # this is beta-hat-j, marginal effect size for every SNP, computed for all SNPs at once
  syy <- sum(yc^2) # total variance in phenotype, need to figure out leftover error next
  rss <- pmax(syy - beta * sxy, 0) # residual sum of squares, floored at 0 in case of rounding error
  se <- sqrt((rss / (length(y) - 2)) / sxx) # standard error of each beta-hat-j, confidence of effect size estimate
  z <- beta / se # Z-score for each SNP, beta-hat-j divided by its own standard error
  z[!is.finite(z)] <- 0 # if SNP has zero variaiton and this divides by zero, just assign Z-score as 0
  list(beta = beta, z = z) # return both raw effect sizes (weights w) and Z-scores (use in T_burden)
}

# T_burden = w'Z / sqrt(w'Rw), from problem statement
t_burden <- function(w, Z, R) {
  denom <- sqrt(as.numeric(t(w) %*% R %*% w))
  if (!is.finite(denom) || denom <= 0) return(0)
  as.numeric((t(w) %*% Z) / denom)
}

# summary-stat RAS scan: for each pivotal SNP, find window with max -log10(p)
scan_ss <- function(b_disc, z_targ, R, mask) {
  n_snps <- length(b_disc)
  sites <- seq(1, n_snps, by = skip1)
  sub_windows  <- c(0, seq(min_window_size, max_window_size, by = skip2))
  if (sub_windows[length(sub_windows)] != max_window_size) sub_windows <- c(sub_windows, max_window_size)
  y_profile <- sapply(sites, function(s) {
    best_p <- 1
    for (ws in sub_windows) {
      win_snps <- if (ws == 0) s else max(1, s - ws):min(n_snps, s + ws)
      win_snps <- win_snps[mask[win_snps]]
      if (length(win_snps) < 1) next
      w_sub <- b_disc[win_snps]
      z_sub <- z_targ[win_snps]
      R_sub <- R[win_snps, win_snps, drop = FALSE]
      num <- sum(w_sub * z_sub)
      denom <- sqrt(as.numeric(w_sub %*% R_sub %*% w_sub))
      if (!is.finite(denom) || denom <= 0) next
      t_val <- num / denom
      best_p <- min(best_p, 2 * pnorm(-abs(t_val)))
    }
    -log10(max(best_p, .Machine$double.xmin))
  })
  list(x = sites, y = y_profile)
}

# NOTE: no num_rep averaging here. method A has fixed disc/targ stats, nothing to resplit
one_rep <- function(true_beta, seed) {
  set.seed(seed)
  X_d <- sim_genotypes(n_disc, R_true)
  X_t <- sim_genotypes(n_targ, R_true)
  y_d <- X_d %*% true_beta + rnorm(n_disc, 0, 3)
  y_t <- X_t %*% true_beta + rnorm(n_targ, 0, 3)
  disc <- get_marginal_stats(X_d, y_d)
  targ <- get_marginal_stats(X_t, y_t)
  scan <- scan_ss(disc$beta, targ$z, R_emp, prune_filter)
  list(scan = scan, X_t = X_t, y_t = y_t, b_disc = disc$beta)
}

# scw=5 (package default) only gave ~17-20% first-pass acceptance on true signal
# because slope re-check window sits inside the flat top of the plateau
# scw=8 fixes this, acceptance goes to ~100%
detect_peaks <- function(scan, window_size, scw = 8) {
  stopifnot(window_size < length(scan$x))
  tryCatch({
    # first pass: changepoint detection
    # use sink() to hide cat() verbose
    null_con <- file(nullfile(), open = "wt")
    sink(null_con, type = "output")
    sink(null_con, type = "message")
    cp <- ras_detect(
      x = scan$x, y = scan$y,
      window_size                    = window_size,
      slope_check_window_size        = scw,
      slope.p.values.threshold.left  = slope_thresh,
      slope.p.values.threshold.right = slope_thresh
    )
    # second pass: local Davies validation
    val <- suppressWarnings(ras_validate(
      this.result = cp, x = scan$x, y = scan$y,
      this.start = 1, this.skip = skip1,
      second_window_size = 15, min_signal = 2.5,
      p.value.threshold = davies_thresh
    ))
    # stop redirecting output
    sink(type = "message")
    sink(type = "output")
    close(null_con)
    list(val = val, candidates = cp$all.changepoints, error = FALSE)
  }, error = function(e) {
    # ensure sink is turned off even if an error occurs
    try({ sink(type = "message"); sink(type = "output"); close(null_con) }, silent = TRUE)
    list(val = list(tau_hats = numeric(0)), candidates = NULL, error = TRUE)
  })
}

# individual-level RAS for comparison (original method)
indiv_scan <- function(scan, X_t, y_t, b_disc) {
  sub_windows <- c(0, seq(min_window_size, max_window_size, by = skip2))
  if (sub_windows[length(sub_windows)] != max_window_size) sub_windows <- c(sub_windows, max_window_size)
  sapply(scan$x, function(s) {
    best_p <- 1
    for (ws in sub_windows) {
      win_snps <- if (ws == 0) s else max(1, s - ws):min(n_snps, s + ws)
      win_snps <- win_snps[prune_filter[win_snps]]
      if (length(win_snps) < 1) next
      lprs <- X_t[, win_snps, drop = FALSE] %*% b_disc[win_snps]
      fit  <- suppressWarnings(summary(lm(y_t ~ lprs)))
      p    <- if (nrow(fit$coefficients) > 1) fit$coefficients[2, 4] else 1
      best_p <- min(best_p, p)
    }
    -log10(max(best_p, .Machine$double.xmin))
  })
}

# did the validated peak land on the causal cluster
in_zone <- function(tau) any(tau >= 130 & tau <= 170)

# the package's individual-level workflow, run as a user would:
# internal 50/50 split, no LD pruning. returns the validated tau_hats.
native_run <- function(true_beta, seed) {
  set.seed(seed)
  X_t <- sim_genotypes(n_targ, R_true)
  y_t <- as.numeric(X_t %*% true_beta + rnorm(n_targ, 0, 3))
  n <- nrow(X_t)
  train <- sort(sample(n, n %/% 2)); hold <- setdiff(seq_len(n), train)
  nc <- file(nullfile(), open = "wt"); sink(nc, type = "output"); sink(nc, type = "message")
  w   <- suppressWarnings(compute_gwas_weights(X_t, y_t[train], train,
                                               data.frame(row.names = seq_len(n)), TRUE))[, 1]
  pgs <- compute_pgs_matrix(X_t, hold, w)
  pv  <- suppressWarnings(screen_forward_max_region(X_t, pgs,
                                                    data.frame(phenotype2 = y_t[hold]), -1, is_continuous = TRUE,
                                                    covariate_formula = "1", skip1 = skip1, skip2 = skip2,
                                                    min_window_size = min_window_size, max_window_size = max_window_size,
                                                    isSimulation = FALSE, isPlot = FALSE))
  sink(type = "message"); sink(type = "output"); close(nc)
  x_grid <- seq(1, ncol(X_t), by = skip1)
  stopifnot(length(x_grid) == length(pv))     # grid must equal package profile
  suppressWarnings(
    detect_peaks(list(x = x_grid, y = pv), best_ws, scw = 8)$val$tau_hats
  )
}
