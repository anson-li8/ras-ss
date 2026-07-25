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
# burden test statistic
t_burden <- function(w, Z, R) {
  denom <- sqrt(as.numeric(t(w) %*% R %*% w)) # quadratic form, sqrt(w'Rw), standard deviation, no need to invert R, so no ridge/regularization needed
  if (!is.finite(denom) || denom <= 0) return(0) # safety check if denominator is somehow broken, just say no signal instead of crash
  as.numeric((t(w) %*% Z) / denom) # actual T_burden, weights dot-multiplied w/ Z-scores in numerator, divided by std dev
}

# summary-stat RAS scan: for each pivotal SNP, find window with minimum p-value
# ras-ss replacement for individual-level LPRS + regression step
scan_ss <- function(b_disc, z_targ, R, prune_filter) {
  n_snps <- length(b_disc) # how many total SNPs we are scanning across
  sites <- seq(1, n_snps, by = skip1) # pivotal SNPS, picked at regular intervals instead of every single SNP, to save time
  sub_windows  <- c(0, seq(min_window_size, max_window_size, by = skip2)) # build candidate window sizes for adaptive window, from 0 up to max
  if (sub_windows[length(sub_windows)] != max_window_size) sub_windows <- c(sub_windows, max_window_size) # make sure biggest candidate window size always included even if step size goes past it
  y_profile <- sapply(sites, function(s) { # for every pivotal SNP s, compute its RAS value, so we can get the time-series data
    best_p <- 1 # start out assuming least significant p-value possible, we will find smallest one
    for (ws in sub_windows) { # loop through every candidate window size for pivotal SNP, adaptive window search
      win_snps <- if (ws == 0) s else max(1, s - ws):min(n_snps, s + ws) # adaptive window, pivotal SNP s plus ws SNPs on either side, clip so we don't exceed the chromosome
      win_snps <- win_snps[prune_filter[win_snps]] # discard any SNPs that were removed by LD pruning earlier
      if (length(win_snps) < 1) next # if pruning wiped out whole window, skip window size and try next candidate
      w_sub <- b_disc[win_snps] # get discovery-cohort weights (w) for SNPs in this window
      z_sub <- z_targ[win_snps] # get target-cohort Z-scores (Z) for SNPs in this window
      R_sub <- R[win_snps, win_snps, drop = FALSE] # get LD correlation submatrix for SNPs in this window
      t_val <- t_burden(w_sub, z_sub, R_sub) # calculate T-burden using our helper function
      best_p <- min(best_p, 2 * pnorm(-abs(t_val))) # turn T_burden into two-tailed p-value and keep smallest one seen so far across window sizes
    }
    -log10(max(best_p, .Machine$double.xmin)) # one every window size has been tried, take -log10 of best p-value found, this is RAS for this pivotal SNP, safe guard to prevent taking log of exactly 0
  })
  list(x = sites, y = y_profile) # return pivotal SNP positions and their RAS values, formatted like time-series data CPD algorithm expects
}

# NOTE: no num_rep averaging here. method A has fixed disc/targ stats, nothing to resplit
one_rep <- function(true_beta, seed) {
  set.seed(seed) # reproducibility
  X_d <- sim_genotypes(n_disc, R_true) # simulate discovery cohort genotypes - Method A discovery dataset
  X_t <- sim_genotypes(n_targ, R_true) # simulate target cohort's genotypes - cohort we are testing
  y_d <- X_d %*% true_beta + rnorm(n_disc, 0, 3) # build discovery cohort phenotype from true causal effect plus random noise
  y_t <- X_t %*% true_beta + rnorm(n_targ, 0, 3) # build target cohort phenotype the same way
  disc <- get_marginal_stats(X_d, y_d) # run GWAS on discovery cohort to give weight vector w, per Method A
  targ <- get_marginal_stats(X_t, y_t) # run GWAS on target cohort to give Z-score vector Z we are testing
  scan <- scan_ss(disc$beta, targ$z, R_emp, prune_filter) # run full ras-ss scan, discovery weights (w), target Z-scores (Z), pruned reference LD matrix (R)
  list(scan = scan, X_t = X_t, y_t = y_t, b_disc = disc$beta) # return RAS time series, plus raw target data for individual-level comparison on exact same cohort
}

# scw=5 (package default) only gave ~17-20% first-pass acceptance on true signal
# because slope re-check window sits inside the flat top of the plateau
# scw=8 fixes this, acceptance goes to ~100%
detect_peaks <- function(scan, window_size, scw = 8) {
  stopifnot(window_size < length(scan$x)) # stop if CPD sliding scan window bigger than whole time series
  tryCatch({
    # first pass: changepoint detection
    cp <- ras_detect( # actual, unmodified original package's changepoint detection, sliding checking window running Davies test at each position
    # second pass: local Davies validation
      x = scan$x, y = scan$y,
      window_size                    = window_size, # width of siding checking window, tuned to 12 from simulation 1
      slope_check_window_size        = scw, # how many points on each side used to check left-rising / right-falling slope condition, set to 8 instead of package default of 5
      slope.p.values.threshold.left  = slope_thresh, # how strict left-slope-rising check has to be before call it real peak
      slope.p.values.threshold.right = slope_thresh # same as above, but for right-slope-falling
    )
    # second pass: Davies validation
    val <- ras_validate( # reconfirming change point with tigther local Davies test and going to nearest local maximum, unmodified package function
      this.result = cp, x = scan$x, y = scan$y,
      this.start = 1, this.skip = skip1, # tells ras_validate how pivotal SNP grid was spaced, so it can convert positions back to real SNP coordinates
      second_window_size = 15, # width of second-pass local check window
      min_signal = 2.5, # RAS has to be at least this big for a changepoint to count as real signal, not just noise
      p.value.threshold = davies_thresh # how strict the second-pass Davies test has to be
    )
    list(val = val, candidates = cp$all.changepoints, error = FALSE) # return validated change points, plus all raw canidates before validation
  }, error = function(e) {
    list(val = list(tau_hats = numeric(0)), candidates = NULL, error = TRUE) # something went wrong
  })
}

# individual-level RAS for comparison (original method)
# original LPRS-based calculation
indiv_scan <- function(scan, X_t, y_t, b_disc) {
  sub_windows <- c(0, seq(min_window_size, max_window_size, by = skip2)) # same candidate window size list as scan_ss
  if (sub_windows[length(sub_windows)] != max_window_size) sub_windows <- c(sub_windows, max_window_size) # same as scan_ss
  sapply(scan$x, function(s) { # for every pivotal SNP position, compute ORIGINAL individual-level RAS instead of summary-stat one
    best_p <- 1 # same as before, start assuming least significant
    for (ws in sub_windows) { # loop thru same adaptive window candidates
      win_snps <- if (ws == 0) s else max(1, s - ws):min(n_snps, s + ws) #same window construction as scan_ss
      win_snps <- win_snps[prune_filter[win_snps]] # same LD pruning mask
      if (length(win_snps) < 1) next # same skip-if-empty
      lprs <- X_t[, win_snps, drop = FALSE] %*% b_disc[win_snps] # actual LPRS formula that the original method used, individual-level genotype dosages used
      fit  <- summary(lm(y_t ~ lprs)) # original method's second regression step, regressing real phenotypes against LPRS scores to get p-value
      p    <- if (nrow(fit$coefficients) > 1) fit$coefficients[2, 4] else 1 # pull p-value for LPRS slope or default to not significant if no fit
      best_p <- min(best_p, p) # same as before, keep smallest p-value across window sizes
    }
    -log10(max(best_p, .Machine$double.xmin)) # same conversion to final RAS value
  })
}

# did the validated peak land on the causal cluster (for power calculation)
in_zone <- function(tau) any(tau >= 130 & tau <= 170)

# the package's individual-level workflow, run only a single 50/50 split though
# no LD pruning. returns the validated tau_hats.
native_run <- function(true_beta, seed) {
  set.seed(seed) # reproducibility
  X_t <- sim_genotypes(n_targ, R_true) # simulate cohort's genotypes, no separate discovery cohort here, bc we are testing ORIGINAL algorithm which does the split
  y_t <- as.numeric(X_t %*% true_beta + rnorm(n_targ, 0, 3)) # build cohort phenotype like before
  n <- nrow(X_t) # total number of ppl in cohort
  train <- sort(sample(n, n %/% 2)) # randomly pick half for training split (50/50)
  hold <- setdiff(seq_len(n), train) # everyone not picked for training becomes testing split
  w   <- compute_gwas_weights(X_t, y_t[train], train,
                                               data.frame(row.names = seq_len(n)), TRUE)[, 1] # calls unmodified package function that runs GWAS on training split to get weights
  pgs <- compute_pgs_matrix(X_t, hold, w) # calls rela package function that builds LPRS matrix for testing split, original package func
  pv  <- screen_forward_max_region(X_t, pgs, # call real package's own scanning function
          data.frame(phenotype2 = y_t[hold]), -1, is_continuous = TRUE,
          covariate_formula = "1", skip1 = skip1, skip2 = skip2,
          min_window_size = min_window_size, max_window_size = max_window_size,
          isSimulation = FALSE, isPlot = FALSE)
  x_grid <- seq(1, ncol(X_t), by = skip1) # rebuild pivotal SNP positions to match what package's scan actually used
  stopifnot(length(x_grid) == length(pv))     # grid must equal package profile
  detect_peaks(list(x = x_grid, y = pv), best_ws, scw = 8)$val$tau_hats # run original individual-level result through exact same CPD detection step used earlier for fair comparison
}

# most pure run based on original RAS function
native_run_via_ras <- function(true_beta, seed) {
  set.seed(seed)
  X_t <- sim_genotypes(n_targ, R_true)
  y_t <- as.numeric(X_t %*% true_beta + rnorm(n_targ, 0, 3))
  covs <- data.frame(dummy_cov = rnorm(n_targ))   # harmless, uncorrelated, avoids empty-formula issue
  result <- ras(
    geno = X_t, phenotype = y_t,
    covariates = covs, covariate_cols = "dummy_cov",
    is_continuous = TRUE,
    num_rep = 5,                       # true default averaging of 5 reps, current limitation of ras-ss
    skip1 = skip1, skip2 = skip2,      # skip2 = 3, not the default 20
    min_window_size = min_window_size, max_window_size = max_window_size,
    cp_window_size = best_ws,          # 12
    cp_slope_check_window = 8,
    cp_slope_left = slope_thresh, cp_slope_right = slope_thresh,
    second_window_size = 15,
    second_p_threshold = davies_thresh,
    min_signal = 2.5,
    run_plots = FALSE,                 # don't generate PDFs for every replicate
    save_dir = tempdir()
  )
  result$detection$tau_hats
}
