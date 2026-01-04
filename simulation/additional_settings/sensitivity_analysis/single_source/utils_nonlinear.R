# packages

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

tanh_f <- function(x) {
  tanh(x)
}

relu <- function(x) {
  pmax(0, x)
}

if (!require(glmnet)) {
  install.packages("glmnet")
  library(glmnet)
}

if (!require(lpSolve)) {
  install.packages("lpSolve")
  library(lpSolve)
}

# find trim operator Q
Q_trim <- function(X,rho=0.5,trim_form='standard',quant=50){
  UDV_list = svd(X)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)
  if (trim_form == 'standard'){
    tau = quantile(D, rho)
    Dtilde = pmin(D, tau)
    Q = diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*%  t(U)
    return(Q)
  }
  if (trim_form == 'fixed'){
    tau =D[quant]
    Dtilde = pmin(D, tau)
    Q = diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*%  t(U)
    return(Q)
  }
  if (trim_form == 'zero'){
    gap_list = (D[1:(round(length(D)/2)-1)] - D[2:round(length(D)/2)])/D[2:round(length(D)/2)]
    q_hat = which.max(gap_list)
    Q = diag(nrow(X)) - U[,1:q_hat] %*%t(U[,1:q_hat])
    return(Q)
  }
  if (trim_form =='no_trim'){
    Q = diag(nrow(X))
    return(Q)
  }
}

# lasso for y~X
Lasso_func <- function(X,y){
  fit = cv.glmnet(x=X, y=y,intercept = FALSE)
  beta = as.matrix((coef(fit,S =fit$lambda.min)[-1]))
  result = list("beta" = beta,"res" = y-X %*% beta)
  return(result)
}

# find Zt = argmin\|X Z_t-y\|_{\inf}
find_Zt <- function(X,y){
  p <- dim(X)[1]
  n <- dim(X)[2]
  
  
  # Objective: Minimize t
  # Decision variables: [beta_1, beta_2, ..., beta_n,beta_1_, beta_2_, ..., beta_n_, t]
  f.obj <- c(rep(0, 2*n), 1)  # Coefficients for the objective function
  
  # Constraints:
  # 1. X * beta - t <= y
  # 2. -X * beta - t <= -y
  
  # Construct A_ub matrix
  A_upper <- cbind(X,-X,-1)    # X * beta+ - X * beta- - t
  A_lower <- cbind(-X,X,-1)   # -X * beta+ + X * beta- - t
  A_ub <- rbind(A_upper, A_lower)
  
  # Construct b_ub vector
  b_ub <- c(y, -y)
  
  # Solve the linear program
  lp_result <- lp(
    direction = "min",
    objective.in = f.obj,
    const.mat = A_ub,
    const.dir = rep("<=", 2 * p),
    const.rhs = b_ub,
    compute.sens = FALSE,
    dense.const = TRUE
  )
  
  # Check if the optimization was successful
  if (lp_result$status == 0) {
    # Extract the solution
    solution <- lp_result$solution
    beta_opt <- solution[1:n] - solution[(n+1):(2*n)]
    t_opt <- solution[2*n + 1]
    
    # cat("Optimal Zt:\n")
    # print(beta_opt)
    # cat("\nMinimum infinity norm:", t_opt, "\n")
  } else {
    cat("Optimization failed. Status code:", lp_result$status, "\n")
    cat("Message:", lp_result$message, "\n")
  }
  return(list('Z_t'=beta_opt,'loss'=t_opt))
}

# data
generate_dataset = function(Psi,phi,beta,n, p, s, q, sigmaE, sigma, pert){
  #create H and Gamma with N(0,1) values and of appropriate size. H can be tuned with pert
  H = pert*matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE)
  
  #value of X independent from H
  E = matrix(rnorm(n*p,mean=0,sd=sigmaE),n,p,byrow = TRUE)

  #f_H <- sin(H)

  #defined in eq. (2), high-dimensional measured covariates
  X = E + 4/(1+exp(-H%*%Psi)) + 4*sin(H%*%Psi)


  # beta = matrix(rep(c(1,0),times = c(s,p-s)),p,1,byrow = TRUE)

  #nx1 matrix with values of mean 0 and SD of sigma, error in Y independent of X
  nu = matrix(rnorm(n*1,mean=0,sd=sigma),n,1,byrow = TRUE)

  #eq. (1), the response of the Structural Equation Model
  Y = X %*% beta + H %*% phi + nu
  return_list = list("X"= X,"Y"= Y,"Hphi"=H%*%phi)
  return(return_list)
}

generate_dataset_1 = function(Psi, phi, beta, n, p, s, q, sigmaE, sigma, pert, nonlin_type='sigmoid_sin', nonlin_strength=4) {
  H = pert * matrix(rnorm(n*q, mean=0, sd=1), n, q, byrow=TRUE)
  E = matrix(rnorm(n*p, mean=0, sd=sigmaE), n, p, byrow=TRUE)
  
  # Compute non-linear term based on type
  HPsi = H %*% Psi
  if (nonlin_type == 'sigmoid_sin') {
    nonlin_term = nonlin_strength / (1 + exp(-HPsi)) + nonlin_strength * sin(HPsi)
  } else if (nonlin_type == 'tanh') {
    nonlin_term = nonlin_strength * tanh_f(HPsi)
  } else if (nonlin_type == 'relu') {
    nonlin_term = nonlin_strength * relu(HPsi)
  } else if (nonlin_type == 'linear') {
    nonlin_term = nonlin_strength * HPsi
  } else if (nonlin_type == 'sigmoid_only') {
    nonlin_term = nonlin_strength / (1 + exp(-HPsi))
  } else {
    stop("Invalid nonlin_type")
  }
  
  X = E + nonlin_term
  nu = matrix(rnorm(n, mean=0, sd=sigma), n, 1, byrow=TRUE)
  Y = X %*% beta + H %*% phi + nu
  return(list("X"=X, "Y"=Y, "Hphi"=H%*%phi))
}

stable_ols <- function(X, y, ridge = 1e-6) {
  p <- ncol(X)
  n <- nrow(X)
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  # if p>n or XtX ill-conditioned, add tiny ridge
  add_ridge <- (p > n) || (rcond(XtX) < 1e-10)
  if (add_ridge) {
    beta <- solve(XtX + ridge * diag(p), Xty)
  } else {
    # try Cholesky for speed; fall back to solve if needed
    ok <- TRUE
    beta <- tryCatch({
      cholX <- chol(XtX)
      backsolve(cholX, forwardsolve(t(cholX), Xty))
    }, error = function(e) { ok <<- FALSE; NULL })
    if (!ok) beta <- solve(XtX, Xty)
  }
  as.numeric(beta)
}

screen_topk <- function(X, y, k) {
  p <- ncol(X)
  k <- max(1, min(k, p))
  # center X (no scaling) for marginal correlation ranking
  xc <- scale(X, scale = FALSE)
  sc <- as.numeric(crossprod(xc, y)) # proportional to corr if y centered
  idx <- order(abs(sc), decreasing = TRUE)[seq_len(k)]
  idx
}

screen_topk_mic <- function(X, y, k) {
  p <- ncol(X)
  k <- max(1, min(k, p))
  
  y <- as.numeric(y)
  
  # 1) Pearson corr with X
  sc_lin <- suppressWarnings(cor(X, y))  # returns length p
  
  # 2) Pearson corr with X^2 (quadratic effect)
  X2 <- X * X
  sc_quad <- suppressWarnings(cor(X2, y))
  
  # 3) Spearman corr (monotone nonlinear)
  sc_spear <- suppressWarnings(cor(X, y, method = "spearman"))
  
  # replace NA with 0 if any (e.g. constant columns)
  sc_lin[is.na(sc_lin)]     <- 0
  sc_quad[is.na(sc_quad)]   <- 0
  sc_spear[is.na(sc_spear)] <- 0
  
  # final score: max of absolute correlations across the three
  score <- pmax(abs(sc_lin), abs(sc_quad), abs(sc_spear))
  
  idx <- order(score, decreasing = TRUE)[seq_len(k)]
  idx
}

screen_topk_nonlinear <- function(X, y, k,
                                  method = c("spearman", "pearson")) {
  method <- match.arg(method)
  p <- ncol(X)
  k <- max(1, min(k, p))
  
  sc <- apply(X, 2, function(x) {
    # robust to constant columns
    if (sd(x) == 0 || sd(y) == 0) return(0)
    suppressWarnings(cor(x, y, method = method))
  })
  
  idx <- order(abs(sc), decreasing = TRUE)[seq_len(k)]
  idx
}

get_beta_init <- function(method,
                          X_s, y_s,
                          X_t = NULL, y_t = NULL,
                          # params for ols_screen
                          k = NULL,        # if NULL, pick a data-driven default
                          screen_domain = c("source", "target", "pooled"),
                          # params for noise
                          noise_sigma = 0.05,
                          noise_base = c("lasso","ols_screen","ols"),
                          seed = NULL,
                          # params for ridge
                          ridge_lambda = 1e-2) {
  # if you want per-call RNG control, uncomment:
  # if (!is.null(seed)) set.seed(seed)
  
  screen_domain <- match.arg(screen_domain)
  noise_base    <- match.arg(noise_base)
  
  p  <- ncol(X_s)
  n_s <- nrow(X_s)
  
  lasso_beta <- function(X, y) {
    Lasso_func(X, y)$beta
  }
  
  if (method == "lasso") {
    return(lasso_beta(X_s, y_s))
  }
  
  if (method == "ols") {
    return(stable_ols(X_s, y_s))
  }
  
  if (method == "ridge") {
    # simple closed-form ridge on source: (X'X + λ I)^{-1} X'y
    XtX <- crossprod(X_s)        # p x p
    Xty <- crossprod(X_s, y_s)   # p x 1
    beta <- solve(XtX + ridge_lambda * diag(p), Xty)
    return(as.numeric(beta))
  }
  
  if (method == "ols_screen_nonlinear") {
    # choose data to screen on (usually source)
    Xd <- switch(screen_domain,
                 source = X_s,
                 target = { stopifnot(!is.null(X_t), !is.null(y_t)); X_t },
                 pooled = { stopifnot(!is.null(X_t), !is.null(y_t)); rbind(X_s, X_t) })
    yd <- switch(screen_domain,
                 source = y_s,
                 target = { stopifnot(!is.null(y_t)); y_t },
                 pooled = { stopifnot(!is.null(y_t)); c(y_s, y_t) })
    
    if (is.null(k)) {
      n_eff <- nrow(Xd)
      k_raw <- floor(n_eff / max(1, log(max(p, 2))))
      k <- min(
        max(10, k_raw),
        floor(p / 2),
        max(1, floor(n_s / 2))
      )
    }
    
    # nonlinear (monotone) screening via Spearman
    idx <- screen_topk_nonlinear(Xd, yd, k, method = "spearman")
    
    beta <- rep(0, p)
    beta[idx] <- stable_ols(X_s[, idx, drop = FALSE], y_s)
    return(beta)
  }
  
  if (method == "ols_screen_mic") {
    # choose data for screening
    Xd <- switch(screen_domain,
                 source = X_s,
                 target = { stopifnot(!is.null(X_t), !is.null(y_t)); X_t },
                 pooled = { stopifnot(!is.null(X_t), !is.null(y_t)); rbind(X_s, X_t) })
    yd <- switch(screen_domain,
                 source = y_s,
                 target = { stopifnot(!is.null(y_t)); y_t },
                 pooled = { stopifnot(!is.null(y_t)); c(y_s, y_t) })
    
    if (is.null(k)) {
      n_eff <- nrow(Xd)
      k_raw <- floor(n_eff / max(1, log(max(p, 2))))
      k <- min(
        max(10, k_raw),
        floor(p / 2),
        max(1, floor(n_s / 2))
      )
    }
    
    idx <- screen_topk_mic(Xd, yd, k)
    
    beta <- rep(0, p)
    beta[idx] <- stable_ols(X_s[, idx, drop = FALSE], y_s)
    return(beta)
  }
  
  if (method == "ols_screen") {
    # choose data for screening (ranking step)
    Xd <- switch(screen_domain,
                 source = X_s,
                 target = { stopifnot(!is.null(X_t), !is.null(y_t)); X_t },
                 pooled = { stopifnot(!is.null(X_t), !is.null(y_t)); rbind(X_s, X_t) })
    yd <- switch(screen_domain,
                 source = y_s,
                 target = { stopifnot(!is.null(y_t)); y_t },
                 pooled = { stopifnot(!is.null(y_t)); c(y_s, y_t) })
    
    if (is.null(k)) {
      # data-driven default:
      # - n_eff for screening complexity
      # - n_s for stability of OLS on X_s
      n_eff  <- nrow(Xd)
      k_raw  <- floor(n_eff / max(1, log(max(p, 2))))   # n / log p style
      k <- min(
        max(10, k_raw),               # at least 10, or n/log p
        floor(p / 2),                 # don't use more than half of all vars
        max(1, floor(n_s / 2))        # keep OLS dimension <= n_s/2
      )
    }
    
    # screen top-k by marginal correlation
    idx <- screen_topk(Xd, yd, k)
    
    # fit OLS on source restricted to screened set, then embed
    beta <- rep(0, p)
    beta[idx] <- stable_ols(X_s[, idx, drop = FALSE], y_s)
    return(beta)
  }
  
  # if (method == "noise") {
  #   # build base init then perturb (currently always from Lasso on source)
  #   lasso_initial <- Lasso_func(X_s, y_s)
  #   base_beta     <- lasso_initial$beta
  #   supp          <- which(abs(base_beta) > 0.1)
  #   beta_noisy    <- base_beta
  #   if (length(supp) > 0) {
  #     beta_noisy[supp] <- base_beta[supp] + rnorm(length(supp), 0, noise_sigma)
  #   }
  #   return(as.numeric(beta_noisy))
  # }
  
  if (method == "noise") {
    # build base init then perturb
    lasso_initial = Lasso_func(X_s,y_s)
    base_beta = lasso_initial$beta
    idx <- sample(seq_len(p), 30, replace = FALSE)
    beta_noisy <- base_beta
    beta_noisy[idx] <- base_beta[idx] + rnorm(30, 0, noise_sigma)
    return(as.numeric(beta_noisy))
  }
  
  
  stop("Unknown init method: ", method)
}


burn_rng_like_original_ours <- function(X_s, y_s, X_t, y_t, trim_form = "standard", rho = 0.5) {
  n_s <- nrow(X_s); n_t <- nrow(X_t)
  Q_s <- Q_trim(X_s, rho = rho, trim_form = trim_form)
  Q_t <- Q_trim(X_t, rho = rho, trim_form = trim_form)
  QX_s <- Q_s %*% X_s; Qy_s <- Q_s %*% y_s
  QX_t <- Q_t %*% X_t; Qy_t <- Q_t %*% y_t
  
  # (1) beta_s (trimmed source)
  invisible(Lasso_func(QX_s, Qy_s)$beta)
  
  # (2) lasso_initial on source (untrimmed)
  lasso_initial <- Lasso_func(X_s, y_s)
  beta_initial  <- lasso_initial$beta
  Z_s           <- lasso_initial$res
  
  # deterministic Z_t mapping (no RNG)
  Z_t <- find_Zt(t(X_t)/n_t, t(X_s) %*% Z_s / n_s)$Z_t
  
  # (3) eta on target around beta_initial (untrimmed)
  invisible(Lasso_func(X_t, y_t - Z_t - X_t %*% beta_initial)$beta)
}

with_preserve_seed <- function(expr) {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  eval.parent(substitute(expr))
}