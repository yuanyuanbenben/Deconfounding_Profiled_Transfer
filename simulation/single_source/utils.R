# packages

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
  result = list("betas" = beta,"res" = y-X %*% beta)
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
    
    cat("Optimal Zt:\n")
    print(beta_opt)
    cat("\nMinimum infinity norm:", t_opt, "\n")
  } else {
    cat("Optimization failed. Status code:", lp_result$status, "\n")
    cat("Message:", lp_result$message, "\n")
  }
  return(list('Z_t'=beta_opt,'loss'=t_opt))
}

# data
generate_dataset = function(Psi,phi,beta,n=300, p=500, s=5, q=3, sigmaE=2, sigma=2, pert=1){
  #create H and Gamma with N(0,1) values and of appropriate size. H can be tuned with pert
  H = pert*matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE)
  #value of X independent from H
  E = matrix(rnorm(n*p,mean=0,sd=sigmaE),n,p,byrow = TRUE)
  
  #defined in eq. (2), high-dimensional measured covariates
  X = E + H %*% Psi 
  
  
  # beta = matrix(rep(c(1,0),times = c(s,p-s)),p,1,byrow = TRUE)
  
  #nx1 matrix with values of mean 0 and SD of sigma, error in Y independent of X
  nu = matrix(rnorm(n*1,mean=0,sd=sigma),n,1,byrow = TRUE)
  
  #eq. (1), the response of the Structural Equation Model
  Y = X %*% beta + H %*% phi + nu
  return_list = list("X"= X,"Y"= Y,"Hphi"=H%*%phi)
  return(return_list)
}



