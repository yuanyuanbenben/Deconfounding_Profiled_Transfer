source("utils.R")
# source("estimate_function.R")
source("estimate_function_fixed_p.R")


# setting
n_s = 600
n_t = 100
p = 10
q = 3
s = 4
s_delta = 1
sigma = 1
# control the cofounder structure
alpha_Psi = 1
# control the effect of confounder to model
alpha_phi = 1
# trim operator
trim_form = 'standard' # 'zero', 'no_trim'

# parameters
beta_s = matrix(rep(c(2,0),times = c(s,p-s)),p,1,byrow = TRUE)
eta = matrix(rep(c(-4,0,2),times = c(s_delta,p-2*s_delta,s_delta)),p,1,byrow = TRUE)
beta_t = beta_s + eta

# repeated experiments
our_eta = matrix(0,p,100)
baseline1_eta = matrix(0,p,100)
baseline2_eta = matrix(0,p,100)

our_beta = matrix(0,p,100)
baseline1_beta = matrix(0,p,100)
baseline2_beta = matrix(0,p,100)



for (repeat_index in 1:100){
  print(repeat_index)
  set.seed(repeat_index+1234)
  
  # confounder structure 
  # we can relax the spike condition for confounder structure
  Psi = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE) / n_s**0.5 * log(n_s) /2 * alpha_Psi
  # Psi = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE) 
  phi = matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE) * alpha_phi
  # alpha=0.8
  # phi = phi_*alpha + matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE) * alpha_phi * (1-alpha)
  ## source data
  source = generate_dataset(Psi,phi,beta_s,n_s,p,s,q,sigmaE=2,sigma=sigma,pert=1)
  X_s = source$X
  y_s = as.vector(source$Y)
  
  ## target data
  target = generate_dataset(Psi,phi,beta_t,n_t,p,s,q,sigmaE=2,sigma=sigma,pert=1)
  X_t = target$X
  y_t = as.vector(target$Y)
  
  p = 10
  q = 3
  s_s = 4
  s_t = 5
  
  # our method 
  # fixed p
  result = transfer_estimation_fixed_p(X_s,y_s,X_t,y_t,s_s,q)
  our_eta[1:p,repeat_index] = result$eta
  our_beta[1:p,repeat_index] = result$beta_t
  # baseline 1
  # directly using target data
  baseline1_result = baseline1_estimation_fixed_p(X_s,y_s,X_t,y_t,s_s,s_t,q)
  baseline1_eta[1:p,repeat_index] = baseline1_result$eta
  baseline1_beta[1:p,repeat_index] = baseline1_result$beta_t
  # baseline 2
  # classical transfer learning
  baseline2_result = baseline2_estimation_fixed_p(X_s,y_s,X_t,y_t,s_s,q)
  baseline2_eta[1:p,repeat_index] = baseline2_result$eta
  baseline2_beta[1:p,repeat_index] = baseline2_result$beta_t
}


loss_eta_l2 = colSums((our_eta-c(eta))**2)**0.5
loss_eta_baseline1_l2 = colSums((baseline1_eta-c(eta))**2)**0.5
loss_eta_baseline2_l2 = colSums((baseline2_eta-c(eta))**2)**0.5
loss_eta_l1 = colSums(abs(our_eta-c(eta)))
loss_eta_baseline1_l1 = colSums(abs(baseline1_eta-c(eta)))
loss_eta_baseline2_l1 = colSums(abs(baseline2_eta-c(eta)))

loss_beta_l2 = colSums((our_beta-c(beta_t))**2)**0.5
loss_beta_baseline1_l2 = colSums((baseline1_beta-c(beta_t))**2)**0.5
loss_beta_baseline2_l2 = colSums((baseline2_beta-c(beta_t))**2)**0.5
loss_beta_l1 = colSums(abs(our_beta-c(beta_t)))
loss_beta_baseline1_l1 = colSums(abs(baseline1_beta-c(beta_t)))
loss_beta_baseline2_l1 = colSums(abs(baseline2_beta-c(beta_t)))

write.csv(data.frame(loss_eta_l2,loss_eta_baseline1_l2,loss_eta_baseline2_l2,loss_eta_l1,loss_eta_baseline1_l1,loss_eta_baseline2_l1,
                     loss_beta_l2,loss_beta_baseline1_l2,loss_beta_baseline2_l2,loss_beta_l1,loss_beta_baseline1_l1,loss_beta_baseline2_l1),'fix_p_single_source_linear_error.csv')

loss_data = read.csv('fix_p_single_source_linear_error.csv')

print(mean(loss_data$loss_eta_l2))
print(mean(loss_data$loss_eta_baseline1_l2))
print(mean(loss_data$loss_eta_baseline2_l2))
print(mean(loss_data$loss_beta_l2))
print(mean(loss_data$loss_beta_baseline1_l2))
print(mean(loss_data$loss_beta_baseline2_l2))

print(mean(loss_data$loss_eta_l1))
print(mean(loss_data$loss_eta_baseline1_l1))
print(mean(loss_data$loss_eta_baseline2_l1))
print(mean(loss_data$loss_beta_l1))
print(mean(loss_data$loss_beta_baseline1_l1))
print(mean(loss_data$loss_beta_baseline2_l1))

print(mean((loss_data$loss_eta_l2 - mean(loss_data$loss_eta_l2))**2)/100**0.5)
print(mean((loss_data$loss_eta_baseline1_l2 - mean(loss_data$loss_eta_baseline1_l2))**2)/100**0.5)
print(mean((loss_data$loss_eta_baseline2_l2 - mean(loss_data$loss_eta_baseline2_l2))**2)/100**0.5)
print(mean((loss_data$loss_beta_l2 - mean(loss_data$loss_beta_l2))**2)/100**0.5)
print(mean((loss_data$loss_beta_baseline1_l2 - mean(loss_data$loss_beta_baseline1_l2))**2)/100**0.5)
print(mean((loss_data$loss_beta_baseline2_l2 - mean(loss_data$loss_beta_baseline2_l2))**2)/100**0.5)


print(mean((loss_data$loss_eta_l1 - mean(loss_data$loss_eta_l1))**2)/100**0.5)
print(mean((loss_data$loss_eta_baseline1_l1 - mean(loss_data$loss_eta_baseline1_l1))**2)/100**0.5)
print(mean((loss_data$loss_eta_baseline2_l1 - mean(loss_data$loss_eta_baseline2_l1))**2)/100**0.5)
print(mean((loss_data$loss_beta_l1 - mean(loss_data$loss_beta_l1))**2)/100**0.5)
print(mean((loss_data$loss_beta_baseline1_l1 - mean(loss_data$loss_beta_baseline1_l1))**2)/100**0.5)
print(mean((loss_data$loss_beta_baseline2_l1 - mean(loss_data$loss_beta_baseline2_l1))**2)/100**0.5)