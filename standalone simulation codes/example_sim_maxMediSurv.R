##### Here the number of predictors (p) is set as 100; 
##### the number of MC simulation runs (num_r) is set as 1;
##### the number of random data ordering (num_rdo) is set as 1, 
##### so as to deliver a quick demonstration on a standalone local computer.
##### To reproduce our numeric results, it's required that p = c(100, 1000, 1e4, 1e5, 1e6), 
##### num_r = 1000, and num_rdo = 1 or 10.
##### We highly suggest that parallel computing techniques should be used.


##### Configuration that can be ignored when for data analysis
num_r = 1;
# Censoring rate: '10%' or '30%'
cr = '10%'
# List of (n,p) values to run simulation
# for full simulation: list(c(500, 100), c(500, 1000), c(500, 1e4), c(500, 1e5), c(500, 1e6))
np_list = list(c(200, 1e5))
# List of correlation values between predictors
rho_val = 0.75
# List of data generating distributions to use when running simulations
model_list = c('A1.IE') #'N.IE','A2.IE','N.DE','A1.DE','A2.DE'


##### Configuration needed 
num_rdo = 1
# Estimation Index
est_index_val = 'whole samples'
# Methods to run
## 'SOSE' is our method; the rest are competing methods though not perfect ones.
meths_vec = c('SOSE') #'Naive_OSE', 'BONF_OSE' 
# Tau
quar = 0.9
# Significance level
alpha_val = 0.05
# elln 
elln_part = c(2); dim_elln = length(elln_part)
# Number of digits used here
num_digits = 7; options(digits = num_digits)


# Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'mvtnorm', 'dplyr', 'here')
for (package in package_list) {
    require(package, character.only = TRUE); library(package, character.only = TRUE) }

# Read in source codes that should be in the same working directory as is this script	
here::i_am('example_sim_maxMediSurv.R')
source('code_maxMediSurv.R')


##### For data analysis
##### Below X = observed time; delta = censoring status; A = exposure/treatment; B = n X p matrix of pre-standardized mediators
##### n: number of subjects and p: number of mediators
dat = 'Input your data whose format is list(X = X, delta = delta, A = A, B = B)'
obj0 <- NumericalStudy$new(input_data = dat)

set.seed(2023)
n = 'Sample size of your data'; p = 'Number of mediators of your data'; rho = rho_val

res = data.frame(method = NA, n = NA, p = NA, quar = NA, elln = NA, 
                 est = NA, se = NA, lb_ci = NA, ub_ci = NA, p_val = NA)
    
out = list()
for (meth in meths_vec) {
  if (meth == 'SOSE') {
  
    #print(meth) 
    SOSE_est = do.call(rbind, lapply(seq(num_rdo), function(rdo_idx){
      inds = sample(1:n)
      dat0 = list(X = dat$X[1:n][inds], delta = dat$delta[1:n][inds], 
                  A = dat$A[1:n][inds], B = dat$B[1:n,1:p][inds,])
      obj1 <- NumericalStudy$new(input_data = dat0)
          
          do.call(rbind, lapply(1:length(elln_part), function(d_idx){
            d = elln_part[d_idx]; elln = ceiling(n/d)
            chunk_size = ceiling((n - elln)/10)
            SOSE_est1 = obj1$Stab_onestep_est(all_obs = 1:n, chunk_size, 
                                              elln, est_index = est_index_val, 
                                              alpha = alpha_val, num_top = 1, 
                                              quar_trunc = quar)
            return( t( c( SOSE_est1$ci, SOSE_est1$est, SOSE_est1$se, 
                          SOSE_est1$p_val, elln ) ) ) } ) )
        } ) )
    out$SOSE = SOSE_est
        
  } else if (meth == 'Naive_OSE') {
        
    #print(meth)
    Naive_OSE_est = obj0$Naive_onestep_est(alpha = alpha_val, 
                                           quar_trunc = quar, num_top = 1)
    out$Naive_OSE = c( Naive_OSE_est$est, Naive_OSE_est$se, Naive_OSE_est$p_val )
        
  } else if (meth == 'BONF_OSE') {
        
    #print(meth)
    BONF_OSE_est = obj0$Bonf_onestep_ests(alpha = alpha_val, 
                                          quar_trunc = quar)
    out$BONF_OSE = c( BONF_OSE_est$est, BONF_OSE_est$se, BONF_OSE_est$p_val )
        
  } else { stop('Invalid method.') }
    
  res = rbind(sim, do.call(rbind,lapply(meths_vec, function(meth) {
      if ( meth == 'SOSE' ) {
        result = data.frame( method = rep(meth, dim_elln), 
                             n = rep(n, dim_elln), p = rep(p, dim_elln),
                             quar = rep(quar, dim_elln), 
                             elln = out[[meth]][,6],
                             est = out[[meth]][,3], se = out[[meth]][,4],
                             lb_ci = out[[meth]][,1], ub_ci = out[[meth]][,2],
                             p_val = out[[meth]][,5] )
      } else {
        result = data.frame( method = meth, n = n, p = p, 
                             quar = quar, elln = NA, 
                             est = out[[meth]][1], se = out[[meth]][2],
                             lb_ci = NA, ub_ci = NA,
                             p_val = out[[meth]][3] )
      }
      return( result ) })))
}

res = res[-1,]
rownames(res) = 1:nrow(res)
# To save the result
save(res, file = 'dataanalysis_example_maxCorrSurv.Rdata')




##### For Simulation
source('sim_data_acquisition.R')

sim = data.frame(method = NA, model = NA, n = NA, p = NA, rho = NA, 
                 quar = NA, elln = NA, est = NA, se = NA, 
                 lb_ci = NA, ub_ci = NA, p_val = NA)

set.seed(2023)
np = np_list[[1]]
n = np[1]; p = np[2]; rho = rho_val

t0 <- proc.time()

for (r in 1:num_r) {  
  for (model in model_list) {

	  sim_dat = Simulate_data(n, p, model, censoring_rate = cr, rho = rho, 
	                      bin_prob = 0.5, A.B_effect = 0.2,
	                      censoring_dist_par = c(0.06, 0.13, 0.24, 0.38))
	  obj0 <- NumericalStudy$new(input_data = sim_dat)
	  #print(c(n, p, rho, model, cr))
	  
	  out = list()
	  for (meth in meths_vec) {
	    if (meth == 'SOSE') {
        
	      print(meth) 
	      SOSE_est = do.call(rbind, lapply(seq(num_rdo), function(rdo_idx){
	        inds = sample(1:n)
	        dat0 = list(X = sim_dat$X[1:n][inds], delta = sim_dat$delta[1:n][inds], 
	                    A = sim_dat$A[1:n][inds], B = sim_dat$B[1:n,1:p][inds,])
	        obj1 <- NumericalStudy$new(input_data = dat0)
	        
	        do.call(rbind, lapply(1:length(elln_part), function(d_idx){
	          d = elln_part[d_idx]; elln = ceiling(n/d)
	          chunk_size = ceiling((n - elln)/10)
	          SOSE_est1 = obj1$Stab_onestep_est(all_obs = 1:n, chunk_size, 
	                                            elln, est_index = est_index_val, 
	                                            alpha = alpha_val, num_top = 1, 
	                                            quar_trunc = quar)
	          return( t( c( SOSE_est1$ci, SOSE_est1$est, SOSE_est1$se, 
	                        SOSE_est1$p_val, elln ) ) ) } ) ) 
	        } ) )
	      out$SOSE = SOSE_est
		        
		  } else if (meth == 'OOSE') {
		          
		    print(meth)
		    OOSE_est = obj0$Oracle_onestep_est(alpha = alpha_val, quar_trunc = quar, 
		                                       idx_taken = 1)
		    out$OOSE = c( OOSE_est$est, OOSE_est$se, OOSE_est$p_val )
		          
		  } else if (meth == 'Naive_OSE') {
		    
		    print(meth)
		    Naive_OSE_est = obj0$Naive_onestep_est(alpha = alpha_val, 
		                                           quar_trunc = quar, num_top = 1)
		    out$Naive_OSE = c( Naive_OSE_est$est, Naive_OSE_est$se, Naive_OSE_est$p_val )
		    
		  } else if (meth == 'BONF_OSE') {
		          
		    print(meth)
		    BONF_OSE_est = obj0$Bonf_onestep_ests(alpha = alpha_val, 
		                                          quar_trunc = quar)
		    out$BONF_OSE = c( BONF_OSE_est$est, BONF_OSE_est$se, BONF_OSE_est$p_val )
		          
		  } else stop('Invalid method.') }
	  
	  sim = rbind(sim, do.call(rbind,lapply(meths_vec, function(meth) {
	    if ( meth == 'SOSE' ) {
	      result = data.frame( method = rep(meth, dim_elln), 
	                           model = rep(model, dim_elln),
	                           n = rep(n, dim_elln), p = rep(p, dim_elln), 
	                           rho = rep(rho, dim_elln), 
	                           quar = rep(quar, dim_elln), 
	                           elln = out[[meth]][,6],
	                           est = out[[meth]][,3], se = out[[meth]][,4],
	                           lb_ci = out[[meth]][,1], ub_ci = out[[meth]][,2],
	                           p_val = out[[meth]][,5] )
	    } else {
	      result = data.frame( method = meth, model = model, n = n, p = p, 
	                           rho = rho, quar = quar, elln = NA, 
	                           est = out[[meth]][1], se = out[[meth]][2],
	                           lb_ci = NA, ub_ci = NA,
	                           p_val = out[[meth]][3] )
	    }
	    return( result ) })))
  }
}

proc.time() - t0

sim = sim[-1,]
rownames(sim) = 1:nrow(sim)
# To save the simulation result
save(sim, file = 'simexample_maxMediSurv.Rdata')


##### sim contains the simulation output. Each row contains:
# method: Method run
# model: Data generating distributions
# n: Sample size
# p: Covariate dimension
# rho: Designed correlation of covariates
# quar: The quantile of outcomes used for truncation
# elln: Value of elln used by our method
# est: Estimate of the parameter
# se: Standard error of the statistic
# lb_ci: Lower bound of the confidence interval of the parameter
# ub_ci: Upper bound of the confidence interval of the parameter
# p-value: P-value of rejecting the null of no correlation
