##### This script is for one single MC simulation run, say r = 1.
##### To reproduce our numeric results, it is required that r = 1,...,1000 and 
##### it is highly suggested that parallel computing techniques should be used.


##### Configuration
# Read in source codes that should be in the same working directory as is this script	
source('data_management.R')
source('sim_data_acquisition.R')
source('code_maxMediSurv.R')
# Estimation Index
est_index_val = 'whole samples'
# Censoring rate
cr = '20%'
# Methods to run
meths_vec = c('SOSE', 'BONF_OSE', 'OOSE') #'SOSE.ext', 'BONF_OSE.ext', 'OOSE.ext'
# List of (n,p) values to run simulation
np_list = list(c(800, 1e2), c(800, 1e3), c(800, 1e4), c(800, 1e5))
#np_list = list(c(800, 1e6))
# List of correlation values between predictors
rho_val = 0.1
# Tau
quar = 0.9
# Significance level
alpha_val = 0.1
# List of data generating distributions to use when running simulations
model_list = c('A1.IE','A2.IE','N.IE') #'A1.IE.ext','A2.IE.ext','N.IE.ext'
as.weights_val = TRUE

# True parameters
trueNIEs = data.frame(model = c('N.IE','A1.IE','A2.IE'), 
                      NIE = c(0,0.2,0.2))
#trueNIEs = data.frame(model = c('N.IE.ext','A1.IE.ext','A2.IE.ext'), 
#                      NIE = c(0,0.4,0.4))


##### In this script, we set r = 1 because it's for a single MC run.
# r = 1
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
r <- as.numeric(slurm_arrayid)
set.seed(r)

# Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'mvtnorm', 'dplyr') #'foreach', 'doParallel'
# new_packages = package_list[!(package_list %in% installed.packages()[,"Package"])]
# if (length(new_packages)) { install.packages(new_packages) } 
for (package in package_list) {
    require(package, character.only = TRUE); library(package, character.only = TRUE) }

# n_cores = detectCores() - 2 # 3 cores detected on local machine; 6 on each node of clusters
# cl = makeCluster(n_cores)
# registerDoParallel(cl) 

# Number of digits used here
num_digits = 7; options(digits = num_digits)


sim = data.frame(method = NA, model = NA, n = NA, p = NA, 
                 quar_trunc = NA, rej = NA, cover = NA, est = NA, se = NA)

for (model in model_list) {
  
  #trueNIEs = get(load('trueNIEs.Rdata'))
  Psi0 = as.numeric(trueNIEs[trueNIEs$model == model, 'NIE'])
  # elln 
  elln_part = ifelse(model %in% c('N.IE', 'N.IE.ext'), 1.2, 1.2); 
  dim_elln = length(elln_part)
  #true_data = Simulate_data(n = 10^5, p = 15, model, censoring_rate = cr, rho = rho_val,
  #                          mc_true_par = T)
  #Psi_vec = true_Psi_k(input_data = true_data, B_index = 1:15)
  #k0 = order(abs(Psi_vec), decreasing = TRUE)[1]
  #Psi0 = Psi_vec[k0]
  
	for (np in np_list) {
	  
	  n = np[1]; p = np[2];
	  dat = Simulate_data(n, p, model, censoring_rate = cr, rho = rho_val)
	  print(c(n, p, rho_val, model, cr))
	  
	  out = list()
	  for (meth in meths_vec) {
	    if (meth == 'SOSE') {
        
	      print(meth) 
	      SOSE_est = do.call(rbind, lapply(1:length(elln_part), function(d_indx){
	        d = elln_part[d_indx]; elln = ceiling(n/d)
	        chunk_size = ceiling((n - elln)/5)
	        SOSE_est = Stab_onestep_est(input_data = dat, all_obs = 1:n, chunk_size, 
	                                    elln, est_index = est_index_val, 
	                                    alpha = alpha_val, num_top = 1, 
	                                    quar_trunc = quar, ext_flag = F,
	                                    as.weights = as.weights_val)
	        return( t( c( SOSE_est$ci, SOSE_est$est, SOSE_est$se, elln ) ) ) 
	      } ) )
	      out$SOSE = SOSE_est
	      
	    } else if (meth == 'SOSE.ext') {
	      
	      print(meth) 
	      SOSEext_est = do.call(rbind, lapply(1:length(elln_part), function(d_indx){
	        d = elln_part[d_indx]; elln = ceiling(n/d)
	        chunk_size = ceiling((n - elln)/5)
	        SOSEext_est = Stab_onestep_est(input_data = dat, all_obs = 1:n, chunk_size, 
	                                       elln, est_index = est_index_val, 
	                                       alpha = alpha_val, num_top = 1, 
	                                       quar_trunc = quar, ext_flag = T,
	                                       as.weights = as.weights_val)
	        return( t( c( SOSEext_est$ci, SOSEext_est$est, SOSEext_est$se, elln ) ) ) 
	      } ) )
	      out$SOSEext = SOSEext_est
	   
		  } else if (meth == 'OOSE') {
		          
		    print(meth)
		    OOSE_est = Oracle_onestep_est(input_data = dat, alpha = alpha_val, 
		                                  quar_trunc = quar, ext_flag = F, idx_taken = 1)
		    out$OOSE = c( OOSE_est$rej, OOSE_est$est, OOSE_est$se, OOSE_est$ci )
		          
		  } else if (meth == 'OOSE.ext') {
		    
		    print(meth)
		    OOSEext_est = Oracle_onestep_est(input_data = dat, alpha = alpha_val, 
		                                     quar_trunc = quar, ext_flag = T, idx_taken = 1)
		    out$OOSEext = c( OOSEext_est$rej, OOSEext_est$est, OOSEext_est$se, OOSEext_est$ci )
		    
		  } else if (meth == 'BONF_OSE') {
		    
		    print(meth)
		    BONF_OSE_est = Bonf_onestep_est(input_data = dat, alpha = alpha_val, 
		                                      quar_trunc = quar, ext_flag = F, num_top = 1)
		    out$BONF_OSE = c( BONF_OSE_est$rej, BONF_OSE_est$est, BONF_OSE_est$se, 
		                      BONF_OSE_est$ci )
		          
		  } else if (meth == 'BONF_OSE.ext') {
		    
		    print(meth)
		    BONF_OSEext_est = Bonf_onestep_est(input_data = dat, alpha = alpha_val, 
		                                         quar_trunc = quar, ext_flag = T, num_top = 1)
		    out$BONF_OSEext = c( BONF_OSEext_est$rej, BONF_OSEext_est$est, 
		                         BONF_OSEext_est$se, BONF_OSEext_est$ci )
		          
		  } else stop('Invalid method.') }
	  gc()
	  
	  sim = rbind(sim, do.call(rbind,lapply(meths_vec, function(meth) {
	    if ( meth %in% c('SOSE', 'SOSE.ext') ) {
	      result = data.frame( method = rep(meth, dim_elln), 
	                           model = rep(model, dim_elln),
	                           n = rep(n, dim_elln), p = rep(p, dim_elln),
	                           quar_trunc = rep(quar, dim_elln),
	                           rej = 1*( out[[meth]][,1] > 0 | out[[meth]][,2] < 0 ),
	                           cover = 1*( out[[meth]][,1] <= Psi0 & Psi0 <= out[[meth]][,2] ),
	                           est = out[[meth]][,3], se = out[[meth]][,4] )
	    } else {
	      result = data.frame( method = meth, model = model, n = n, p = p, 
	                           quar_trunc = quar, 
	                           rej = 1*(out[[meth]][1] > 0), 
	                           cover = 1*( out[[meth]][4] <= Psi0 & Psi0 <= out[[meth]][5] ),
	                           est = out[[meth]][2], se = out[[meth]][3] )
	    }
	    return( result ) })))
  }
}

sim = sim[-1,]
rownames(sim) = 1:nrow(sim)
# Used to save the simulation result in the r-th MC run (here r = 1)
save(sim, file = paste('sim_maxMediSurv_cr', gsub("\\%", "", cr), '_', r, '.Rdata', sep = ''))


##### sim contains the simulation output. Each row contains:
# method: Method run
# model: Data generating distributions
# n: Sample size
# p: Covariate dimension
# quar_trunc: The quantile of outcomes used for truncation
# rej: Whether rejecting the null of no mediation effects
# cover: Whether the c.i. covers the true parameter
# est: Estimate of the parameter
# se: Standard error of the statistic

