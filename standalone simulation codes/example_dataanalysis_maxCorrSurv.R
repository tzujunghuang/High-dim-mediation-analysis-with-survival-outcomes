##### Here the number of predictors (p) is set as 100; 
##### the number of MC simulation runs (num_r) is set as 1;
##### the number of random data ordering (num_rdo) is set as 1, 
##### so as to deliver a quick demonstration on a standalone local computer.
##### To reproduce our numeric results, it's required that p = c(100, 1000, 1e4, 1e5, 1e6), 
##### num_r = 1000, and num_rdo = 1 or 10.
##### We highly suggest that parallel computing techniques should be used.


##### Configuration needed 
num_rdo = 1
# Estimation Index
est_index_val = 'whole samples'
# Methods to run
## 'SOSE' is our method; the rest are competing methods though not perfect ones.
meths_vec = c('SOSE')  #'Naive_OSE', 'BONF_OSE'
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
here::i_am('example_sim_maxCorrSurv.R')
source('code_maxCorrSurv.R')


##### For data analysis
##### Below X = observed time; delta = censoring status; A = exposure/treatment; B = n X p matrix of pre-standardized mediators
##### n: number of subjects and p: number of mediators
dat = 'Input your data whose format is list(X = X, delta = delta, A = A, B = B)'
#dat = get(load('test.RData'))
dat$A = as.vector(dat$A)
obj0 <- NumericalStudy$new(input_data = dat)

set.seed(2023)
n = 'Sample size of your data'; p = 'Number of mediators of your data';
#n = dim(dat$B)[1]; p = dim(dat$B)[2]

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

