##### This script is for the real data analysis that uses 'BONF_OSE' and 'SOSE'
##### with one random ordering of data, that is, r = 1. 
##### Further values of r depend on how many random orderings will be taken to generate the results for 'SOSE',
##### say X-fold stabilized one-step that uses r = X. In the real data application of this article, X = 100.


##### Configuration needed 
num_rdo = 1
# Estimation Index
est_index_val = 'whole samples'
# Tau
quar = 0.9
# Significance level
alpha_val = 0.1
# elln 
elln_part = c(1.2, 1.5, 1.8, 2, 2.5); dim_elln = length(elln_part)
# Number of digits used here
num_digits = 7; options(digits = num_digits)

# Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'mvtnorm', 'dplyr', 'here')

for (package in package_list) {
  require(package, character.only = TRUE); library(package, character.only = TRUE) }

# In this script, we set r = 1.
# r = 1
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
r <- as.numeric(slurm_arrayid)
set.seed(r)

# Read in source codes and your desired data set (and change the path accordingly) 
# that should be in the same working directory as is this script	
here::i_am('TCGA_analysis.R')
source('data_management.R')
source('code_maxMediSurv.R')
data_path = 'TCGA lungcancer data.RData'

# Methods to run
if (r == 1) { meths_vec = c('SOSE', 'BONF_OSE', 'SOSE.ext', 'BONF_OSE.ext')
} else { meths_vec = c('SOSE', 'SOSE.ext') }

# Detected censoring rate
cr = '60%'

dat <- get(load(data_path))
n <- nrow(dat$B); p <- ncol(dat$B)
print(c(n, p))

dat$U <- dat$U %>% select(age, stage_reduced_f)

res = data.frame(method = NA, n = NA, p = NA, quar = NA, elln = NA, 
                 est = NA, se = NA, lb_ci = NA, ub_ci = NA, p_val = NA)

res_idx = data.frame(method = NA, n = NA, p = NA, quar = NA, elln = NA,
                     kest = NA)

out = list(); est_idx = list()
for (meth in meths_vec) { print(meth) 
  if (meth == 'SOSE') {
    
    SOSE_est = do.call(rbind, lapply(seq(num_rdo), function(rdo_idx){
      inds = sample(1:n)
      dat0 = list(X = dat$X[1:n][inds], delta = dat$delta[1:n][inds], 
                  A = dat$A[1:n][inds], B = dat$B[1:n,1:p][inds,])
      
      temp = do.call(rbind, lapply(1:length(elln_part), function(d_indx){ print(d_indx)
        d = elln_part[d_indx]; elln = ceiling(n/d)
        chunk_size = ceiling((n - elln)/5)
        SOSE_est1 = Stab_onestep_est(input_data = dat0, all_obs = 1:n, 
                                     chunk_size, elln, 
                                     est_index = est_index_val, 
                                     alpha = alpha_val, num_top = 1, 
                                     quar_trunc = quar, ext_flag = F,
                                     as.weights = T)
        
        return( t( c( SOSE_est1$ci, SOSE_est1$est, SOSE_est1$se, SOSE_est1$p_val, elln ) ) ) 
      } ) )
      return( temp )
    } ) )
    out$SOSE = SOSE_est
    
    kest = do.call(rbind, lapply(seq(num_rdo), function(rdo_idx){
      inds = sample(1:n)
      dat0 = list(X = dat$X[1:n][inds], delta = dat$delta[1:n][inds], 
                  A = dat$A[1:n][inds], B = dat$B[1:n,1:p][inds,])
      
      do.call(rbind, lapply(1:length(elln_part), function(d_indx){ print(d_indx)
        d = elln_part[d_indx]; elln = ceiling(n/d)
        chunk_size = ceiling((n - elln)/5)
        kest1 = Stab_onestep_kestseq(input_data = dat0, all_obs = 1:n, chunk_size, 
                                     elln, num_top = 1, quar_trunc = quar, ext_flag = F)
        
        return( data.frame( 'elln' = elln, 'est_idx' = kest1 ) ) 
      } ) )
    } ) )
    est_idx$SOSE = kest
    
  } else if (meth == 'SOSE.ext') {
    
    SOSE.ext_est = do.call(rbind, lapply(seq(num_rdo), function(rdo_idx){
      inds = sample(1:n)
      dat0 = list(X = dat$X[1:n][inds], delta = dat$delta[1:n][inds], 
                  A = dat$A[1:n][inds], B = dat$B[1:n,1:p][inds,],
                  U = dat$U[1:n,][inds,])
      
      do.call(rbind, lapply(1:length(elln_part), function(d_indx){ print(d_indx)
        d = elln_part[d_indx]; elln = ceiling(n/d)
        chunk_size = ceiling((n - elln)/5)
        SOSE.ext_est1 = Stab_onestep_est(input_data = dat0, all_obs = 1:n, chunk_size, 
                                         elln, est_index = est_index_val, 
                                         alpha = alpha_val, num_top = 1, 
                                         quar_trunc = quar, ext_flag = T,
                                         as.weights = T)
        return( t( c( SOSE.ext_est1$ci, SOSE.ext_est1$est, SOSE.ext_est1$se, 
                      SOSE.ext_est1$p_val, elln ) ) ) 
      } ) )
    } ) )
    out$SOSE.ext = SOSE.ext_est
    
    kest.ext = do.call(rbind, lapply(seq(num_rdo), function(rdo_idx){
      inds = sample(1:n)
      dat0 = list(X = dat$X[1:n][inds], delta = dat$delta[1:n][inds], 
                  A = dat$A[1:n][inds], B = dat$B[1:n,1:p][inds,],
                  U = dat$U[1:n,][inds,])
      
      do.call(rbind, lapply(1:length(elln_part), function(d_indx){ print(d_indx)
        d = elln_part[d_indx]; elln = ceiling(n/d)
        chunk_size = ceiling((n - elln)/5)
        kest1 = Stab_onestep_kestseq(input_data = dat0, all_obs = 1:n, chunk_size, 
                                     elln, num_top = 1, quar_trunc = quar, ext_flag = T)
        
        return( data.frame( 'elln' = elln, 'est_idx' = kest1 ) ) 
      } ) )
    } ) )
    est_idx$SOSE.ext = kest.ext
    
  } else if (meth == 'BONF_OSE') {
    
    BONF_OSE_est = Bonf_onestep_est(input_data = dat, alpha = alpha_val, 
                                      quar_trunc = quar, ext_flag = F, num_top = 1)
    out$BONF_OSE = c( BONF_OSE_est$ci, BONF_OSE_est$est, BONF_OSE_est$se, 
                      BONF_OSE_est$p_val )
    
  } else if (meth == 'BONF_OSE.ext') {
    
    BONF_OSE.ext_est = Bonf_onestep_est(input_data = dat, alpha = alpha_val, 
                                          quar_trunc = quar, ext_flag = T, num_top = 1)
    out$BONF_OSE.ext = c( BONF_OSE.ext_est$ci, BONF_OSE.ext_est$est, 
                          BONF_OSE.ext_est$se, BONF_OSE.ext_est$p_val )
  }  else { stop('Invalid method.') }
  
  
  if ( meth %in% c('SOSE', 'SOSE.ext') ) {
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
                           est = out[[meth]][3], se = out[[meth]][4],
                           lb_ci = out[[meth]][1], ub_ci = out[[meth]][2],
                           p_val = out[[meth]][5] )
  }
  res = rbind(res, result)
  
  if (meth %in% c('SOSE', 'SOSE.ext')) {
      result = data.frame( method = rep(meth, dim_elln), 
                           n = rep(n, dim_elln), p = rep(p, dim_elln),
                           quar = rep(quar, dim_elln), 
                           elln = est_idx[[meth]][,1],
                           kest = est_idx[[meth]][,2] )
      res_idx = rbind(res_idx, result) }
}

res = res[-1,]
rownames(res) = 1:nrow(res)

res_idx = res_idx[-1,]
rownames(res_idx) = 1:nrow(res_idx)

save(res, file = paste('TCGA_maxMediSurv_rdo', r, '.Rdata', sep = ''))

save(res_idx, file = paste('TCGA_maxMediSurv_predictors_rdo', r, '.Rdata', sep = ''))
