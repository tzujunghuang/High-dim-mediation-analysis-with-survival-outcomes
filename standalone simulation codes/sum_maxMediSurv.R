##### This script is used to collect the results from the parallel computing implemented on 1000 nodes.

##### Configuration
# Estimation Index
est_index_val = 'whole samples'
# Censoring rate
cr = '20%'; cr_num = 20
# Methods to run
meths_vec = c('BONF_OSE', 'OOSE', 'SOSE')
# List of (n,p) values to run simulation
np_list = list(c(800, 1e5),c(800, 1e4), c(800, 1e3), c(800, 1e2))
n_char = 800; p_char = 5
#np_list = list(c(800, 1e6)); n_char = 800; p_char = 6
# elln
elln_char = '08'
# List of correlation values between predictors
rho_val = 0.1; rho_char = '01'
# Significance level
alpha_val = 0.1; alpha_char = '01'
# List of data generating distributions to use when running simulations
model_list = c('N.IE','A1.IE','A2.IE') #'N.DE','N1.DE','N2.DE','A1.DE','A2.DE'

num_r = 1000
orig_seq = 1:num_r;

final = do.call(rbind, lapply(orig_seq, function(r_idx){
  file = paste('sim_maxMediSurv_cr', cr_num, '_', r_idx, '.Rdata', sep = '')
  load(file)
  if (r_idx %in% c(250, 500, 750, 1000)) { print(r_idx) } 
  
  temp_final = NULL   
  for (meth in meths_vec) {
    for (model in model_list) {
      for (np in np_list) {
        
        n = np[1]; p = np[2]    
        if (meth %in% c('SOSE')) {
          #print(meth)
          condition1 = ( sim$method == meth & sim$model == model & sim$n == n & 
                         sim$p == p )

          est = sim[condition1, 'est']; se =  sim[condition1, 'se']
          rej =  sim[condition1, 'rej']; cover =  sim[condition1, 'cover']
          results = data.frame(model = model, meth = meth, 
                               n = n, p = p,
                               est = est, se = se, rej = rej, cover = cover, 
                               replic_num = r_idx)

        } else if (meth %in% c('BONF_OSE2')) { 
          #print(meth)
          condition1 = ( sim$method == meth & sim$model == model & sim$n == n & 
                         sim$p == p )
          est = sim[condition1, 'est']; se = sim[condition1, 'se']
          rej =  sim[condition1, 'rej']
          results = data.frame(model = model, meth = meth, 
                               n = n, p = p, 
                               est = est, se = se, rej = rej, cover = NA, 
                               replic_num = r_idx)
 
        } else {
          #print(meth)
          condition1 = ( sim$method == meth & sim$model == model & sim$n == n &
                         sim$p == p )
          est = sim[condition1, 'est']; se = sim[condition1, 'se']
          rej =  sim[condition1, 'rej']; cover =  sim[condition1, 'cover']
          results = data.frame(model = model, meth = meth, 
                               n = n, p = p, 
                               est = est, se = se, rej = rej, cover = cover,
                               replic_num = r_idx)
        }    
        temp_final = rbind(temp_final, results)
      }      
    }        
  }  
  return(temp_final) }))


# if (anyNA(final$rej)) { print(r_idx) }  
final$est_index = est_index_val

save(final, file = paste('results_n', n_char, 'p', p_char, '_cr', cr_num, '_', est_index_val, 
                         '_alpha', alpha_char, '_rho', rho_char, '_elln', elln_char, 
                         '.Rdata', sep = ''))

