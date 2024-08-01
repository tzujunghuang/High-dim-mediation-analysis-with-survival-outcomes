num_digits = 7
source('data_management.R')


# Simulate data, according to the distributions described in the article
Simulate_data = function(n, p, model, censoring_rate, rho, bin_prob = 0.5,
                         mc_true_par = FALSE, block_size = 1e4){
  
  # Simulate data of predictor structure
  A = rbinom(n, 1, bin_prob)
  U = rbinom(n, 1, 0.4)

  if (model %in% c('N.IE', 'N.IE.ext')) {
    a = uniroot( function(x){ 1 - 0.5 - (x - sqrt((1 - x^2)/(p - 1)))^2 }, c(0,1), 
                 tol = .Machine$double.eps )$root
    b = sqrt((1 - a^2)/(p - 1))
    tmp = matrix(rnorm(n*p), ncol = p)
    B = matrix( rep(rowSums(tmp), p), ncol = p )*b + tmp*(a - b)
    rm(tmp)
    #B = matrix(rnorm(n*p, 0, 1), ncol = p)

  } else {
    a = uniroot( function(x){ 1 - rho - (x - sqrt((1 - x^2)/(p - 1)))^2 }, c(0,1), 
                 tol = .Machine$double.eps )$root
    b = sqrt((1 - a^2)/(p - 1))
    tmp = matrix(rnorm(n*p), ncol = p)
    E = matrix( rep(rowSums(tmp), p), ncol = p )*b + tmp*(a - b)
    rm(tmp)
    
    B = cbind( 1*A + matrix(rnorm(n*1), ncol = 1), 
               matrix(rep(0.6*A, 4), ncol = 4) + matrix(rnorm(n*4), ncol = 4),
               matrix(rep(0.3*A, 5), ncol = 5) + matrix(rnorm(n*5), ncol = 5), 
               E[,11:p] )
  }
  
  # Simulate data of error structure
  if (grepl('IE', model, fixed = TRUE)) {
    error_IE = rnorm(n, 0, 1)
  } else {  
    error_DE = rnorm(n, 0, 1.2*abs(B[,1])) } 
  
  # Simulate data of survival time based on the specified AFT model
  if (model == 'N.IE') {
    log_T = 0.2*A + error_IE
  } else if (model == 'A1.IE') {
    log_T = 0.4*A + 0.2*B[,1] + error_IE; 
  } else if (model == 'A2.IE') {
    log_T = c( as.matrix(cbind(A,B[,1:10])) %*% matrix(c(0.4,rep(0.2,5),rep(-0.1,5)), ncol = 1) ) + error_IE
  } else if (model == 'N.IE.ext') {
    log_T = 0.2*A - 0.1*U + error_IE
  } else if (model == 'A1.IE.ext') {
    log_T = 0.4*A - 0.1*U + 0.2*B[,1] + error_IE; 
  } else if (model == 'A2.IE.ext') {
    log_T = c( as.matrix(cbind(A,U,B[,1:10])) %*% matrix(c(0.4,-0.1,rep(0.2,5),rep(-0.1,5)), 
                                                         ncol = 1) ) + error_IE
  } else stop('Invalid model choice.')
  
  if (!(mc_true_par)) {
    # Simulate data of censoring time
    if (censoring_rate == '0%') {
      C = 1e5
    }
    else {
      if (model %in% c('N.IE', 'N.IE.ext')) {
        censoring_dist_par = c(0.06, 0.14, 0.26, 0.4)
      } else if (model %in% c('A1.IE', 'A1.IE.ext')) { 
        censoring_dist_par = c(0.05, 0.11, 0.2, 0.3) 
      } else if (model %in% c('A2.IE', 'A2.IE.ext')) { 
        censoring_dist_par = c(0.04, 0.08, 0.15, 0.25) 
      }
      
      for (i in 1:length(censoring_dist_par)) { 
        assign( paste('par', i, sep = ''), censoring_dist_par[i] ) }
      cr_par = par1*(censoring_rate == '10%') + par2*(censoring_rate == '20%') + 
               par3*(censoring_rate == '30%') + par4*(censoring_rate == '40%')
      C = log(rexp(n, cr_par)) }
    
    X = round(pmin(log_T, C), num_digits); 
    delta = 1*(log_T <= C);
    
    # Mediator pre-standardization and saved for further analyses
    B_colSDs = colVars(B, sd_use = TRUE)
    B = colStandardization(B, B_colSDs)
    # B = (B - matrix(rep(colMeans(B), n), nrow=n, byrow=TRUE))/matrix(rep(B_colSDs, n), nrow=n, byrow=TRUE)
    input_data = list(X = X, delta = delta, A = A, B = B, U = U)
  } else {
    B_colSDs = colVars(B, sd_use = TRUE)
    B = colStandardization(B, B_colSDs)
    # B = (B - matrix(rep(colMeans(B), n), nrow=n, byrow=TRUE))/matrix(rep(B_colSDs, n), nrow=n, byrow=TRUE)
    log_T = round(log_T, num_digits)
    input_data = list(timeT = log_T, A = A, B = B, U = U)
  }
  return( input_data )
}    


true_Psi_k = function(input_data, B_index){
  timeT = input_data$timeT; A = input_data$A; 
  selectB = input_data$B[,B_index]; mean_A = mean(A)
  
  if (is.matrix(selectB)) {
    
    nn = nrow(selectB[A == 1,])
    a_term = mean(timeT[A == 1]) 
    
    mean_selectB = colMeans(selectB[A == 1,]);
    var_selectB = (colMeans((selectB[A == 1,])^2) - (mean_selectB)^2)*nn/(nn - 1)
    cov_Y_selectB = cov(selectB[A == 1,], timeT[A == 1])
    b_term = cov_Y_selectB / var_selectB; 
    b_term[is.nan(b_term) | is.na(b_term)] = 0
    
    mean_full_mat = matrix(rep(mean_selectB, dim(selectB)[1]), byrow = T, ncol = length(B_index))
    c_term = selectB - mean_full_mat
    
    Ehat_givenAB = a_term + matrix(rep(b_term, dim(selectB)[1]), byrow = T, ncol = length(B_index))*c_term
    eta_10 = colMeans((1 - A)*Ehat_givenAB/(1 - mean_A))
    
    return( mean(timeT[A == 1]) - eta_10 )
    
  } else {
    
    a_term = mean(timeT[A == 1]) 
    mean_selectB = mean(selectB[A == 1]); var_selectB = var(selectB[A == 1])
    cov_Y_selectB = mean((selectB[A == 1] - mean_selectB)*timeT[A == 1])
    b_term = cov_Y_selectB / var_selectB; 
    b_term[is.nan(b_term) | is.na(b_term)] = 0
    c_term = selectB - mean_selectB
    
    Ehat_givenAB = a_term + b_term*c_term
    eta_10 = mean((1 - A)*Ehat_givenAB/(1 - mean_A))
    
    return( mean(timeT[A == 1]) - eta_10 ) }
}

