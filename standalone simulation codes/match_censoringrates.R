get_censoringrates = function(n, p, model, censoring_rate, rho, bin_prob = 0.5,
                              A.B_effect = 0.2,
                              censoring_dist_par = c(0.06, 0.13, 0.24, 0.38), 
                              just.first100 = FALSE, mc_true_par = FALSE, 
                              block_size = 1e4){
  
  # Simulate data of predictor structure
  A = rbinom(n, 1, bin_prob)
  
  # For any p, especially p > 1000
  if (rho == 0) {
    B = matrix(rnorm(n*p), ncol = p) - matrix(A.B_effect*rep(A,p), ncol = p)
  } else {
    a = uniroot( function(x){ 1 - rho - (x - sqrt((1 - x^2)/(p - 1)))^2 }, c(0,1), 
                 tol = .Machine$double.eps )$root
    b = sqrt((1 - a^2)/(p - 1))
    p0 = ifelse(just.first100, 100, p)
    tmp = matrix(rnorm(n*p0), ncol = p0)
    B = (matrix( rep(rowSums(tmp), p0), ncol = p0 )*b + tmp[,1:p0]*(a - b)
         - matrix(A.B_effect*rep(A,p), ncol = p))
    rm(tmp) }
  
  # Updating p by the column number of U to accommodate the case of just.first100 = TRUE
  p0 = dim(B)[2]
  
  # Simulate data of error structure
  if (grepl('IE', model, fixed = TRUE)) {
    error_IE = rnorm(n, 0, 1)
  } else {  
    error_DE = rnorm(n, 0, 0.7*(abs(B[,1]) + 0.7)) } 
  
  # Simulate data of survival time based on the specified AFT model
  if (model == 'N.IE') {
    log_T = error_IE
  } else if (model == 'A1.IE') {
    log_T = 0.4*A + 0.25*B[,1] + error_IE; 
  } else if (model == 'A2.IE') {
    log_T = c( as.matrix(cbind(A,B[,1:10])) %*% matrix(c(0.4, rep(0.15,5),rep(-0.1,5)), ncol = 1) ) + error_IE
  } else if (model == 'N.DE') {
    log_T = error_DE
  } else if (model == 'A1.DE') {
    log_T = 0.4*A + 0.25*B[,1] + error_DE
  } else if (model == 'A2.DE') {
    log_T = c( as.matrix(cbind(A,B[,1:10])) %*% matrix(c(0.4, rep(0.15,5),rep(-0.1,5)), ncol = 1) ) + error_DE
  } else stop('Invalid model choice.')
  
  # Simulate data of censoring time
  if (censoring_rate == '0%') {
      C = 1e5
  } else{
     for (i in 1:length(censoring_dist_par)) { assign( paste('par', i, sep = ''), censoring_dist_par[i] )}
       cr_par = par1*(censoring_rate == '10%') + par2*(censoring_rate == '20%') + 
         par3*(censoring_rate == '30%') + par4*(censoring_rate == '40%')
       C = log(rexp(n, cr_par)) }
    
    X = round(pmin(log_T, C), num_digits); 
    delta = 1*(log_T <= C);
    
    return( round(1 - mean(delta), 2) )
}    

inputs = expand.grid(c('N.IE', 'A1.IE', 'A2.IE', 'N.DE', 'A1.DE', 'A2.DE'),
                      c('10%', '20%', '30%', '40%'), 
                      c(0, 0.75))

results = do.call(rbind, lapply(1:nrow(inputs), function(r){ 
  
  cr = get_censoringrates(n = 10000, p = 1000, model = inputs[r,1],
                          censoring_rate = inputs[r,2], rho = inputs[r,3])
  return( c(inputs[r,], cr) )
}))

results

