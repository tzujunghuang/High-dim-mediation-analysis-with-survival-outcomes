num_digits = 7
source('data_management.R')

NumericalStudy <- setRefClass( "NumericalStudy",
 fields = list(
    input_data = "list",
    n = "numeric",
    p = "numeric"
  ),
  
 methods = list(
  initialize = function(input_data){
    input_data <<- input_data
    n <<- dim(input_data$B)[1]
    p <<- dim(input_data$B)[2]
  },
  
  subset_data = function(input_obs, input_Bindex){
    # input_obs is a vector of observation/individual indices
    # input_Bindex is the label of the selected mediator 
    df = cbind(input_data$X[input_obs], input_data$delta[input_obs], 
               input_data$A[input_obs], input_data$B[input_obs, input_Bindex])
    
    if (length(input_obs) > 1) {
      df = df[order(df[,1]),]
      X = df[,1]; delta = df[,2]; A = df[,3]; selectB = df[,4:dim(df)[2]]
    } else{ X = df[1]; delta = df[2]; A = df[3]; selectB = df[4:length(df)] }
    
    res = list(X = X, delta = delta, A = A, selectB = selectB)
    return( res )
  },
    
  KM_weight_func = function(obs){ # obs is a vector of observation/individual indices 
    data_km = data.frame(X = input_data$X[obs], delta = input_data$delta[obs])
    data_km = data_km[order(data_km$X),]
    prod = c(1, cumprod( sapply(1:(n - 1), function(i){ 
      ((n - i)/(n - i + 1))^(data_km$delta[i]) }) ) )
      kmwts = sapply(1:n, function(i){ data_km$delta[i]*prod[i]/(n - i + 1) } )
    return( kmwts )
  },  
    
  KM_SurF = function(t, obs){ # t is a scalar; obs is a vector of observation/individual indices 
    data_km = data.frame(X = input_data$X[obs], delta = input_data$delta[obs])
    km = survfit(Surv(X, 1-delta)~1, data = data_km)
    rm(data_km)
    
    survest = cbind(km$time, km$surv)
    if ( length(which(survest[, 1] <= t)) > 0 ) {
        return( survest[max(which(survest[, 1] <= t)), 2] )  
    } else { return( 1 ) }
  },  
  
  KM_SurF_self = function(obs, quar_trunc){ # obs is a vector of observation/individual indices 
    data_km = data.frame(X = input_data$X[obs], delta = input_data$delta[obs])
    tau = quantile(data_km$X, probs = quar_trunc)
    tau_surv = KM_SurF(tau, obs)
    
    km = survfit(Surv(X, 1-delta)~1, data = data_km)
    rm(data_km)
    
    km$surv[km$time > tau] = tau_surv
    survest = cbind(km$time, km$surv)
    return( survest )
  },
  
  Inverse_weight_func = function(x0, obs, quar_trunc, err_msg){
    # (x0, delta0) are the replication of subjects to be predicted
    # obs is a vector of observation/individual indices of x0
    tau = quantile(input_data$X[1:n], probs = quar_trunc)
    tau_surv = KM_SurF(tau, obs)

    if (length(x0) == 1) {
      inverse_weight = ifelse(x0 < tau, KM_SurF(x0, obs), tau_surv)
    } else {    
      if (length(unique(x0)) == length(x0)) { 
        KM_table = KM_SurF_self(obs, quar_trunc) 
      } else {
        KM_table0 = data.frame(KM_SurF_self(obs, quar_trunc))
        names(KM_table0) = c("time", "surv_prob")
        n_occur = data.frame(table(x0)); names(n_occur) = c("time", "freq")
        n_freq = n_occur$freq; n1 = n_freq[n_freq > 1]; 
        n1[is.nan(n1) | is.na(n1) | is.null(n1)] = 0
        
        tryCatch({ temp = data.frame(cbind(
          rep(KM_table0[KM_table0$time %in% n_occur$time[n_freq > 1],]$time, n1),
          rep(KM_table0[KM_table0$time %in% n_occur$time[n_freq > 1],]$surv_prob, n1))) },
          error = function(err_msg){
            message("Original error message:"); message(paste(err_msg,"\n", sep = ''))
            save(list(x0), file = paste('errordata_', cr, '_', r, '.Rdata', sep = ''))
            stop(err_msg) })
        names(temp) = names(KM_table0)
        # times_once = as.numeric( as.character(n_occur$time[n_freq == 1]) )
        KM_table = as.matrix(rbind(temp, 
          KM_table0[KM_table0$time %in% n_occur$time[n_freq == 1],])) }
        
      diff = setdiff(x0, KM_table[,1])  
      if (length(diff) > 0) {
        for (item in diff) { 
          KM_table = rbind(KM_table, 
                           c(item, KM_table[max(which(item >= KM_table[,1])), 2])) } }
      KM_table[KM_table[,1] > tau, 2] = tau_surv
      inverse_weight = as.vector(KM_table[order(KM_table[,1]),][,2]) }  
    
    if (length(inverse_weight) != length(x0)) {
      save(list(x0), file = paste('errordata_', cr, '_', r, '.Rdata', sep = ''))
      stop('Dimension Problem!')  
    } else { return( inverse_weight ) }
  },
  
  Est_Psi_k = function(B_index, obs, quar_trunc){
    # B_index is a scalar or a vector of mediator indices; 
    # obs is a vector of observation/individual indices used for estimation 
    df = subset_data(input_obs = obs, input_Bindex = B_index)
    X = df$X; delta = df$delta; A = df$A; selectB = df$selectB
    rm(df)
    
    inverse_weight = Inverse_weight_func(x0 = X, obs, quar_trunc, 
                                         err_msg = 'Error in Est_Psi0_d Inverse_weight KM_table')
    Y_Ghat = X*delta / inverse_weight; Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0
    
    Psi_k_1 = cov(selectB[A == 1], Y_Ghat[A == 1]) / var(selectB[A == 1])
    Psi_k_2 = mean(selectB[A == 1]) - mean(selectB[A == 0])
  
    return( Psi_k_1*Psi_k_2 )
  },
  
  get_Ehat_givenAX = function(A_val, B_val, obs_all, obs0, B_index, quar_trunc){
    # (A_val, B_val) are the given value of (A,B) this function to be evaluated at
    # obs_all is a vector of observation/individual indices used for estimation of G
    # obs0 is a vector of observation/individual indices used for estimation of the rest
    # B_index is a scalar or a vector of mediator indices
    X_all = input_data$X[obs_all];
    
    df0 = subset_data(input_obs = obs0, input_Bindex = B_index)
    X0 = df0$X; delta0 = df0$delta; A0 = df0$A; selectB0 = df0$selectB; rm(df0)
    
    tau = quantile(X_all, probs = quar_trunc)
    inverse_weight = Inverse_weight_func(x0 = X0, obs = obs0, quar_trunc, 
                                         err_msg = 'Error in IF_star_self Inverse_weight KM_table')
    
    time_comparison = outer(X0, X0, '>=')
    n_s = colSums(time_comparison*(A0 == A_val)) 
    Y_Ghat_mat = time_comparison*(A0 == A_val)*((X0*delta0) / inverse_weight) 
    Y_Ghat_mat[is.nan(Y_Ghat_mat) | is.na(Y_Ghat_mat)] = 0
    mean_YGhat_mat = sapply(1:length(n_s), function(j){
      ifelse( n_s[j] == 0, 0, 
              sum(Y_Ghat_mat[,j])/n_s[j] )}) 
    
    if ( is.matrix(selectB0) ) {
      
      selectB0_A = selectB0*(A0 == A_val)
      Ehat_pointwise_vals = do.call(rbind, lapply(1:length(X0), function(i){
        a_term = mean_YGhat_mat[i];
        if (i == length(X0)) { 
          selectB = t(as.matrix(selectB0_A[length(X0),]))
        } else { selectB = selectB0_A[i:length(X0),] }
        
      mean_selectB = sapply(1:ncol(selectB), function(col){
        ifelse( sum(selectB[,col] > 0) == 0, 0, 
                sum(selectB[,col])/sum(selectB[,col] > 0) )})
      var_selectB = sapply(1:ncol(selectB), function(col) {
        ifelse( sum(selectB[,col] > 0) == 0, 1, 
                sum((selectB[,col] - mean_selectB[col])^2)/sum(selectB[,col] > 0) )})
      cov_Y_selectB = sapply(1:ncol(selectB), function(col) {
        ifelse( sum(selectB[,col] > 0) == 0, 0, 
                sum((selectB[,col] - mean_selectB[col])*Y_Ghat_mat[i:length(X0),i])/sum(selectB[,col] > 0) )})
        
      b_term = cov_Y_selectB / var_selectB; 
      b_term[is.nan(b_term) | is.na(b_term)] = 0
      c_term = B_val - mean_selectB
      rm(selectB); rm(mean_selectB); rm(var_selectB); rm(cov_Y_selectB)
      return( a_term + b_term*c_term ) } ))
      Ehat_pointwise_vals[is.nan(Ehat_pointwise_vals) | is.na(Ehat_pointwise_vals)] = 0
    
    } else {
      selectB0_A = selectB0*(A0 == A_val)
      Ehat_pointwise_vals = sapply( 1:length(X0), function(i){
        a_term = mean_YGhat_mat[i];
        if (i == length(X0)) { 
          selectB = selectB0_A[length(X0)]
        } else{ selectB = selectB0_A[i:length(X0)] }
        
        mean_selectB = ifelse( sum(selectB > 0) == 0, 0, 
                               sum(selectB)/sum(selectB > 0) )
        
        if (length(selectB) > 1) { 
          var_selectB = ifelse( sum(selectB > 0) == 0, 1, 
                                sum((selectB - mean_selectB)^2)/sum(selectB > 0) ); 
          cov_Y_selectB = ifelse( sum(selectB > 0) == 0, 0, 
                                  sum((selectB - mean_selectB)*Y_Ghat_mat[i:length(X0),i])/sum(selectB > 0) )                          
          b_term = cov_Y_selectB / var_selectB
          b_term[is.nan(b_term) | is.na(b_term)] = 0
        } else { var_selectB = 'NA'; cov_Y_selectB = 0; b_term = 0 } 
        
        c_term = B_val - mean_selectB
        rm(selectB); rm(mean_selectB); rm(var_selectB); rm(cov_Y_selectB)
        return( a_term + b_term*c_term ) } )
      Ehat_pointwise_vals[is.nan(Ehat_pointwise_vals) | is.na(Ehat_pointwise_vals)] = 0 
    }
    return( Ehat_pointwise_vals )
  },
  
  get_Ehat_givenA = function(A_val, obs_all, obs0, obs1, B_index, quar_trunc){
    # A_val is the given value of A this function to be evaluated at
    # obs_all is a vector of observation/individual indices used for estimation of G
    # obs0 is a vector of observation/individual indices used for estimation of the rest
    # obs1 is a vector of observation/individual indices to be predicted
    # B_index is a scalar or a vector of mediator indices
    Xall = input_data$X[obs_all];
    
    df0 = subset_data(input_obs = obs0, input_Bindex = B_index)
    X0 = df0$X; delta0 = df0$delta; A0 = df0$A; selectB0 = df0$selectB; rm(df0)
    
    selectB1 = matrix(input_data$B[obs1, B_index], nrow = length(obs1), 
                      ncol = length(B_index))
    
    tau = quantile(Xall, probs = quar_trunc)
    inverse_weight = Inverse_weight_func(x0 = X0, obs = obs0, quar_trunc, 
                                         err_msg = 'Error in IF_star_self Inverse_weight KM_table')
    
    Y_Ghat = (X0*delta0 / inverse_weight)[A0 == A_val] 
    Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0; mean_YGhat = mean(Y_Ghat)
    a_term = mean_YGhat
    
    if ( is.matrix(selectB0) ) {
      selectB = selectB0[A0 == A_val,]; n_s = nrow(selectB)
      mean_selectB = colMeans(selectB)
      var_selectB = apply(selectB, 2, var)
      cov_Y_selectB = colMeans((selectB - matrix(rep(mean_selectB, n_s), 
                                                 byrow = T, nrow = n_s))*Y_Ghat)
      b_term = cov_Y_selectB / var_selectB; 
      b_term[is.nan(b_term) | is.na(b_term)] = 0
      b_term = matrix(b_term, byrow = T, ncol = ncol(selectB))
      c_term = do.call(cbind, lapply(1:ncol(selectB), 
                          function(col){ selectB1[,col] - mean_selectB[col] }))
      Ehat_pointwise_vals = a_term + do.call(cbind, lapply(1:ncol(selectB), 
                                        function(col){ b_term[,col]*c_term[,col] }))
      Ehat_pointwise_vals[is.nan(Ehat_pointwise_vals) | is.na(Ehat_pointwise_vals)] = 0
      rm(selectB); rm(mean_selectB); rm(var_selectB); rm(cov_Y_selectB);
    } else {
      selectB = selectB0[A0 == A_val]
      mean_selectB = mean(selectB); var_selectB = var(selectB)
      cov_Y_selectB = mean((selectB - mean_selectB)*Y_Ghat)
      b_term = cov_Y_selectB / var_selectB; 
      b_term[is.nan(b_term) | is.na(b_term)] = 0
      c_term = selectB1 - mean_selectB
      rm(selectB); rm(mean_selectB); rm(var_selectB); rm(cov_Y_selectB);
      Ehat_pointwise_vals = a_term + b_term*c_term
      Ehat_pointwise_vals[is.nan(Ehat_pointwise_vals) | is.na(Ehat_pointwise_vals)] = 0
    }
    return( Ehat_pointwise_vals )
  },
  
  get_martingale = function(obs, obs1, quar_trunc){
    # obs0 is a vector of observation/individual indices used for estimation of the rest
    # obs1 is a vector of observation/individual indices to be predicted
    X = input_data$X[obs]; delta = input_data$delta[obs]
    X1 = input_data$X[obs1]; delta1 = input_data$delta[obs1]
    tau = quantile(X, probs = quar_trunc)
    time_comparison1 = outer(X, X, '==') 
    time_comparison2 = outer(X, X, '>=')
    event_nums = colSums(time_comparison1*(1 - delta)); 
    risk_set = colSums(time_comparison2)
    hazard = event_nums/risk_set; 
    hazard[is.nan(hazard) | is.na(hazard)] = 0
    return( (X1 <= tau & delta1 == 0) - (X1 <= tau)*hazard )
  },
  
  IF_star_self = function(m, B_index, obs, quar_trunc){
    # m is a scalar; B_index is a scalar or a vector of mediator indices; 
    # obs is the vector of old/whole individuals indices both for estimation and to be predicted
    df = subset_data(input_obs = obs, input_Bindex = B_index)
    X = df$X; delta = df$delta; A = df$A; selectB = df$selectB; rm(df)
    
    m_mat = matrix(rep(m, length(obs)), byrow = TRUE, ncol = p)
    
    Psi_pars = sapply(B_index, function(k){ 
      Est_Psi_k(B_index = k, obs, quar_trunc) })
    
    tau = quantile(X, probs = quar_trunc)
    inverse_weight = Inverse_weight_func(x0 = X, obs, quar_trunc, 
                                         err_msg = 'Error in IF_star_self Inverse_weight KM_table')

    mean_A = mean(A)    
    Y_Ghat = (X*delta) / inverse_weight; Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0
    a = Y_Ghat*A/mean_A; a_mat = matrix(rep(a, p), ncol = p)
      
    Ehat_givenA = get_Ehat_givenA(A_val = 1, obs_all = obs, obs0 = obs, 
                                  obs1 = obs, B_index, quar_trunc)
    mart_X = get_martingale(obs, obs1 = obs, quar_trunc)
    
    if ( is.matrix(selectB) ) {
      reversed_oddsratios <- do.call(cbind, lapply(1:ncol(selectB), function(col){
        logit_model <- glm(A ~ selectB[,col], family = binomial(link = 'logit'))
        reversed_oddsratio <- sapply(selectB[,col], function(B){ 
          exp( -sum(coef(logit_model)*c(1,B)) ) })
        return( reversed_oddsratio ) }))
      
      b = do.call(cbind, lapply(1:ncol(selectB), function(col){ 
             A*reversed_oddsratios[,col]*(Y_Ghat - Ehat_givenA[,col])/(1 - mean_A) }))
      
      c = do.call(cbind, lapply(1:ncol(selectB), function(col){ 
             (1 - A)*Ehat_givenA[,col]/(1 - mean_A) }))
      
      get_integral = function(row){ 
        a_val = row[1]; b_val = row[2]; k_val = row[3]
        integrand = get_Ehat_givenAX(A_val = a_val, B_val = b_val, obs_all = obs, 
                                     obs0 = obs, B_index = k_val, quar_trunc)
        return( sum(integrand*mart_X) ) }
      
      ##### EXPENSIVE !!!
      d_integral = do.call(cbind, lapply(B_index, function(k){
        AB_mat = cbind(A, selectB[,k], k)
        integral = apply(AB_mat, 1, 'get_integral')
        return( integral ) }))
      
      d = do.call(cbind, lapply(1:ncol(selectB), function(col){ 
             return( (A*(1 - reversed_oddsratios[,col]*mean_A/(1 - mean_A))*
                      d_integral[,col])/mean_A ) }))
      
      Psi_pars_mat = matrix(rep(Psi_pars, length(obs)), byrow = TRUE, ncol = p)
      
      return( round( m_mat*(a_mat - b - c + d - Psi_pars_mat), 3) ) 
      
    } else { 
      logit_model <- glm(A ~ selectB, family = binomial(link = 'logit'))
      reversed_oddsratios <- sapply(selectB, function(B){ 
          exp( -sum(coef(logit_model)*c(1,B)) ) })
      
      b = as.vector( A*reversed_oddsratios*(Y_Ghat - Ehat_givenA)/(1 - mean_A) )
      c = as.vector( (1 - A)*Ehat_givenA/(1 - mean_A) )
      
      d_integral = do.call(rbind, lapply(1:length(X), function(j){ 
        # A rows by 1 cols
        integrand = get_Ehat_givenAX(A_val = A[j], B_val = selectB[j], obs_all = obs, 
                                     obs0 = obs, B_index, quar_trunc) 
        return( sum(integrand*mart_X) ) }))
      
      d = as.vector( A*(1 - reversed_oddsratios*mean_A/(1 - mean_A))*d_integral/mean_A )
      
      return( round( m*(a - b - c + d - Psi_pars), 3) ) 
    }  
  },
  
  IF_star = function(m, B_index, obs_all, obs0, obs1, quar_trunc){
    # m is a scalar; B_index is a scalar or a vector of mediator indices; 
    # (obs_all, obs0, obs1) are the vectors of whole, old and new individual indices
    df_all = subset_data(input_obs = obs_all, input_Bindex = B_index)
    X_all = df_all$X; delta_all = df_all$delta; rm(df_all)
    
    df0 = subset_data(input_obs = obs0, input_Bindex = B_index)
    X0 = df0$X; delta0 = df0$delta; A0 = df0$A; selectB0 = df0$selectB; rm(df0)
    
    df1 = subset_data(input_obs = obs1, input_Bindex = B_index)
    X1 = df1$X; delta1 = df1$delta; A1 = df1$A; selectB1 = df1$selectB; rm(df1)
      
    par_k = Est_Psi_k(B_index, obs = obs0, quar_trunc)
    
    inverse_weight1 = Inverse_weight_func(x0 = X1, obs = obs1, quar_trunc, 
                                          err_msg = 'Error in IF_star Inverse_weight KM_table')
    inverse_weight0 = Inverse_weight_func(x0 = X0, obs = obs0, quar_trunc, 
                                          err_msg = 'Error in IF_star Inverse_weight KM_table')
    mean_A0 = mean(A0)    
    Y_Ghat1 = (X1*delta1) / inverse_weight1; 
    Y_Ghat1[is.nan(Y_Ghat1) | is.na(Y_Ghat1)] = 0
    a = Y_Ghat1*A1/mean_A0;
    
    Y_Ghat0 = (X0*delta0) / inverse_weight0; 
    Y_Ghat0[is.nan(Y_Ghat0) | is.na(Y_Ghat0)] = 0
    b = mean(Y_Ghat0[A0 == 1])*A1/mean_A0
  
    logit_model0 <- glm(A0 ~ selectB0, family = binomial(link = 'logit'))
    reversed_oddsratios_selectB1 <- sapply(selectB1, function(B){ 
      exp( -sum(coef(logit_model0)*c(1,B)) ) })
    Ehat_givenA_selectB1 = get_Ehat_givenA(A_val = 1, obs_all, obs0, 
                                           obs1, B_index, quar_trunc)
    
    c = A1*reversed_oddsratios_selectB1*(Y_Ghat1 - Ehat_givenA_selectB1)/(1 - mean_A0)
    
    mart_X = get_martingale(obs = obs_all, obs1 = obs0, quar_trunc)
    
    integrand = do.call(rbind, lapply(1:length(X1), function(j){ 
      # A rows by 1 cols
      get_Ehat_givenAX(A_val = A1[j], B_val = selectB1[j], obs_all, 
                       obs0, B_index, quar_trunc) }))
    d_integral = rowSums( integrand*matrix(rep(mart_X, length(X1)), 
                                           byrow = TRUE, nrow = length(X1)) )
    d = A1*(1 - reversed_oddsratios_selectB1*mean_A0/(1 - mean_A0))*d_integral/mean_A0 
    
    e = (1 - A1)*(Ehat_givenA_selectB1 + (par_k - mean(Y_Ghat0[A0 == 1])))/(1 - mean_A0)
    
    return( round(m*(a - b - c + d - e), 3) )  
  },
    
  Stab_onestep_est = function(all_obs, chunk_size, elln, est_index, alpha, num_top = 1, quar_trunc){ 
    # chunk_size, elln, alpha are scalars.
    mt = rowMeans( sapply(0:(ceiling((n - elln)/chunk_size) - 1), function(i){ print(i)
                  # print(i)
                  old_obs = all_obs[1:(elln + i*chunk_size)]
                  if (i < ceiling((n - elln)/chunk_size) - 1) {
                     new_obs = all_obs[(elln + i*chunk_size + 1):(elln + (i + 1)*chunk_size)]
                  } else {
                     new_obs = all_obs[(elln + i*chunk_size + 1):n]  }
                  
                  Psi_pars = sapply(1:p, function(k){ 
                               Est_Psi_k(B_index = k, obs = old_obs, quar_trunc) })
                  sgn0 = rep(1, p)
                     
                  ### EXPENSIVE!!!
                  if (p <= 1e4) {
                    IF_star_mat = IF_star_self(m = sgn0, B_index = 1:p, obs = old_obs, quar_trunc)
                    utility = (Psi_pars + my_colMeans(IF_star_mat))
                  } else { utility = Psi_pars }   
                     
                  if (num_top < p) {
                    k0 = order(abs(utility), decreasing = TRUE)[seq(num_top)]
                  } else { k0 = seq(p) }
        
                  Psi_k0 = Psi_pars[k0]
                  sgn_k0 = 2*(utility[k0] >= 0) - 1
        
                  Bk0 = input_data$B[old_obs, k0, drop = FALSE]
        
                  mu_Bk0 = colMeans(Bk0)
                  sd_Bk0 = apply(Bk0, 2, sd)
            
                  if (est_index == 'subsamples') { obs_all_val = old_obs; obs0_val = old_obs
                  } else if (est_index == 'partial subsamples') { obs_all_val = all_obs; obs0_val = old_obs
                  } else if (est_index == 'whole samples') { obs_all_val = all_obs; obs0_val = all_obs }
      
         curr_sigma_inv0 = 1/sd(rowMeans(IF_star(m = sgn_k0, B_index = k0, obs_all = obs_all_val, 
                                                 obs0 = obs0_val, obs1 = old_obs, quar_trunc)))
         
         if (length(new_obs) > 1) {
           est0 = mean( (mean(sgn_k0*Psi_k0) 
                       + rowMeans(IF_star(m = sgn_k0, B_index = k0, obs_all = obs_all_val, 
                                          obs0 = obs0_val, obs1 = new_obs, quar_trunc)))
                       * curr_sigma_inv0 )
         } else {
           est0 = mean( (mean(sgn_k0*Psi_k0) 
                         + mean(IF_star(m = sgn_k0, B_index = k0, obs_all = obs_all_val, 
                                        obs0 = obs0_val, obs1 = new_obs, quar_trunc)))
                        * curr_sigma_inv0 )
         }  
         return( c(est0, curr_sigma_inv0) )  } ) )
       est = mt[1]/mt[2]
       se = 1/(mt[2] * sqrt(n - elln))
       ci = c( est - qnorm(1 - alpha/2)/(mt[2] * sqrt(n - elln)), 
               est + qnorm(1 - alpha/2)/(mt[2] * sqrt(n - elln)) )
       rej = 1*( ci[1] > 0 | ci[2] < 0 )
       p_val = 2*pnorm(abs(est/se), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
    
    return( list(rej = rej, Sn = abs(est/se), est = est, se = se, ci = ci, p_val = p_val) )  
  },
  
  Stab_onestep_kest = function(obs = 1:n, num_top = 1, quar_trunc){
    
    Psi_pars = sapply(1:p, function(k){ Est_Psi_k(B_index = k, obs, quar_trunc) })
    ### EXPENSIVE!!!
    if (p <= 1e4) {
      IF_star_mat = IF_star_self(m = rep(1, p), B_index = 1:p, obs, quar_trunc)
      utility = (Psi_pars + my_colMeans(IF_star_mat))
    } else { utility = Psi_pars }   
    
    if (num_top < p) {
      k_est = order(abs(utility), decreasing = TRUE)[seq(num_top)]
    } else { k_est = seq(p) }
    
    return( k_est )
  },  
  
  Oracle_onestep_est = function(alpha, quar_trunc, idx_taken = 1){
    # Using 1:2 here is to compile the vectorized syntax in codes, even if we only need the results of the first predictor 
    Psi_par0 = Est_Psi_k(B_index = idx_taken, obs = 1:n, quar_trunc)
    IF_star0 = IF_star_self(m = 1, B_index = idx_taken, obs = 1:n, quar_trunc)
    sigma_IF_star0 = sd(IF_star0)
    est = (Psi_par0 + mean(IF_star0))
    rej = 1*( sqrt(n)*abs(est/sigma_IF_star0) > qnorm(1 - alpha/2, 0, 1, lower.tail = TRUE, log.p = FALSE) )
    p_val = 2*pnorm( sqrt(n)*abs(est/sigma_IF_star0), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE )

    if (is.nan(est) | is.na(sigma_IF_star0)) {
      save(input_data, file = paste('NAdata_singleOSE_model', model, '_n', n, '_r', r, '.Rdata', sep = '') )
    }
    return( list( rej = rej, Sn = sqrt(n)*abs(est/sigma_IF_star0), est = est, se = sigma_IF_star0, 
                  p_val = p_val, k_est = 1 ) )
  },
  
  Naive_onestep_est = function(alpha, quar_trunc, num_top=1){
    
    Psi_pars = sapply(1:p, function(k){ 
                      Est_Psi_k(B_index = k, obs = 1:n, quar_trunc) });
    ### EXPENSIVE!!!
    if (p <= 1e4) {
      IF_star_mat0 = IF_star_self(m = rep(1, p), B_index = 1:p, obs = 1:n, quar_trunc)
      utility = (Psi_pars + my_colMeans(IF_star_mat0))
    } else { utility = Psi_pars }   
    
    if (num_top < p) {
      k0 = order(abs(utility), decreasing = TRUE)[seq(num_top)]
    } else { k0 = seq(p) }
    
    Psi_k0 = Psi_pars[k0]
    sgn_k0 = 2*(utility[k0] >= 0) - 1
    
    IF_star_mat = IF_star_self(m = sgn_k0, B_index = k0, obs = 1:n, quar_trunc)
    est = (Psi_k0 + mean(IF_star_mat))
    sigma_IF_star = sd(IF_star_mat)
    rej = 1*( sqrt(n)*abs(est/sigma_IF_star) 
              > qnorm(1 - alpha/2, 0, 1, lower.tail = TRUE, log.p = FALSE) )
    abs_test_stat_val = sqrt(n)*abs(est/sigma_IF_star)
    p_val = min(1, 2*pnorm(abs_test_stat_val, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
    
    return( list( rej = rej, Sn = abs_test_stat_val,
                  est = est, se = sigma_IF_star, p_val = p_val, k_est = k0 ) )
  },
  
  Bonf_onestep_ests = function(alpha, quar_trunc){
    
    res = do.call(rbind, lapply(1:p, function(k){
      est_Psi_k = Est_Psi_k(B_index = k, obs = 1:n, quar_trunc) 
      est_IF_star_func = IF_star_self(m = 1, B_index = k, obs = 1:n, quar_trunc)
      est = est_Psi_k + mean(est_IF_star_func)
      sigma_IF_star = sd(est_IF_star_func)
      res_k = data.frame('est' = est, 'sigma_IF_star' = sigma_IF_star)
      return( res_k ) }))
    
    test_stats = res[,'est']/res[,'sigma_IF_star']
    rej = 1*( sqrt(n)*max(abs(test_stats)) 
                            > qnorm(1 - alpha/(2*p), 0, 1, lower.tail = TRUE, log.p = FALSE) )
    abs_test_stat_vals = sqrt(n)*abs(test_stats)
    p_vals = 2*pnorm(abs_test_stat_vals, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
    p_val = min(1, p*2*pnorm(max(abs_test_stat_vals), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
    
    return( list( rej = rej, Sn = sqrt(n)*max(abs(test_stats)),
                  est = res[which.max(abs_test_stat_vals), 'est'], 
                  se = res[which.max(abs_test_stat_vals), 'sigma_IF_star'],
                  p_val = p_val, k_est = which.max(abs_test_stat_vals), 
                  k_ests = which(p_vals <= alpha/p) ) )
  }
 )
)