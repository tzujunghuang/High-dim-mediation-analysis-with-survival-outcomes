num_digits = 7
source('data_management.R')

subset_data = function(input_data, input_obs, input_Bindex){
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
}
    

KM_weight_func = function(input_data, obs){ # obs is a vector of observation/individual indices 
    n = length(input_data$X)
    data_km = data.frame(X = input_data$X[obs], delta = input_data$delta[obs])
    data_km = data_km[order(data_km$X),]
    prod = c(1, cumprod( sapply(1:(n - 1), function(i){ 
      ((n - i)/(n - i + 1))^(data_km$delta[i]) }) ) )
    kmwts = sapply(1:n, function(i){ data_km$delta[i]*prod[i]/(n - i + 1) } )
    return( kmwts )
}  
    

KM_SurF = function(input_data, t, obs){ # t is a scalar; obs is a vector of observation/individual indices 
    data_km = data.frame(X = input_data$X[obs], delta = input_data$delta[obs])
    km = survfit(Surv(X, 1-delta)~1, data = data_km)
    rm(data_km)
    
    survest = cbind(km$time, km$surv)
    if ( length(which(survest[, 1] <= t)) > 0 ) {
      return( survest[max(which(survest[, 1] <= t)), 2] )  
    } else { return( 1 ) }
}  

  
KM_SurF_self = function(input_data, obs, quar_trunc){ # obs is a vector of observation/individual indices 
    data_km = data.frame(X = input_data$X[obs], delta = input_data$delta[obs])
    tau = quantile(data_km$X, probs = quar_trunc)
    tau_surv = KM_SurF(input_data, tau, obs)
    
    km = survfit(Surv(X, 1-delta)~1, data = data_km)
    rm(data_km)
    
    km$surv[km$time > tau] = tau_surv
    survest = cbind(km$time, km$surv)
    return( survest )
}
  

Inverse_weight_func = function(input_data, x0, obs, quar_trunc, err_msg){
    # (x0, delta0) are the replication of subjects to be predicted
    # obs is a vector of observation/individual indices of x0
    n = length(input_data$X)
    tau = quantile(input_data$X[1:n], probs = quar_trunc)
    tau_surv = KM_SurF(input_data, tau, obs)
    
    res = sapply(1:length(x0), function(i){ 
      ifelse(x0[i] < tau, KM_SurF(input_data, x0[i], obs), tau_surv) })
    return( res )
}
  
  
Inverse_weight_func_self = function(input_data, x0, obs, quar_trunc, err_msg){
    # (x0, delta0) are the replication of subjects to be predicted
    # obs is a vector of observation/individual indices of x0
    n = length(input_data$X)
    tau = quantile(input_data$X[1:n], probs = quar_trunc)
    tau_surv = KM_SurF(input_data, tau, obs)
    
    if (length(x0) == 1) {
      inverse_weight = ifelse(x0 < tau, KM_SurF(input_data, x0, obs), tau_surv)
    } else {    
      if (length(unique(x0)) == length(x0)) { 
        KM_table = KM_SurF_self(input_data, obs, quar_trunc) 
      } else {
        KM_table0 = data.frame(KM_SurF_self(input_data, obs, quar_trunc))
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
}
 
 
Est2_Psi_k = function(input_data, B_index, obs, quar_trunc, block_size = 1e4){
    # B_index is a scalar or a vector of mediator indices; 
    # obs is a vector of observation/individual indices used for estimation
    num_portions = ceiling(length(B_index) / block_size) - 1
    
    if (num_portions == 0) {
      
      df = subset_data(input_data, input_obs = obs, input_Bindex = B_index)
      X = df$X; delta = df$delta; A = df$A; selectB = df$selectB
      rm(df)
      
      inverse_weight = Inverse_weight_func_self(input_data, x0 = X, obs, quar_trunc, 
                                                err_msg = 'Error in Est_Psi0_d Inverse_weight KM_table')
      Y_Ghat = X*delta / inverse_weight; Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0
      mean_A = mean(A)
      
      Ehat_givenAB = get_Ehat_givenAB(input_data, A_val = 1, obs0 = obs, obs1 = obs, B_index, quar_trunc)
      
      eta_10 = colMeans((1 - A)*Ehat_givenAB/(1 - mean_A))
      results = as.numeric(mean(Y_Ghat[A == 1]) - eta_10)
    } else {
      results = sapply(0:num_portions, function(i) {
        if (i < num_portions) {
          p0_vec = (1 + i*block_size):((i + 1)*block_size)
        } else {
          p0_vec = (1 + i*block_size):dim(input_data$B)[2]
        }
        
        df = subset_data(input_data, input_obs = obs, input_Bindex = p0_vec)
        X = df$X; delta = df$delta; A = df$A; selectB = df$selectB
        rm(df)
        
        inverse_weight = Inverse_weight_func_self(input_data, x0 = X, obs, quar_trunc, 
                                                  err_msg = 'Error in Est_Psi0_d Inverse_weight KM_table')
        Y_Ghat = X*delta / inverse_weight; Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0
        mean_A = mean(A)
        
        Ehat_givenAB = get_Ehat_givenAB(input_data, A_val = 1, obs0 = obs, obs1 = obs, 
                                        B_index = p0_vec, quar_trunc)
        
        eta_10 = colMeans((1 - A)*Ehat_givenAB/(1 - mean_A))
        return( as.numeric(mean(Y_Ghat[A == 1]) - eta_10) ) })
    }  
    return( results )
}


Est_Psi_k = function(input_data, B_index, obs, quar_trunc, ext_flag, 
                     block_size = 1e4){
  # B_index is a scalar or a vector of mediator indices; 
  # obs is a vector of observation/individual indices used for estimation
  X = input_data$X[obs]; delta = input_data$delta[obs]; A = input_data$A[obs]; 
  
  inverse_weight = Inverse_weight_func(input_data, x0 = X, obs, quar_trunc, 
                                       err_msg = 'Error in Est_Psi0_d Inverse_weight KM_table')
  Y_Ghat = X*delta / inverse_weight; Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0
  mean_A = mean(A)
  
  #Y_hat <- predict(lm(Y_Ghat ~ A))
  
  if (length(B_index) == 1) {
  
    selectB = input_data$B[obs, B_index]
    #Psi_k_1 = as.numeric(cov(selectB, Y_Ghat - Y_hat) / var(selectB))
    
    if (!ext_flag) {
      Psi_k_1 = as.numeric(coef(lm(Y_Ghat ~ A + selectB))[3])
    } else { 
      if ( dim(input_data$U)[2] == 1 ) {
        U = input_data$U[obs]
        Psi_k_1 = as.numeric(coef(lm(Y_Ghat ~ A + selectB + U))[3])
      } else {
        U = input_data$U[obs,]
        string_variables = paste(colnames(U), collapse = ' + ')
        string_formula = paste0("Y_Ghat ~ A + selectB + ", string_variables)
        
        fitdata = data.frame(Y_Ghat, A, selectB, U)
        Psi_k_1 = as.numeric(coef(lm(string_formula, data = fitdata))[3])
      } }
    Psi_k_2 = mean(selectB[A == 1]) - mean(selectB[A == 0])
    results = Psi_k_1*Psi_k_2

  } else {  
  
  num_portions = ceiling(length(B_index) / block_size) - 1
  
  if (num_portions == 0) {
    
    selectB = input_data$B[obs, B_index]
    #Psi_k_1 = as.numeric(colCovs(mat = selectB, y = Y_Ghat - Y_hat) / colVars(selectB, sd_use = F))
    
    if (!ext_flag) {
      Psi_k_1 = sapply(1:ncol(selectB), function(col) { 
        res = .lm.fit(as.matrix(cbind(A, selectB[,col])), Y_Ghat) 
        return( coefficients(res)[2] ) })
    } else { 
      
      if ( dim(input_data$U)[2] == 1 ) {
        U = input_data$U[obs]
      } else {
        U = input_data$U[obs,] }
        
      Psi_k_1 = sapply(1:ncol(selectB), function(col) { 
          res = .lm.fit(as.matrix(cbind(A, selectB[,col], U)), Y_Ghat) 
          return( coefficients(res)[2] ) }) }
    
    Psi_k_2 = colMeans(selectB[A == 1,]) - colMeans(selectB[A == 0,])
    results = Psi_k_1*Psi_k_2
      
  } else {
    results = sapply(0:num_portions, function(i) {
      if (i < num_portions) {
        p0_vec = (1 + i*block_size):((i + 1)*block_size)
      } else {
        p0_vec = (1 + i*block_size):dim(input_data$B)[2]
      }
      
      selectB = input_data$B[obs, p0_vec]
      #Psi_k_1 = as.numeric(colCovs(mat = selectB, y = Y_Ghat - Y_hat) / colVars(selectB, sd_use = F))
      if (!ext_flag) {
        Psi_k_1 = sapply(1:ncol(selectB), function(col) { 
          res = .lm.fit(as.matrix(cbind(A, selectB[,col])), Y_Ghat) 
          return( coefficients(res)[2] ) })
      } else { 
        
        if ( dim(input_data$U)[2] == 1 ) {
          U = input_data$U[obs]
        } else {
          U = input_data$U[obs,] }
        
        Psi_k_1 = sapply(1:ncol(selectB), function(col) { 
          res = .lm.fit(as.matrix(cbind(A, selectB[,col], U)), Y_Ghat) 
          return( coefficients(res)[2] ) })
      }
      Psi_k_2 = colMeans(selectB[A == 1,]) - colMeans(selectB[A == 0,])
      return( as.numeric(Psi_k_1*Psi_k_2) ) }) }
  }
  return( results )
}
 

get_Ehat_givenABX = function(input_data, A_val, obs_all, obs0, obs1, B_index, quar_trunc){
    # (A_val, B_val) are the given value of (A,B) this function to be evaluated at
    # obs0 is a vector of observation/individual indices used for estimation of the rest
    # obs1 is a vector of observation/individual indices to be predicted
    # B_index is a scalar or a vector of mediator indices
    X_all = input_data$X[obs_all]
    
    df0 = subset_data(input_data, input_obs = obs0, input_Bindex = B_index)
    X0 = df0$X; delta0 = df0$delta; A0 = df0$A; selectB0 = df0$selectB; rm(df0)
    
    df1 = subset_data(input_data, input_obs = obs1, input_Bindex = B_index)
    X1 = df1$X; A1 = df1$A; selectB1 = df1$selectB; rm(df1)
    
    selectB1 = matrix(selectB1, nrow = length(obs1), ncol = length(B_index))
    
    inverse_weight = Inverse_weight_func_self(input_data, x0 = X0, obs = obs0, quar_trunc, 
                                              err_msg = 'Error in IF_star_self Inverse_weight KM_table')
    
    time_comparison = outer(X0, X_all, '>=')
    n_s = colSums(time_comparison*(A0 == A_val)) 
    Y_Ghat_mat = time_comparison*(A0 == A_val)*((X0*delta0) / inverse_weight) 
    Y_Ghat_mat[is.nan(Y_Ghat_mat) | is.na(Y_Ghat_mat)] = 0
    mean_YGhat_mat = colMeans(Y_Ghat_mat)
    
    if ( is.matrix(selectB0) ) {
      
      Ehat_pointwise_vals = lapply(1:length(X_all), function(i){
        a_term = mean_YGhat_mat[i];
        Y_Ghat = ((X0*delta0) / inverse_weight)[A0 == A_val & X0 >= X_all[i]]
        Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0
        selectB = matrix(selectB0[A0 == A_val & X0 >= X_all[i],], ncol = length(B_index))
        
        n_s = nrow(selectB)
        if (n_s == 0) { val = matrix(0, byrow = T, ncol = ncol(selectB))
        } else { mean_selectB = colMeans(selectB)
        if (n_s == 1) { var_selectB = 1
        } else { var_selectB = apply(selectB, 2, var) }
        cov_Y_selectB = colMeans((selectB - matrix(rep(mean_selectB, n_s), 
                                                   byrow = T, nrow = n_s))*Y_Ghat)
        b_term = cov_Y_selectB / var_selectB; 
        b_term[is.nan(b_term) | is.na(b_term)] = 0
        b_term = matrix(b_term, byrow = T, ncol = ncol(selectB))
        c_term = do.call(cbind, lapply(1:ncol(selectB), function(col){ 
          selectB1[,col] - mean_selectB[col] }))
        
        val = a_term + do.call(cbind, lapply(1:ncol(selectB), function(col){ 
          b_term[,col]*c_term[,col] }))
        rm(selectB); rm(mean_selectB); rm(var_selectB); rm(cov_Y_selectB) }
        
        return( val ) } )
      
    } else {
      
      Ehat_pointwise_vals = do.call(cbind, lapply(1:length(X_all), function(i){
        a_term = mean_YGhat_mat[i];
        Y_Ghat = ((X0*delta0) / inverse_weight)[A0 == A_val & X0 >= X_all[i]]
        Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0
        selectB = selectB0[A0 == A_val & X0 >= X_all[i]]; 
        n_s = length(selectB)
        
        if (n_s == 0) { val = 0
        } else { mean_selectB = mean(selectB)
        
        if (length(selectB) > 1) { 
          var_selectB = var(selectB)
          cov_Y_selectB = mean((selectB - mean_selectB)*Y_Ghat)
          b_term = cov_Y_selectB / var_selectB; 
          b_term[is.nan(b_term) | is.na(b_term)] = 0
        } else { var_selectB = 'NA'; cov_Y_selectB = 0; b_term = 0 } 
        
        c_term = selectB1 - mean_selectB
        rm(selectB); rm(mean_selectB); rm(var_selectB); rm(cov_Y_selectB)
        
        val = a_term + b_term*c_term; val[is.nan(val) | is.na(val)] = 0 }
        
        return( val ) } ) )
    }
    return( Ehat_pointwise_vals )
}
  

get_Ehat_givenAB = function(input_data, A_val, obs0, obs1, B_index, quar_trunc){
    # A_val is the given value of A this function to be evaluated at
    # obs0 is a vector of observation/individual indices used for estimation of the rest
    # obs1 is a vector of observation/individual indices to be predicted
    # B_index is a scalar or a vector of mediator indices
    df0 = subset_data(input_data, input_obs = obs0, input_Bindex = B_index)
    X0 = df0$X; delta0 = df0$delta; A0 = df0$A; selectB0 = df0$selectB; rm(df0)
    
    df1 = subset_data(input_data, input_obs = obs1, input_Bindex = B_index)
    X1 = df1$X; A1 = df1$A; selectB1 = df1$selectB; rm(df1)
    
    selectB1 = matrix(selectB1, nrow = length(obs1), ncol = length(B_index))
    
    inverse_weight = Inverse_weight_func_self(input_data, x0 = X0, obs = obs0, quar_trunc, 
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
}
  

get_martingale = function(input_data, obs, obs1, quar_trunc){
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
    
    time_comparison3 = outer(X1, X, '>=')
    return( (X1 <= tau & delta1 == 0) - (X1 <= tau)*time_comparison3*hazard )
}
  

IF_star_self = function(input_data, m, B_index, obs, quar_trunc){
    # m is a scalar; B_index is a scalar or a vector of mediator indices; 
    # obs is the vector of old/whole individuals indices both for estimation and to be predicted
    n = nrow(input_data$B); p = ncol(input_data$B)
    df = subset_data(input_data, input_obs = obs, input_Bindex = B_index)
    X = df$X; delta = df$delta; A = df$A; selectB = df$selectB; rm(df)
    
    m_mat = matrix(rep(m, length(obs)), byrow = TRUE, ncol = length(B_index))
    
    #Psi_pars = unlist( Est_Psi_k(input_data, B_index, obs, quar_trunc, ext_flag) )
    
    tau = quantile(X, probs = quar_trunc)
    inverse_weight = Inverse_weight_func_self(input_data, x0 = X, obs, quar_trunc, 
                                              err_msg = 'Error in IF_star_self Inverse_weight KM_table')
    
    mean_A = mean(A)    
    Y_Ghat = (X*delta) / inverse_weight; Y_Ghat[is.nan(Y_Ghat) | is.na(Y_Ghat)] = 0
    
    Ehat_givenAB = get_Ehat_givenAB(input_data, A_val = 1, obs0 = obs, obs1 = obs, B_index, quar_trunc)
    mart_X = get_martingale(input_data, obs, obs1 = obs, quar_trunc)
    
    if ( is.matrix(selectB) ) {
      reversed_oddsratios = do.call(cbind, lapply(1:ncol(selectB), function(col){
        logit_model = glm(A ~ selectB[,col], family = binomial(link = 'logit'))
        coefs = as.numeric(coef(logit_model))
        
        reversed_oddsratio = sapply(selectB[,col], function(B){ 
          return( exp(-sum(coefs*c(1,B))) ) }) 
        return( reversed_oddsratio ) }))
      
      #density_ratios = do.call(cbind, lapply(1:ncol(selectB), function(col){
      #  logit_model = glm(selectB[,col] ~ A, family = binomial(link = 'logit'))
      #  coefs = as.numeric(coef(logit_model))
        
      #  density_ratio = sapply(selectB[,col], function(B){ 
      #    if (B == 0) { res = (1 + exp(-sum(coefs)))/(exp(-coefs[2]) + exp(-sum(coefs))) } 
      #    else if (B == 1) { res = (1 + exp(-sum(coefs)))/(1 + exp(-coefs[1])) } 
      #    return( res ) })
      #  return( density_ratio ) }))
      
      #density_ratios = do.call(cbind, lapply(1:ncol(selectB), function(col){
      #  mu0 = mean(selectB[A == 0, col]); sd0 = sd(selectB[A == 0, col])
      #  mu1 = mean(selectB[A == 1, col]); sd1 = sd(selectB[A == 1, col])  
      
      #  density_ratio = sapply(selectB[,col], function(B){ 
      #     dnorm(x = B, mean = mu0, sd = sd0)/dnorm(x = B, mean = mu1, sd = sd1) })
      #  return( density_ratio ) }))
      
      #Psi_pars_mat = matrix(rep(Psi_pars, length(obs)), byrow = TRUE, ncol = length(B_index))
      #eta_10_mat = mean(Y_Ghat[A == 1]) - Psi_pars_mat
      eta_10_mat = mean((1 - A)*Ehat_givenAB/(1 - mean_A))
      
      fk_a = -(1 - A)*(Ehat_givenAB - eta_10_mat)/(1 - mean_A) 
      fk_b = A*(Y_Ghat - mean(Y_Ghat[A == 1]))/mean_A
      fk_c = -A*reversed_oddsratios*(Y_Ghat - Ehat_givenAB)/(1 - mean_A)
      #fk_c = -A*density_ratios*(Y_Ghat - Ehat_givenAB)/mean_A
      fk = fk_a + fk_b + fk_c
      
      integrand = get_Ehat_givenABX(input_data, A_val = 1, obs_all = 1:n, obs0 = obs, obs1 = obs,
                                    B_index, quar_trunc)
      
      integral_mat = do.call(cbind, lapply(1:ncol(selectB), function(col){
        weighted_integrand = do.call(cbind, lapply(1:length(obs), function(i){
          return( integrand[[i]][,col]*mart_X[i,] ) }))
        return( rowMeans(weighted_integrand) ) }))
      
      fk_car = -A*integral_mat*(1 - reversed_oddsratios*mean_A/(1 - mean_A))/mean_A
      #fk_car = -A*integral_mat*(1 - density_ratios)/mean_A
      
      return( round( m_mat*(fk - fk_car), 3) ) 
      
    } else { 
      logit_model = glm(A ~ selectB, family = binomial(link = 'logit'))
      coefs = as.numeric(coef(logit_model))
      
      reversed_oddsratios = sapply(selectB, function(B){ 
        return( exp(-sum(coefs*c(1,B))) ) })
      
      #logit_model = glm(selectB ~ A, family = binomial(link = 'logit'))
      #coefs = as.numeric(coef(logit_model))
      
      #density_ratio = sapply(selectB, function(B){ 
      #  if (B == 0) { res = (1 + exp(-sum(coefs)))/(exp(-coefs[2]) + exp(-sum(coefs))) } 
      #  else if (B == 1) { res = (1 + exp(-sum(coefs)))/(1 + exp(-coefs[1])) } 
      #  return( res ) })
      
      #mu0 = mean(selectB[A == 0]); sd0 = sd(selectB[A == 0])
      #mu1 = mean(selectB[A == 1]); sd1 = sd(selectB[A == 1])
      
      #density_ratio = sapply(selectB, function(B){ 
      #   dnorm(x = B, mean = mu0, sd = sd0)/dnorm(x = B, mean = mu1, sd = sd1) })
      
      #eta_10 = mean(Y_Ghat[A == 1]) - Psi_pars
      eta_10 = mean((1 - A)*Ehat_givenAB/(1 - mean_A))
      
      fk_a = -(1 - A)*(Ehat_givenAB - eta_10)/(1 - mean_A) 
      fk_b = A*(Y_Ghat - mean(Y_Ghat[A == 1]))/mean_A
      fk_c = -A*reversed_oddsratios*(Y_Ghat - Ehat_givenAB)/(1 - mean_A)
      #fk_c = -A*density_ratio*(Y_Ghat - Ehat_givenAB)/mean_A
      fk = fk_a + fk_b + fk_c
      
      integrand = get_Ehat_givenABX(input_data, A_val = 1, obs_all = 1:n, obs0 = obs, obs1 = obs,
                                    B_index, quar_trunc)
      
      weighted_integrand = do.call(cbind, lapply(1:length(obs), function(i){
        return( integrand[i,]*mart_X[i,] ) }))
      integral = rowMeans(weighted_integrand)
      
      fk_car = -A*integral*(1 - reversed_oddsratios*mean_A/(1 - mean_A))/mean_A
      #fk_car = -A*integral*(1 - density_ratio)/mean_A
      
      return( round( m*(fk - fk_car), 3) )
    }  
}

  
IF_star = function(input_data, m, B_index, obs_all, obs0, obs1, quar_trunc, as.weights){
    # m is a scalar; B_index is a scalar or a vector of mediator indices; 
    # (obs_all, obs0, obs1) are the vectors of whole, old and new individual indices
    n = nrow(input_data$B); p = ncol(input_data$B)
    df = subset_data(input_data, input_obs = obs_all, input_Bindex = B_index)
    X = df$X; delta = df$delta; A = df$A; selectB = df$selectB; rm(df)
    
    df0 = subset_data(input_data, input_obs = obs0, input_Bindex = B_index)
    X0 = df0$X; delta0 = df0$delta; A0 = df0$A; selectB0 = df0$selectB; rm(df0)
    
    df1 = subset_data(input_data, input_obs = obs1, input_Bindex = B_index)
    X1 = df1$X; delta1 = df1$delta; A1 = df1$A; selectB1 = df1$selectB; rm(df1)
    
    inverse_weight0 = Inverse_weight_func_self(input_data, x0 = X0, obs = obs0, quar_trunc, 
                                               err_msg = 'Error in IF_star Inverse_weight KM_table0')
    inverse_weight1 = Inverse_weight_func(input_data, x0 = X1, obs = obs_all, quar_trunc, 
                                          err_msg = 'Error in IF_star Inverse_weight KM_table1')
    mean_A0 = mean(A0)
    Y_Ghat0 = (X0*delta0) / inverse_weight0; 
    Y_Ghat0[is.nan(Y_Ghat0) | is.na(Y_Ghat0)] = 0
    
    Y_Ghat1 = (X1*delta1) / inverse_weight1; 
    Y_Ghat1[is.nan(Y_Ghat1) | is.na(Y_Ghat1)] = 0
    
    Ehat_givenAB = get_Ehat_givenAB(input_data, A_val = 1, obs0, obs1, B_index, quar_trunc)
    
    ##### NEW estimates of reversed oddsratios at selectB1
    W0 = sapply(A0, function(val){ 
      return( (val == 1)*mean_A0 + (val == 0)*(1 - mean_A0) ) })
    
    logit_model0 = glm(A0 ~ selectB0, weight = round(1/W0), family = binomial(link = 'logit'))
    coefs = as.numeric(coef(logit_model0))
    
    reversed_oddsratios_selectB1 = sapply(selectB1, function(B){ 
      return( exp(-sum(coefs*c(1,B))) ) })
    
    #logit_model0 = glm(selectB0 ~ A0, family = binomial(link = 'logit'))
    #coefs = as.numeric(coef(logit_model0))
    
    #density_ratio = sapply(selectB1, function(B){ 
    #    if (B == 0) { res = (1 + exp(-sum(coefs)))/(exp(-coefs[2]) + exp(-sum(coefs))) } 
    #    else if (B == 1) { res = (1 + exp(-sum(coefs)))/(1 + exp(-coefs[1])) } 
    #    return( res ) })
    
    #mu0 = mean(selectB0[A0 == 0]); sd0 = sd(selectB0[A0 == 0])
    #mu1 = mean(selectB0[A0 == 1]); sd1 = sd(selectB0[A0 == 1])
    #density_ratio = sapply(selectB1, function(B){ 
    #  ( dnorm(x = B, mean = mu0, sd = sd0)/dnorm(x = B, mean = mu1, sd = sd1) ) })
    
    mart_X = get_martingale(input_data, obs = obs_all, obs1, quar_trunc)
    
    #par_k = Est_Psi_k(B_index, obs = obs0, quar_trunc, ext_flag)
    #eta_10 = mean(Y_Ghat0[A0 == 1]) - par_k
    eta_10 = mean((1 - A1)*(Ehat_givenAB)/(1 - mean_A0))
    
    fk_a = -(1 - A1)*(Ehat_givenAB - eta_10)/(1 - mean_A0) 
    
    if (as.weights) {
      fk_b = A1*(Y_Ghat1 - mean(Y_Ghat0[A0 == 1]))/mean_A0
    } else { fk_b = A1*(Y_Ghat1 - mean(Y_Ghat1[A1 == 1]))/mean_A0 }
    
    Ehat_givenAB1 = get_Ehat_givenAB(input_data, A_val = 1, obs0 = obs1, obs1, B_index, quar_trunc)
    #sgn_YGhat1 = 2*(Y_Ghat1 >= 0) - 1
    
    if (as.weights) {
      fk_c = -A1*reversed_oddsratios_selectB1*(Y_Ghat1 - Ehat_givenAB)/(1 - mean_A0)
      # fk_c = -A1*density_ratio*(Y_Ghat1 - Ehat_givenAB)/mean_A0
    } else { fk_c = -A1*reversed_oddsratios_selectB1*(Y_Ghat1 - Ehat_givenAB1)/(1 - mean_A0)
    # fk_c = -A1*density_ratio*(Y_Ghat1 - Ehat_givenAB1)/mean_A0
    }
    fk = fk_a + fk_b + fk_c
    temp0 = data.frame(A1, selectB1, Y_Ghat1, Ehat_givenAB, Ehat_givenAB1, 
                       reversed_oddsratios_selectB1, fk_a, fk_b, fk_c, fk)
    
    integrand = get_Ehat_givenABX(input_data, A_val = 1, obs_all = 1:n, obs0, obs1, B_index, quar_trunc)
    
    weighted_integrand = do.call(rbind, lapply(1:length(obs1), function(i){
      return( integrand[i,]*mart_X[i,] ) }))
    integral = rowMeans(weighted_integrand)
    
    fk_car = -A1*integral*(1 - reversed_oddsratios_selectB1*mean_A0/(1 - mean_A0))/mean_A0
    #fk_car = -A1*integral*(1 - density_ratio)/mean_A0
    
    return( round( m*(fk - fk_car), 3) )  
}
 
 
Stab_onestep_est = function(input_data, all_obs, chunk_size, elln, est_index, alpha, 
                            num_top = 1, quar_trunc, ext_flag, as.weights){ 
    n = nrow(input_data$B); p = ncol(input_data$B)
    # chunk_size, elln, alpha are scalars.
    mt = rowMeans( sapply(0:(ceiling((n - elln)/chunk_size) - 1), function(i){
      # print(i)
      old_obs = all_obs[1:(elln + i*chunk_size)]
      if (i < ceiling((n - elln)/chunk_size) - 1) {
        new_obs = all_obs[(elln + i*chunk_size + 1):(elln + (i + 1)*chunk_size)]
      } else {
        new_obs = all_obs[(elln + i*chunk_size + 1):n]  }
      
      Psi_pars = unlist( Est_Psi_k(input_data, B_index = 1:p, obs = old_obs, quar_trunc, ext_flag) )
      sgn0 = rep(1, p)
      
      ### EXPENSIVE!!!
      #if (p <= 5e1) {
      #  IF_star_mat = IF_star_self(input_data, m = sgn0, B_index = 1:p, obs = old_obs, 
      #                             quar_trunc)
      #  utility = (Psi_pars + colMeans(IF_star_mat))
      #} else { 
      utility = Psi_pars
      
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
      
      curr_sigma_inv0 = 1/sd(rowMeans(IF_star(input_data, m = sgn_k0, B_index = k0, obs_all = obs_all_val, 
                                              obs0 = obs0_val, obs1 = old_obs, quar_trunc,
                                              as.weights = TRUE)))
      
      if (length(new_obs) > 1) {
        est0 = (sgn_k0*Psi_k0 
                + mean(IF_star(input_data, m = sgn_k0, B_index = k0, obs_all = obs_all_val, 
                               obs0 = obs0_val, obs1 = new_obs, quar_trunc, as.weights)))*curr_sigma_inv0
        #est0 = mean( (mean(sgn_k0*Psi_k0) 
        #              + rowMeans(IF_star(m = sgn_k0, B_index = k0, obs_all = obs_all_val, 
        #                                 obs0 = obs0_val, obs1 = new_obs, quar_trunc,
        #                                 as.weights)))
        #             * curr_sigma_inv0 )
      } else {
        est0 = (sgn_k0*Psi_k0 
                + IF_star(input_data, m = sgn_k0, B_index = k0, obs_all = obs_all_val, 
                          obs0 = obs0_val, obs1 = new_obs, quar_trunc, as.weights))*curr_sigma_inv0
        #est0 = mean( (mean(sgn_k0*Psi_k0) 
        #              + mean(IF_star(m = sgn_k0, B_index = k0, obs_all = obs_all_val, 
        #                             obs0 = obs0_val, obs1 = new_obs, quar_trunc,
        #                             as.weights)))
        #             * curr_sigma_inv0 )
      }
      
      if (is.na(est0) | is.na(curr_sigma_inv0)) { return( c(0,0) )
      } else { return( c(est0, curr_sigma_inv0) ) }  } ) )
      #return( c(est0, curr_sigma_inv0) )  } ) )
    
    est = mt[1]/mt[2]
    se = 1/(mt[2] * sqrt(n - elln))
    ci = c( est - qnorm(1 - alpha/2)/(mt[2] * sqrt(n - elln)), 
            est + qnorm(1 - alpha/2)/(mt[2] * sqrt(n - elln)) )
    rej = 1*( ci[1] > 0 | ci[2] < 0 )
    p_val = 2*pnorm(abs(est/se), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
    
    return( list(rej = rej, Sn = abs(est/se), est = est, se = se, ci = ci, p_val = p_val) )  
}


Stab_onestep_kestseq = function(input_data, all_obs, chunk_size, elln,
                                num_top = 1, quar_trunc, ext_flag){
  n = nrow(input_data$B); p = ncol(input_data$B)
  # chunk_size, elln, alpha are scalars.
  kestseq = do.call(rbind, lapply(0:(ceiling((n - elln)/chunk_size) - 1), function(i){
    # print(i)
    old_obs = all_obs[1:(elln + i*chunk_size)]
    if (i < ceiling((n - elln)/chunk_size) - 1) {
      new_obs = all_obs[(elln + i*chunk_size + 1):(elln + (i + 1)*chunk_size)]
    } else {
      new_obs = all_obs[(elln + i*chunk_size + 1):n]  }
    
    Psi_pars = unlist( Est_Psi_k(input_data, B_index = 1:p, obs = old_obs, quar_trunc, ext_flag) )
    sgn0 = rep(1, p)
    
    ### EXPENSIVE!!!
    #if (p <= 5e1) {
    #  IF_star_mat = IF_star_self(input_data, m = sgn0, B_index = 1:p, obs = old_obs, 
    #                             quar_trunc)
    #  utility = (Psi_pars + colMeans(IF_star_mat))
    #} else { 
    utility = Psi_pars
    
    if (num_top < p) {
      k0 = order(abs(utility), decreasing = TRUE)[seq(num_top)]
    } else { k0 = seq(p) }
    
    return( data.frame('k0' = k0) ) }))
}
 

Stab_onestep_kest = function(input_data, obs, num_top = 1, quar_trunc, ext_flag){
    
    n = nrow(input_data$B); p = ncol(input_data$B)
    Psi_pars = unlist( Est_Psi_k(input_data, B_index = 1:p, obs, quar_trunc, ext_flag) )
    ### EXPENSIVE!!!
    #if (p <= 5e1) {
    #  IF_star_mat = IF_star_self(input_data, m = rep(1, p), B_index = 1:p, obs, 
    #                             quar_trunc)
    #  utility = (Psi_pars + colMeans(IF_star_mat))
    #} else { 
    utility = Psi_pars
    
    if (num_top < p) {
      k_est = order(abs(utility), decreasing = TRUE)[seq(num_top)]
    } else { k_est = seq(p) }
    
    return( k_est )
}  
  

Oracle_onestep_est = function(input_data, alpha, quar_trunc, ext_flag, idx_taken = 1){
    # Using 1:2 here is to compile the vectorized syntax in codes, even if we only need the results of the first predictor 
    n = nrow(input_data$B); p = ncol(input_data$B)
    Psi_par0 = Est_Psi_k(input_data, B_index = idx_taken, obs = 1:n, quar_trunc, ext_flag)
    sgn_Psi_par0 = 2*(Psi_par0 >= 0) - 1
    IF_star0 = IF_star_self(input_data, m = sgn_Psi_par0, B_index = idx_taken, obs = 1:n, 
                            quar_trunc)
    sigma_IF_star0 = sd(IF_star0)
    est = (Psi_par0 + mean(IF_star0))
    rej = 1*( sqrt(n)*abs(est/sigma_IF_star0) > qnorm(1 - alpha/2, 0, 1, lower.tail = TRUE, log.p = FALSE) )
    ci = c( est - qnorm(1 - alpha/2)*sigma_IF_star0/sqrt(n), 
            est + qnorm(1 - alpha/2)*sigma_IF_star0/sqrt(n) )
    p_val = 2*pnorm( sqrt(n)*abs(est/sigma_IF_star0), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE )

    if (is.nan(est) | is.na(sigma_IF_star0)) {
      save(input_data, file = paste('NAdata_singleOSE_model', model, '_n', n, '_r', r, '.Rdata', sep = '') )
    }
    return( list( rej = rej, Sn = sqrt(n)*abs(est/sigma_IF_star0), est = est, se = sigma_IF_star0, 
                  ci = ci, p_val = p_val, k_est = 1 ) )
}
  

Bonf_onestep_est = function(input_data, alpha, quar_trunc, ext_flag, num_top=1){
    
    n = nrow(input_data$B); p = ncol(input_data$B)
    Psi_pars = unlist( Est_Psi_k(input_data, B_index = 1:p, obs = 1:n, quar_trunc, ext_flag) );
    ### EXPENSIVE!!!
    #if (p <= 5e1) {
    #  IF_star_mat0 = IF_star_self(input_data, m = rep(1, p), B_index = 1:p, obs = 1:n, 
    #                              quar_trunc)
    #  utility = (Psi_pars + colMeans(IF_star_mat0))
    #} else { 
    utility = Psi_pars
    
    if (num_top < p) {
      k0 = order(abs(utility), decreasing = TRUE)[seq(num_top)]
    } else { k0 = seq(p) }
    
    Psi_k0 = Psi_pars[k0]
    sgn_k0 = 2*(utility[k0] >= 0) - 1
    
    IF_star_mat = IF_star_self(input_data, m = sgn_k0, B_index = k0, obs = 1:n, 
                               quar_trunc)
    est = (Psi_k0 + mean(IF_star_mat))
    sigma_IF_star = sd(IF_star_mat)
    rej = 1*( sqrt(n)*abs(est/sigma_IF_star) 
              > qnorm(1 - alpha/(2*p), 0, 1, lower.tail = TRUE, log.p = FALSE) )
    abs_test_stat_val = sqrt(n)*abs(est/sigma_IF_star)
    ci = c( est - qnorm(1 - alpha/(2*p))*sigma_IF_star/sqrt(n), 
            est + qnorm(1 - alpha/(2*p))*sigma_IF_star/sqrt(n) )
    p_val = min(1, p*2*pnorm(abs_test_stat_val, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
    
    return( list( rej = rej, Sn = abs_test_stat_val, est = est, 
                  se = sigma_IF_star, ci = ci, p_val = p_val, k_est = k0 ) )
}
  

Bonf_onestep_est2 = function(input_data, alpha, quar_trunc, ext_flag, block_size = 1e4){
    
    n = nrow(input_data$B); p = ncol(input_data$B)
    num_portions = ceiling(p / block_size) - 1
    
    if (num_portions == 0) {
      
      est_Psis = unlist( Est_Psi_k(input_data, B_index = 1:p, obs = 1:n, quar_trunc, ext_flag) )
      # sgn_est_Psis = 2*(est_Psis >= 0) - 1
      est_IF_star_mat = IF_star_self(input_data, m = rep(1,p), B_index = 1:p, 
                                     obs = 1:n, quar_trunc)
      est = est_Psis + colMeans(est_IF_star_mat)
      sigma_IF_star = c(colVars(est_IF_star_mat, sd_use = TRUE))
      test_stats = est/sigma_IF_star
      
      results = data.frame('est' = est, 'sigma_IF_star' = sigma_IF_star,
                           'abs_test_stats' = sqrt(n)*abs(test_stats),
                           'p_vals' = 2*pnorm(abs(test_stats), mean = 0, sd = 1, 
                                              lower.tail = FALSE, log.p = FALSE))
    } else {
      
      results = do.call(rbind, lapply(0:num_portions, function(i) {
        if (i < num_portions) {
          sub_mat = input_data$B[,(1 + i*block_size):((i + 1)*block_size)]
          p0_vec = (1 + i*block_size):((i + 1)*block_size)
        } else {
          sub_mat = input_data$B[,(1 + i*block_size):dim(input_data$B)[2]]  
          p0_vec = (1 + i*block_size):dim(input_data$B)[2]
        }
        est_Psis = unlist( Est_Psi_k(input_data, B_index = p0_vec, obs = 1:n, quar_trunc, ext_flag) )
        # sgn_est_Psis = 2*(est_Psis >= 0) - 1
        est_IF_star_mat = IF_star_self(input_data, m = rep(1,length(p0_vec)), B_index = p0_vec, 
                                       obs = 1:n, quar_trunc)
        est = est_Psis + colMeans(est_IF_star_mat)
        sigma_IF_star = colVars(est_IF_star_mat, sd_use = TRUE)
        test_stats = est/sigma_IF_star
        
        return( data.frame('est' = est, 'sigma_IF_star' = sigma_IF_star,
                           'abs_test_stats' = sqrt(n)*abs(test_stats),
                           'p_vals' = 2*pnorm(abs(test_stats), mean = 0, sd = 1, 
                                              lower.tail = FALSE, log.p = FALSE)) ) }))
    }
    p_vals = 2*pnorm(results$abs_test_stats, mean = 0, sd = 1, lower.tail = FALSE, 
                     log.p = FALSE)
    max_idx = which.max(results$abs_test_stats)
    max_abs_test_stats = max(results$abs_test_stats)  
    rej = 1*( max_abs_test_stats > qnorm(1 - alpha/(2*p), 0, 1, lower.tail = TRUE, log.p = FALSE) )
    p_val = min(1, p*2*pnorm(max_abs_test_stats, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
    
    return( list( rej = rej, Sn = max_abs_test_stats,
                  est = results$est[max_idx], se = results$sigma_IF_star[max_idx],
                  p_val = p_val, k_est = max_idx, 
                  k_ests = which(p_vals <= alpha/p) ) )
}
