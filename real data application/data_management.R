num_digits = 7

mostfrequent_percent = function(col){
  tbl = table(col)
  res = cbind(tbl, round(prop.table(tbl),2))
  colnames(res) = c('Count','Percentage')
  res[which.max(res[,2]), 2]
}

check_data_type = function(col){
  if (length(unique(col)) < 2) { data_type = 'void'
  } else if (length(unique(col)) == 2 & mostfrequent_percent(col) >= 0.95) {
    data_type = 'void binary'
  } else if (length(unique(col)) == 2 & mostfrequent_percent(col) < 0.95) {
    data_type = 'binary'
  } else if (length(unique(col)) > 2 & all((col == floor(col)) & (col >= 0))) {
    data_type = 'count/categorical'
  } else if (length(unique(col)) > 2 & all(col > floor(col))) {
    data_type = 'continuous'
  }
  return(data_type)
}

check_void = function(df, maineffect=TRUE){
  df_info = data.frame(names(df), apply(df, 2, function(x){check_data_type(x)}), 
                       apply(df, 2, function(x){mostfrequent_percent(x)}))
  
  if (maineffect) {
    names(df_info) = c('GeneticVariable', 'DataType', 'LargestProb_of_Success') # 65.58%
  } else{ names(df_info) = c('Interaction', 'DataType', 'LargestProb_of_Success') # 59.94% 
  }    
  
  df_info$Void = 1*grepl('void', as.vector(df_info$DataType))
  n_binary = dim(df_info[grepl('binary', df_info$DataType) | grepl('void', df_info$DataType),])[1]
  n_void = dim(df_info[grepl('void', df_info$DataType),])[1]
  return( list(df_info = df_info, void_percent = round(100*n_void/n_binary, 2)) ) 
}

colRound <- function(mat, num_digits=7, block_size=1e5) {
  
  num_portions = ceiling(dim(mat)[2] / block_size) - 1
      
  if (num_portions == 0) {
      results = round(mat, num_digits)
  } else {
      results = do.call(cbind, lapply(0:num_portions, function(i) {
                   if (i < num_portions) {
                       sub_mat = mat[,(1 + i*block_size):((i + 1)*block_size)]
                   } else {
                       sub_mat = mat[,(1 + i*block_size):dim(mat)[2]]  }
          
                   round(sub_mat, num_digits) })) }
  return(results)
} 

colSubtract <- function(mat, colMeans_mat, block_size=1e4) {
  
  num_portions = ceiling(dim(mat)[2] / block_size) - 1
  
  if (num_portions == 0) {
    n0 = dim(mat)[1];
    results = (mat - matrix(rep(colMeans_mat, n0), byrow = TRUE, nrow = n0))
  } else {
    results <- do.call(cbind, lapply(0:num_portions, function(i) {
      print(i)
      if (i < num_portions) {
        n0 = dim(mat)[1];
        sub_mat = mat[,(1 + i*block_size):((i + 1)*block_size)]
        sub_mean_mat = matrix(rep(colMeans_mat[(1 + i*block_size):((i + 1)*block_size)], n0), 
                              byrow = TRUE, nrow = n0)
        result = sub_mat - sub_mean_mat
      } else {
        n0 = dim(mat)[1];
        sub_mat = mat[,(1 + i*block_size):dim(mat)[2]]
        sub_mean_mat = matrix(rep(colMeans_mat[(1 + i*block_size):dim(mat)[2]], n0), byrow = TRUE, nrow = n0)
        result = sub_mat - sub_mean_mat }
      return( result ) })) }
  return( results )
}  

colIF_ipw <- function(mat, colMeans_mat, colVars_mat, colCovs_mat, m, 
                      central_Y_Ghat = Y_Ghat - mean_YGhat, block_size = 1e4) {
  
  num_portions = ceiling(dim(mat)[2] / block_size) - 1
  
  if (num_portions == 0) {
    n0 = dim(mat)[1] 
    mu_selectU_0 = matrix(rep(colMeans_mat, n0), byrow = TRUE, nrow = n0)
    var_selectU_0 = matrix(rep(colVars_mat, n0), byrow = TRUE, nrow = n0)
    cov_selectU0_Y = matrix(rep(colCovs_mat, n0), byrow = TRUE, nrow = n0)
    
    a = (mat - mu_selectU_0)*central_Y_Ghat / var_selectU_0
    b = ((mat - mu_selectU_0)^2)*cov_selectU0_Y / (var_selectU_0)^2
    results = matrix(rep(m,n0), byrow = TRUE, nrow = n0)*(a - b)
    
  } else {
    results <- do.call(cbind, lapply(0:num_portions, function(i) {
      if (i < num_portions) {
        n0 = dim(mat)[1] 
        sub_mat = mat[,(1 + i*block_size):((i + 1)*block_size)]
        sub_mu_selectU_0 = matrix(rep(colMeans_mat[(1 + i*block_size):((i + 1)*block_size)], n0), 
                                  byrow = TRUE, nrow = n0)
        sub_var_selectU_0 = matrix(rep(colVars_mat[(1 + i*block_size):((i + 1)*block_size)], n0), 
                                   byrow = TRUE, nrow = n0)
        sub_cov_selectU0_Y = matrix(rep(colCovs_mat[(1 + i*block_size):((i + 1)*block_size)], n0), 
                                    byrow = TRUE, nrow = n0)
        
        a = (sub_mat - sub_mu_selectU_0)*central_Y_Ghat / sub_var_selectU_0
        b = ((sub_mat - sub_mu_selectU_0)^2)*sub_cov_selectU0_Y / (sub_var_selectU_0)^2
        result = matrix(rep(m[(1 + i*block_size):((i + 1)*block_size)], n0), 
                        byrow = TRUE, nrow = n0)*(a - b)
        
      } else {
        n0 = dim(mat)[1]
        sub_mat = mat[,(1 + i*block_size):dim(mat)[2]]  
        sub_mu_selectU_0 = matrix(rep(colMeans_mat[(1 + i*block_size):dim(mat)[2]], n0), 
                                  byrow = TRUE, nrow = n0)
        sub_var_selectU_0 = matrix(rep(colVars_mat[(1 + i*block_size):dim(mat)[2]], n0), 
                                   byrow = TRUE, nrow = n0)
        sub_cov_selectU0_Y = matrix(rep(colCovs_mat[(1 + i*block_size):dim(mat)[2]], n0), 
                                    byrow = TRUE, nrow = n0)
        
        a = (sub_mat - sub_mu_selectU_0)*central_Y_Ghat / sub_var_selectU_0
        b = ((sub_mat - sub_mu_selectU_0)^2)*sub_cov_selectU0_Y / (sub_var_selectU_0)^2
        result = matrix(rep(m[(1 + i*block_size):dim(mat)[2]], n0), byrow = TRUE, nrow = n0)*(a - b) }
      
      return( result ) })) }
  
  return( results )
} 

colIF_CAR <- function(mat, colMeans_mat, colVars_mat, mart_val, m, block_size=1e4) {
  
  num_portions = ceiling(dim(mat)[2] / block_size) - 1
  
  if (num_portions == 0) {
    n0 = dim(mat)[1] 
    mu_selectU_0 = matrix(rep(colMeans_mat, n0), byrow = TRUE, nrow = n0)
    var_selectU_0 = matrix(rep(colVars_mat, n0), byrow = TRUE, nrow = n0)
    results = ( matrix(rep(m,n0), byrow = TRUE, nrow = n0)
                * (mat - mu_selectU_0)*mart_val / var_selectU_0 )
    
  } else {
    results <- do.call(cbind, lapply(0:num_portions, function(i) {
      if (i < num_portions) {
        n0 = dim(mat)[1] 
        sub_mat = mat[,(1 + i*block_size):((i + 1)*block_size)]
        sub_mu_selectU_0 = matrix(rep(colMeans_mat[(1 + i*block_size):((i + 1)*block_size)], n0), 
                                  byrow = TRUE, nrow = n0)
        sub_var_selectU_0 = matrix(rep(colVars_mat[(1 + i*block_size):((i + 1)*block_size)], n0), 
                                   byrow = TRUE, nrow = n0)
        result = ( matrix(rep(m[(1 + i*block_size):((i + 1)*block_size)], n0), 
                          byrow = TRUE, nrow = n0)
                    * (sub_mat - sub_mu_selectU_0)*mart_val / sub_var_selectU_0 )
        
      } else {
        n0 = dim(mat)[1]
        sub_mat = mat[,(1 + i*block_size):dim(mat)[2]]  
        sub_mu_selectU_0 = matrix(rep(colMeans_mat[(1 + i*block_size):dim(mat)[2]], n0), byrow = TRUE, nrow = n0)
        sub_var_selectU_0 = matrix(rep(colVars_mat[(1 + i*block_size):dim(mat)[2]], n0), byrow = TRUE, nrow = n0)
        result = ( matrix(rep(m[(1 + i*block_size):dim(mat)[2]], n0), byrow = TRUE, nrow = n0)
                   * (sub_mat - sub_mu_selectU_0)*mart_val / sub_var_selectU_0 ) }
      return( result ) })) }
  
  return( results )
} 

my_colMeans <- function(mat, block_size=1e5) {
  
  num_portions = ceiling(dim(mat)[2] / block_size) - 1
  
  if (num_portions == 0) {
    colmean = t(colMeans(mat, na.rm = TRUE))
  } else {
    colmean = do.call(cbind, lapply(0:num_portions, function(i) {
      if (i < num_portions) {
        sub_mat = mat[,(1 + i*block_size):((i + 1)*block_size)]
      } else {
        sub_mat = mat[,(1 + i*block_size):dim(mat)[2]]  }
      
      t(colMeans(sub_mat, na.rm = TRUE)) }))  }
  return(colmean)
}

colStandardization <- function(mat, colSD, block_size=1e4) {
  
  num_portions = ceiling(dim(mat)[2] / block_size) - 1
  
  if (num_portions == 0) {
    n0 = dim(mat)[1]
    results = as.matrix((mat - matrix(rep(my_colMeans(mat), n0), nrow = n0, byrow = TRUE))/matrix(rep(colSD, n0), 
                                                                                                  nrow = n0, byrow = TRUE))
  } else {
    results <- do.call(cbind, lapply(0:num_portions, function(i) {
                  if (i < num_portions) {
                      sub_mat = mat[,(1 + i*block_size):((i + 1)*block_size)]
                      sub_colSD = colSD[(1 + i*block_size):((i + 1)*block_size)]
                  } else {
                      sub_mat = mat[,(1 + i*block_size):dim(mat)[2]]  
                      sub_colSD = colSD[(1 + i*block_size):dim(mat)[2]] }
        
                  nn = dim(sub_mat)[1]
                  result = (sub_mat - matrix(rep(my_colMeans(sub_mat), nn), nrow = nn, byrow = TRUE))/matrix(rep(sub_colSD, nn), 
                                                                                                           nrow = nn, byrow = TRUE)
                  as.matrix(result) })) }
  return(results)
}  

colVars <- function(mat, sd_use, block_size=1e4) {
  
  if (is.null(dim(mat))) {
    results = t(rep(0, length(mat))) 
  } else {
    if (dim(mat)[1] <= 1) {
      results = t(rep(0, length(mat)))
    } else {
      num_portions = ceiling(dim(mat)[2] / block_size) - 1
      
      if (num_portions == 0) {
        nn = colSums(!is.na(mat))
        colVar = (colMeans((mat)^2, na.rm = TRUE) - (colMeans(mat, na.rm = TRUE))^2)*nn/(nn - 1)
        if (sd_use) { results = t(sqrt(colVar)) } else { results = t(colVar) }
        
      } else {
        results = do.call(cbind, lapply(0:num_portions, function(i) {
                      if (i < num_portions) {
                          sub_mat = mat[,(1 + i*block_size):((i + 1)*block_size)]
                      } else {
                          sub_mat = mat[,(1 + i*block_size):dim(mat)[2]]  }
          
                      nn = colSums(!is.na(sub_mat))
                      colVar = (colMeans((sub_mat)^2, na.rm = TRUE) - (colMeans(sub_mat, na.rm = TRUE))^2)*nn/(nn - 1)
          
                      if (sd_use) { result = sqrt(colVar) } else { result = colVar }
                      t(result) })) } 
      }
    }
  return(results)
} 

colCovs <- function(mat, y, block_size=1e5, na.rm=TRUE) {
  
  num_portions = ceiling(dim(mat)[2] / block_size) - 1
  
  if (num_portions == 0) {
      colCov = t(cov(mat, y))
  } else {
      colCov = do.call(cbind, lapply(0:num_portions, function(i) {
                  if (i < num_portions) {
                      sub_mat = mat[,(1 + i*block_size):((i + 1)*block_size)]
                  } else {
                      sub_mat = mat[,(1 + i*block_size):dim(mat)[2]]  }
          
                  t(cov(sub_mat, y)) }) ) }
  return(colCov)
} 

