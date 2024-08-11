library(dplyr)
library(tibble)
library(AHM)
library(mvtnorm)


tmat=function(X,f,k=max(f)){
  ###tapply but to columns of matrix
  ## X nxp matrix, including intercept if appropriate
  ## f a factor variable of length n with k classes, stored as an integer
  N = length(f)
  if(length(X)==N){
    X=matrix(X,ncol=1)
  }
  #outmat <- matrix(0,k,ncol(X))
  d=as_tibble(cbind(f,X))
  psum = d %>%
    group_by(f) %>%
    summarize_all(sum)
  outmat<-data.matrix(psum[,-1])
  outmat
}

sum_z_z_T = function(Z,f){
  k = ncol(Z)
  Z_Z_t= NULL
  for(i in 1:k){
    for(j in 1:k){
      Z_Z_t= cbind(Z_Z_t, Z[,i]*Z[,j])
    }
  }
  return(tmat(Z_Z_t,f))
}

factor_analysis_with_just_client_intercept = function(y,f1,f2,X = NULL, initial_v = NULL, initial_beta = NULL, initial_sigma2e = NULL,  max_iter = 200, conv_thres = 10^(-5), test_data = FALSE, n_test = NULL , X_test = NULL, y_test = NULL, f1_test = NULL, f2_test = NULL){
  start.time = Sys.time()
  ngl = c()
  prediction_error = c()
  training_error = c()
  Changes_in_beta = c()
  Changes_in_sigma2e = c()
  Changes_in_v = c()
  Changes_in_f_v = c()
  
  N=length(f1)
  R = length(unique(f1))
  C = length(unique(f2))
  n_is=tmat(rep(1,N),f1)
  m_js=tmat(rep(1,N),f2)
  #print(initial_beta==NULL)
  ## initialise V to be from Standard Normal distribution and covariance matrix Psi from chi square distribution, 
  ## beta to be ols.
  if(!is.null(X)){
    q = ncol(X)
    if(!is.null(initial_beta)){
      old_beta = initial_beta 
    }else{
      fit0 = lm(y~X-1)
      old_beta=fit0$coef
    }
    f_x = X%*% old_beta
    Changes_in_f_x = c()
  }else{
    f_x = 0
  }
  
 
  if(!is.null(initial_v)){
    old_v = initial_v
  }else{
    old_v=rnorm(1)
  }
  if(!is.null(initial_sigma2e)){
    old_sigma2e = initial_sigma2e   
  }else{
    old_sigma2e = 1
  }
  change = 1
  iter = 1
  f_v = 0
  
  while(iter<=max_iter && change>=conv_thres){
    if(!(is.null(X))){
      r_tilde = y-X%*%old_beta
    }else{
      r_tilde = y
    }
    
        
    r_tilde_sum = tmat(r_tilde, f1)
    var_u_i = (1 + n_is*old_v^2/old_sigma2e)^(-1)
    E_u_i = var_u_i * r_tilde_sum*old_v/old_sigma2e
    E_u2_i = E_u_i^2 + var_u_i
    sum_r_tilde_square = sum(r_tilde^2)
    negative_log_likelihood = N*log(old_sigma2e) + sum_r_tilde_square/as.numeric(old_sigma2e) + sum(log(1 + n_is*old_v^2/ old_sigma2e)) - sum(((old_v*r_tilde_sum)^2)/(old_sigma2e*(n_is*old_v^2 + old_sigma2e)))
 
    
    new_v = (sum(E_u_i[f1]*r_tilde))/sum(E_u2_i[f1])

    new_sigma2e = (sum_r_tilde_square - sum(E_u2_i[f1])*new_v^2)/N
    
    ## beta update
    new_f_v = E_u_i[f1]*new_v
    change_fv = sqrt(sum((new_f_v - f_v)^2))/sqrt(sum((f_v)^2))
    Changes_in_f_v = c(Changes_in_f_v, change_fv)
    f_v = new_f_v
    r_tilde = y - f_v
    change_sigma2e = abs(new_sigma2e/old_sigma2e - 1)
    
    
    Changes_in_sigma2e = c(Changes_in_sigma2e, change_sigma2e)
    change_v = abs(new_v/old_v - 1)
    
    change = max(change_fv, change_v, change_sigma2e)
    old_sigma2e = as.numeric(new_sigma2e)
    
    old_v = new_v

    ngl = c(ngl, negative_log_likelihood)
    if(!is.null(X)){
      new_beta = lm(r_tilde~X-1)$coef
      new_f_x = X%*% new_beta
      change_f_x = sqrt(sum((new_f_x - f_x)^2))/sqrt(sum((f_x)^2))
      Changes_in_f_x = c(Changes_in_f_x, change_f_x)
      change_beta = sqrt(sum((new_beta - old_beta)^2))/sqrt(sum(old_beta^2))
      change = max(change, change_beta,change_f_x)
      Changes_in_beta = c(Changes_in_beta, change_beta)
      old_beta = new_beta
      f_x = X%*% old_beta
      train_error = mean((y - f_v - f_x)^2)
      training_error = c(training_error, train_error)
      residuals = y_test - X_test %*% old_beta
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
          }else{
            U = E_u_i[f1_test[i]]
          }

          residuals[i] = residuals[i] - new_v * U
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng beta = ", change_beta,  "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng beta = ", change_beta, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "\n") 
      }
    }else{
      train_error = mean((y-f_v)^2)
      training_error = c(training_error, train_error)
      residuals = y_test
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
          }else{
            U = E_u_i[f1_test[i]]
          }
          
          residuals[i] = residuals[i] - U*new_v
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "\n") 
      }
    }
    iter = iter + 1    
  }
  end.time = Sys.time()
  time_taken = end.time - start.time
  if(is.null(X)){
    return(list(time_taken = time_taken,  v = new_v, sigma2e = new_sigma2e, U = E_u_i, ngl = ngl, training_error = training_error, prediction_error = prediction_error, Changes_in_f_v = Changes_in_f_v, Changes_in_sigma2e = Changes_in_sigma2e))
    
  }else{
    return(list(time_taken = time_taken, beta = new_beta, v = new_v, sigma2e = new_sigma2e, U = E_u_i, ngl = ngl, training_error = training_error, prediction_error = prediction_error, Changes_in_beta = Changes_in_beta,  Changes_in_f_v = Changes_in_f_v, Changes_in_f_x = Changes_in_f_x, Changes_in_sigma2e = Changes_in_sigma2e))
  }
}


factor_analysis_with_just_client_item_intercepts = function(y,f1,f2,X = NULL, initial_beta = NULL, initial_A = NULL,initial_B = NULL, initial_sigma2e = NULL, initial_sigma2a = NULL, initial_sigma2b = NULL, club_a = FALSE, club_b = FALSE, max_iter = 200, conv_thres = 10^(-5), test_data = FALSE, n_test = NULL , X_test = NULL, y_test = NULL, f1_test = NULL, f2_test = NULL){
  start.time = Sys.time()
  prediction_error = c()
  training_error = c()
  Changes_in_beta = c()
  Changes_in_sigma2e = c()
  Changes_in_sigma2a = c()
  Changes_in_sigma2b = c()
  Changes_in_f_a = c()
  Changes_in_f_b = c()
  
  
  N=length(f1)
  q = ncol(X)
  R = length(unique(f1))
  C = length(unique(f2))
  n_is=tmat(rep(1,N),f1)
  m_js=tmat(rep(1,N),f2)
  #print(initial_beta==NULL)
  ## initialise V to be from Standard Normal distribution and covariance matrix Psi from chi square distribution, 
  ## beta to be ols.
  if(!is.null(X)){
    if(!(is.null(initial_beta))){
      old_beta = initial_beta 
      f_x = X%*%old_beta
    }else{
      fit0 = lm(y~X-1)
      old_beta=fit0$coef
      f_x = X%*%old_beta
      new_f_x = f_x
      Changes_in_f_x = c()
    } 
  }else{
    f_x = 0
    new_f_x = f_x
  }
  
  
  if(!is.null(initial_sigma2e)){
    old_sigma2e = as.numeric(initial_sigma2e)   
  }else{
    old_sigma2e = 1
  }
  if(!(is.null(initial_sigma2a))){
    old_sigma2a = as.numeric(initial_sigma2a)   
  }else{
    old_sigma2a = 1
  }
  if(!is.null(initial_A)){
    old_A = initial_A  
  }else{
    old_A = rnorm(R)
  }
  if(!is.null(initial_sigma2b)){
    old_sigma2b = as.numeric(initial_sigma2b)  
  }else{
    old_sigma2b = 1
  }
  if(!is.null(initial_B)){
    old_B = initial_B  
  }else{
    old_B = rnorm(C)
  }
  change = 1
  iter = 1
  f_a = old_A[f1]
  f_b = old_B[f2]
  
  while(iter<=max_iter && change>=conv_thres){
    
    var_a_given_y = 1/(n_is/(old_sigma2e) + 1/ (old_sigma2a))
    var_b_given_y = 1/(m_js/(old_sigma2e) + 1/ (old_sigma2b))
    
    
    

    ## a update
    
    r_a = y - new_f_x - f_b
    
    E_a = (tmat(r_a,f1))/(n_is + (old_sigma2e)/(old_sigma2a))
    
    new_f_a = E_a[f1]
    
    new_sigma2a = mean((E_a)^2) + mean(var_a_given_y)
    

    ## b update
    
    r_b = y - new_f_x - new_f_a
    
    E_b = (tmat(r_b,f2))/(m_js + (old_sigma2e)/(old_sigma2b))
    
    new_f_b = E_b[f2]
    
    new_sigma2b = mean((E_b)^2) + mean(var_b_given_y)  
    
    r_tilde = y-new_f_x - new_f_a - new_f_b
    
    new_sigma2e = (sum(var_a_given_y[f1] + var_b_given_y[f2]) + sum(r_tilde^2))/N
    

    change_fa = sqrt(sum((new_f_a - f_a)^2))/sqrt(sum(f_a^2))
    Changes_in_f_a = c(Changes_in_f_a, change_fa)
    change_fb = sqrt(sum((new_f_b - f_b)^2))/sqrt(sum(f_b^2))
    Changes_in_f_b = c(Changes_in_f_b, change_fb)
    f_a = new_f_a
    f_b = new_f_b
    change_sigma2e = abs(new_sigma2e/old_sigma2e - 1)
    change_sigma2a = abs(new_sigma2a/old_sigma2a - 1)
    change_sigma2b = abs(new_sigma2b/old_sigma2b - 1)
    change = max(change_fa,change_sigma2e, change_sigma2a, change_fb, change_sigma2b)
    
    Changes_in_sigma2e = c(Changes_in_sigma2e, change_sigma2e)
    Changes_in_sigma2a = c(Changes_in_sigma2a, change_sigma2a)
    Changes_in_sigma2b = c(Changes_in_sigma2b, change_sigma2b)
    
    old_sigma2e = as.numeric(new_sigma2e)
    old_sigma2a = as.numeric(new_sigma2a)
    old_sigma2b = as.numeric(new_sigma2b)
    
    if(!is.null(X)){
      r_tilde = y - new_f_a - new_f_b
      new_beta = lm(r_tilde~X-1)$coef
      change_beta = sqrt(sum((new_beta - old_beta)^2))/sqrt(sum(old_beta^2))
      Changes_in_beta = c(Changes_in_beta, change_beta)
      old_beta = new_beta
      new_f_x = X%*% new_beta
      change_f_x = sqrt(sum((new_f_x - f_x)^2))/sqrt(sum(f_x^2))
      f_x = new_f_x
      Changes_in_f_x = c(Changes_in_f_x, change_f_x)
      change = max(change, change_beta,change_f_x)
      train_error = mean((y-X%*%old_beta - E_a[f1] - E_b[f2])^2)
      training_error = c(training_error, train_error)
      if(test_data){
        residuals = y_test - X_test %*% old_beta
        for(i in 1:n_test){
          if(f1_test[i]==0){
            a = 0
          }else{
            a = E_a[f1_test[i]]
          }
          if(f2_test[i]==0){
            b = 0
          }else{
            b = E_b[f2_test[i]]
          }
          residuals[i] = residuals[i] - a - b
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng beta = ", change_beta, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng sigma2a", change_sigma2a, "chng fa", change_fa,"chng sigma2b", change_sigma2b, "chng fb", change_fb,  "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng beta = ", change_beta, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng sigma2a", change_sigma2a, "chng fa", change_fa,"chng sigma2b", change_sigma2b, "chng fb", change_fb, "training_error = ",train_error, "\n") 
      }
    }else{
      train_error = mean((y - E_a[f1] - E_b[f2] )^2)
      training_error = c(training_error, train_error)
      if(test_data){
        residuals = y_test 
        for(i in 1:n_test){
          if(f1_test[i]==0){
            a = 0
          }else{
            a = E_a[f1_test[i]]
          }
          if(f2_test[i]==0){
            b = 0
          }else{
            b = E_b[f2_test[i]]
          }
          residuals[i] = residuals[i] - a - b
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng sigma2e",change_sigma2e, "chng sigma2a", change_sigma2a, "chng fa", change_fa, "chng sigma2b", change_sigma2b, "chng fb", change_fb,  "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng sigma2e",change_sigma2e, "chng sigma2a", change_sigma2a, "chng fa", change_fa, "chng sigma2b", change_sigma2b, "chng fb", change_fb,  "training_error = ",train_error, "\n") 
      }
    }
    iter = iter + 1    
  }
  end.time = Sys.time()
  time_taken = end.time - start.time
  if(is.null(X)){
    return(list(time_taken = time_taken, A = E_a, B = E_b, sigma2e = new_sigma2e,sigma2a = new_sigma2a, Changes_in_f_a = Changes_in_f_a, Changes_in_sigma2a = Changes_in_sigma2a, Changes_in_f_b = Changes_in_f_b, Changes_in_sigma2b = Changes_in_sigma2b, training_error = training_error, prediction_error = prediction_error, Changes_in_sigma2e = Changes_in_sigma2e))
    
  }else{
    return(list(time_taken = time_taken, Changes_in_f_a = Changes_in_f_a, A = E_a, B = E_b,  Changes_in_sigma2a = Changes_in_sigma2a,Changes_in_f_b = Changes_in_f_b, Changes_in_sigma2b = Changes_in_sigma2b,beta = new_beta, sigma2e = new_sigma2e, training_error = training_error, prediction_error = prediction_error, Changes_in_beta = Changes_in_beta,Changes_in_f_x = Changes_in_f_x, Changes_in_sigma2e = Changes_in_sigma2e))
  }
}



factor_analysis_without_intercepts = function(y,f1,f2,K,X = NULL, initial_V = NULL, initial_beta = NULL, initial_sigma2e = NULL, club_V = FALSE, club_U = FALSE, max_iter = 200, conv_thres = 10^(-5), test_data = FALSE, n_test = NULL , X_test = NULL, y_test = NULL, f1_test = NULL, f2_test = NULL){
  start.time = Sys.time()
  ngl = c()
  prediction_error = c()
  training_error = c()
  Changes_in_beta = c()
  Changes_in_sigma2e = c()
  Changes_in_V = c()
  Changes_in_f_v = c()
  
  N=length(f1)
  R = length(unique(f1))
  C = length(unique(f2))
  n_is=tmat(rep(1,N),f1)
  m_js=tmat(rep(1,N),f2)
  #print(initial_beta==NULL)
  ## initialise V to be from Standard Normal distribution and covariance matrix Psi from chi square distribution, 
  ## beta to be ols.
  if(!is.null(X)){
    q = ncol(X)
    if(!is.null(initial_beta)){
      old_beta = initial_beta 
    }else{
      fit0 = lm(y~X-1)
      old_beta=fit0$coef
    }
    f_x = X%*% old_beta
    Changes_in_f_x = c()
  }else{
    f_x = 0
  }
  
  if(!is.null(initial_V)){
    old_V = initial_V
  }else{
    old_V=matrix(rnorm(C*K),C,K)
  }
  if(!is.null(initial_sigma2e)){
    old_sigma2e = initial_sigma2e   
  }else{
    old_sigma2e = 1
  }
  change = 1
  iter = 1
  f_v = 0
  
  while(iter<=max_iter && change>=conv_thres){
    sum_vvT =  sum_z_z_T(old_V[f2,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f1)
    if(!(is.null(X))){
      r_tilde = y-X%*%old_beta
    }else{
      r_tilde = y
    }
    
    if(club_U){
      project_on_V = function(r){
        
        mean_term_r = (tmat((data.frame(old_V[f2,,drop = FALSE])*(r)/as.numeric(old_sigma2e)),f1))
        
        project_i = function(i){
          Variance_i = solve(matrix(sum_vvT[i,],K,K) + diag(K))
          return(solve(Variance_i)%*%mean_term_r[i,])    
        }
        
        E_U_r = sapply(1:R,project_i)
        
        return(colSums(t(E_U_r)[f1,]*(old_V)[f2,,drop = FALSE]))
      }
      
      
      X_residual = X
      
      for(i in 1:q){
        X_residual[,i] = X_residual[,i] - project_on_V(X[,i])
      }
      
      new_beta = solve(t(X)%*% X_residual) %*% t(X_residual) %*% y
      r_tilde = y-X%*%new_beta
      new_f_x = X%*%new_beta
    }
    
    
    negative_log_likelihood = N*(log(old_sigma2e)) + sum((r_tilde^2)/old_sigma2e)
    E_U = matrix(0, R, K)
    E_U_U_t = matrix(0, R, K^2)
    ## E-step
    #sum_vvT =  sum_z_z_T(old_V[f2,]/sqrt(old_Psi[f2]),f1)
    sum_vvT =  sum_z_z_T(old_V[f2,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f1)
    #mean_term = (tmat((data.frame(old_V[f2,])*r_tilde/old_Psi[f2]),f1))
    mean_term = (tmat((data.frame(old_V[f2,,drop = FALSE])*(r_tilde)/as.numeric(old_sigma2e)),f1))
    for(i in 1:R){
      variance_u_i = solve(matrix(sum_vvT[i,],K,K) + diag(K))
      E_U[i,] = variance_u_i %*% mean_term[i,]
      E_U_U_t[i,] = as.vector(variance_u_i + E_U[i,]%*% t(E_U[i,]))
      negative_log_likelihood = negative_log_likelihood + log(det(matrix(sum_vvT[i,],K,K) + diag(K))) - t(mean_term[i,])%*%E_U[i,]           
      
    }
    Variance = tmat(E_U_U_t[f1,],f2)
    
    if(club_V){
      project_on_U = function(r){
        
        mean_term_r = tmat((data.frame(E_U[f1,])*r),f2)
        
        project_j = function(j){
          Variance_j = matrix(Variance[j,],K,K)
          return(solve(Variance_j)%*%mean_term_r[j,])    
        }
        
        V_r = sapply(1:C,project_j)
        return(colSums(E_U[f1,]*t(V_r)[f2,]))
        
        
      }
      
      
      X_residual = X
      
      for(i in 1:q){
        X_residual[,i] = X_residual[,i] - project_on_U(X[,i])
      }
      
      new_beta = solve(t(X)%*% X_residual) %*% t(X_residual) %*% y
      r_tilde = y - X%*% new_beta  
      new_f_x = X%*% new_beta
    }
    sum_r_tilde_square_i = tmat((r_tilde^2),f2)
    mean_term = tmat((data.frame(E_U[f1,])*r_tilde),f2)
    
    ## V and Psi update
    
    
    new_V = old_V
    #new_Psi = old_Psi
    new_sigma2e = 0
    for(j in 1:C){
      Variance_j = matrix(Variance[j,],K,K)
      new_V[j,] = solve(Variance_j)%*%mean_term[j,]
      new_sigma2e = new_sigma2e + sum_r_tilde_square_i[j] - t(new_V[j,])%*%Variance_j%*%new_V[j,]
      #new_Psi[j] = sum_r_tilde_square_i[j]/m_js[j] - t(new_V[j,])%*%Variance_j%*%new_V[j,]/m_js[j]
    }
    new_sigma2e = new_sigma2e/N
    ## beta update
    new_f_v = rowSums(E_U[f1,,drop = FALSE]*new_V[f2,,drop = FALSE])
    change_fv = sqrt(sum((new_f_v - f_v)^2))/sqrt(sum((f_v)^2))
    Changes_in_f_v = c(Changes_in_f_v, change_fv)
    f_v = new_f_v
    r_tilde = y - rowSums(E_U[f1,,drop = FALSE]*new_V[f2,,drop = FALSE])
    change_V = sqrt(sum((new_V - old_V)^2))/sqrt(sum(old_V^2))
    change_sigma2e = abs(new_sigma2e/old_sigma2e - 1)
    
    
    Changes_in_V = c(Changes_in_V, change_V)
    Changes_in_sigma2e = c(Changes_in_sigma2e, change_sigma2e)
    
    
    change = max(change_fv, change_V, change_sigma2e)
    old_sigma2e = as.numeric(new_sigma2e)
    old_V = new_V
    
    ngl = c(ngl, negative_log_likelihood)
    if(!is.null(X)){
      new_beta = lm(r_tilde~X-1)$coef
      new_f_x = X%*% new_beta
      change_f_x = sqrt(sum((new_f_x - f_x)^2))/sqrt(sum((f_x)^2))
      Changes_in_f_x = c(Changes_in_f_x, change_f_x)
      change_beta = sqrt(sum((new_beta - old_beta)^2))/sqrt(sum(old_beta^2))
      change = max(change, change_beta,change_f_x)
      Changes_in_beta = c(Changes_in_beta, change_beta)
      old_beta = new_beta
      f_x = X%*% old_beta
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V[f2,,drop = FALSE])-X%*%old_beta)^2)
      training_error = c(training_error, train_error)
      residuals = y_test - X_test %*% old_beta
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
          }else{
            U = E_U[f1_test[i],]
          }
          if(f2_test[i]==0){
            V = 0
          }else{
            V = old_V[f2_test[i],]
          }
          residuals[i] = residuals[i] - sum(U*V)
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "\n") 
      }
    }else{
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V[f2,,drop = FALSE]))^2)
      training_error = c(training_error, train_error)
      residuals = y_test
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
          }else{
            U = E_U[f1_test[i],]
          }
          if(f2_test[i]==0){
            V = 0
          }else{
            V = old_V[f2_test[i],]
          }
          residuals[i] = residuals[i] - sum(U*V)
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "\n") 
      }
    }
    iter = iter + 1    
  }
  end.time = Sys.time()
  time_taken = end.time - start.time
  if(is.null(X)){
    return(list(time_taken = time_taken, V = new_V, sigma2e = new_sigma2e, U = E_U, ngl = ngl, training_error = training_error, prediction_error = prediction_error, Changes_in_f_v = Changes_in_f_v, Changes_in_V = Changes_in_V, Changes_in_sigma2e = Changes_in_sigma2e))
    
  }else{
    return(list(time_taken = time_taken, beta = new_beta, V = new_V, sigma2e = new_sigma2e, U = E_U, ngl = ngl, training_error = training_error, prediction_error = prediction_error, Changes_in_beta = Changes_in_beta, Changes_in_V = Changes_in_V, Changes_in_f_v = Changes_in_f_v, Changes_in_f_x = Changes_in_f_x, Changes_in_sigma2e = Changes_in_sigma2e))
  }
}


factor_analysis_without_intercepts_both_sides = function(y,f1,f2,K,X = NULL, initial_V = NULL, initial_Lam = NULL, initial_beta = NULL, initial_sigma2e = NULL, club_V = FALSE, club_U = FALSE, max_iter = 200, conv_thres = 10^(-5), test_data = FALSE, n_test = NULL , X_test = NULL, y_test = NULL, f1_test = NULL, f2_test = NULL){
  start.time = Sys.time()
  prediction_error = c()
  training_error = c()
  Changes_in_beta = c()
  Changes_in_sigma2e = c()
  Changes_in_V = c()
  Changes_in_f_v = c()
  Changes_in_f_lam = c()
  Changes_in_Lam = c()
  
  N=length(f1)
  R = length(unique(f1))
  C = length(unique(f2))
  n_is=tmat(rep(1,N),f1)
  m_js=tmat(rep(1,N),f2)
  #print(initial_beta==NULL)
  ## initialise V to be from Standard Normal distribution and covariance matrix Psi from chi square distribution, 
  ## beta to be ols.
  if(!is.null(X)){
    q = ncol(X)
    if(!is.null(initial_beta)){
      old_beta = initial_beta 
    }else{
      fit0 = lm(y~X-1)
      old_beta=fit0$coef
    }
    f_x = X%*% old_beta
    Changes_in_f_x = c()
  }else{
    f_x = 0
  }
  
  if(!is.null(initial_V)){
    old_V = initial_V
  }else{
    old_V=matrix(rnorm(C*K),C,K)
  }
  
  if(!is.null(initial_Lam)){
    old_Lam = initial_lam
  }else{
    old_Lam = matrix(rnorm(R*K),R,K)
  }
  if(!is.null(initial_sigma2e)){
    old_sigma2e = initial_sigma2e   
  }else{
    old_sigma2e = 1
  }
  change = 1
  iter = 1
  f_v = 0
  f_lam = 0
  
  while(iter<=max_iter && change>=conv_thres){
    sum_vvT =  sum_z_z_T(old_V[f2,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f1)
    if(!(is.null(X))){
      r_tilde = y-X%*%old_beta - f_lam
    }else{
      r_tilde = y - f_lam
    }
    
    
    E_U = matrix(0, R, K)
    E_U_U_t = matrix(0, R, K^2)
    ## E-step
    #sum_vvT =  sum_z_z_T(old_V[f2,]/sqrt(old_Psi[f2]),f1)
    sum_vvT =  sum_z_z_T(old_V[f2,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f1)
    #mean_term = (tmat((data.frame(old_V[f2,])*r_tilde/old_Psi[f2]),f1))
    mean_term = (tmat((data.frame(old_V[f2,,drop = FALSE])*(r_tilde)/as.numeric(old_sigma2e)),f1))
    for(i in 1:R){
      variance_u_i = solve(matrix(sum_vvT[i,],K,K) + diag(K))
      E_U[i,] = variance_u_i %*% mean_term[i,]
      E_U_U_t[i,] = as.vector(variance_u_i + E_U[i,]%*% t(E_U[i,]))
      
    }
    Variance_U = tmat(E_U_U_t[f1,],f2)
    
    
    sum_r_tilde_square_i = tmat((r_tilde^2),f2)
    mean_term = tmat((data.frame(E_U[f1,])*r_tilde),f2)
    
    ## V and Psi update
    
    
    new_V = old_V
    #new_Psi = old_Psi
    new_sigma2e = 0
    for(j in 1:C){
      Variance_U_j = matrix(Variance_U[j,],K,K)
      new_V[j,] = solve(Variance_U_j)%*%mean_term[j,]
      new_sigma2e = new_sigma2e + t(new_V[j,])%*%Variance_U_j%*%new_V[j,]
      #new_Psi[j] = sum_r_tilde_square_i[j]/m_js[j] - t(new_V[j,])%*%Variance_j%*%new_V[j,]/m_js[j]
    }
    #new_sigma2e = new_sigma2e/N
    ## beta update
    new_f_v = rowSums(E_U[f1,,drop = FALSE]*new_V[f2,,drop = FALSE])
    change_fv = sqrt(sum((new_f_v - f_v)^2))/sqrt(sum((f_v)^2))
    Changes_in_f_v = c(Changes_in_f_v, change_fv)
    f_v = new_f_v
    change_V = sqrt(sum((new_V - old_V)^2))/sqrt(sum((old_V)^2))
    Changes_in_V = c(Changes_in_V, change_V)
    
    sum_lam_lam_T =  sum_z_z_T(old_Lam[f1,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f2)
    if(!(is.null(X))){
      r_tilde = y-X%*%old_beta - f_v
    }else{
      r_tilde = y - f_v
    }
    
    
    E_gam = matrix(0, C, K)
    E_gam_gam_t = matrix(0, C, K^2)
    ## E-step
    #sum_vvT =  sum_z_z_T(old_V[f2,]/sqrt(old_Psi[f2]),f1)
    sum_lam_lamT =  sum_z_z_T(old_Lam[f1,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f2)
    #mean_term = (tmat((data.frame(old_V[f2,])*r_tilde/old_Psi[f2]),f1))
    mean_term = (tmat((data.frame(old_Lam[f1,,drop = FALSE])*(r_tilde)/as.numeric(old_sigma2e)),f2))
    for(i in 1:C){
      variance_gam_i = solve(matrix(sum_lam_lamT[i,],K,K) + diag(K))
      E_gam[i,] = variance_u_i %*% mean_term[i,]
      E_gam_gam_t[i,] = as.vector(variance_gam_i + E_gam[i,]%*% t(E_gam[i,]))
      
    }
    Variance_gam = tmat(E_gam_gam_t[f2,],f1)
    
    
    mean_term = tmat((data.frame(E_gam[f2,])*r_tilde),f1)
    
    ## V and Psi update
    
    
    new_Lam = old_Lam
    #new_Psi = old_Psi
    
    
    #new_Psi = old_Psi
    for(j in 1:R){
      Variance_j = matrix(Variance_gam[j,],K,K)
      new_Lam[j,] = solve(Variance_j)%*%mean_term[j,]
      new_sigma2e = new_sigma2e + t(new_Lam[j,])%*%Variance_j%*%new_Lam[j,]
    }
    #new_sigma2e = new_sigma2e/N
    ## beta update
    new_f_lam = rowSums(E_gam[f2,,drop = FALSE]*new_Lam[f1,,drop = FALSE])
    change_f_lam = sqrt(sum((new_f_lam - f_lam)^2))/sqrt(sum((f_lam)^2))
    Changes_in_f_lam = c(Changes_in_f_lam, change_f_lam)
    f_lam = new_f_lam
    
    change_Lam = sqrt(sum((new_Lam - old_Lam)^2))/sqrt(sum(old_Lam^2))
    Changes_in_lam = c(Changes_in_Lam, change_Lam)
    
    
    r_tilde = r_tilde - f_lam
    
    
    new_sigma2e = (new_sigma2e + sum(r_tilde^2))/N
    
    
    change_sigma2e = abs(new_sigma2e/old_sigma2e - 1)
    
    
    Changes_in_sigma2e = c(Changes_in_sigma2e, change_sigma2e)
    
    
    change = max(change_fv, change_V, change_sigma2e, change_Lam, change_f_lam)
    old_sigma2e = as.numeric(new_sigma2e)
    old_V = new_V
    
    if(!is.null(X)){
      new_beta = lm(r_tilde~X-1)$coef
      new_f_x = X%*% new_beta
      change_f_x = sqrt(sum((new_f_x - f_x)^2))/sqrt(sum((f_x)^2))
      Changes_in_f_x = c(Changes_in_f_x, change_f_x)
      change_beta = sqrt(sum((new_beta - old_beta)^2))/sqrt(sum(old_beta^2))
      change = max(change, change_beta,change_f_x)
      Changes_in_beta = c(Changes_in_beta, change_beta)
      old_beta = new_beta
      f_x = X%*% old_beta
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V[f2,,drop = FALSE])-X%*%old_beta - f_lam)^2)
      training_error = c(training_error, train_error)
      residuals = y_test - X_test %*% old_beta
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
            new_lam = 0
          }else{
            U = E_U[f1_test[i],]
            new_lam = new_Lam[f1_test[i]]  
          }
          if(f2_test[i]==0){
            V = 0
            new_gam = 0
          }else{
            V = old_V[f2_test[i],]
            new_gam = E_gam[f2_test[i],]
          }
          residuals[i] = residuals[i] - sum(U*V) - sum(new_lam*new_gam)
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "training_error = ",train_error, "\n") 
      }
    }else{
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V[f2,,drop = FALSE]) - f_lam)^2)
      training_error = c(training_error, train_error)
      residuals = y_test
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
            new_lam = 0
            U = E_U[f1_test[i],]
            new_lam = new_Lam[f1_test[i]]  
          }
          if(f2_test[i]==0){
            V = 0
            new_gam = 0
          }else{
            V = old_V[f2_test[i],]
            new_gam = E_gam[f2_test[i],]
          }
          residuals[i] = residuals[i] - sum(U*V) - sum(new_lam*new_gam)
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e,  "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "training_error = ",train_error, "\n") 
      }
    }
    iter = iter + 1    
  }
  end.time = Sys.time()
  time_taken = end.time - start.time
  if(is.null(X)){
    return(list(time_taken = time_taken, V = new_V, sigma2e = new_sigma2e, U = E_U, training_error = training_error, prediction_error = prediction_error, Changes_in_f_v = Changes_in_f_v, Changes_in_V = Changes_in_V, Changes_in_sigma2e = Changes_in_sigma2e))
    
  }else{
    return(list(time_taken = time_taken, beta = new_beta, V = new_V, sigma2e = new_sigma2e, U = E_U,  training_error = training_error, prediction_error = prediction_error, Changes_in_beta = Changes_in_beta, Changes_in_V = Changes_in_V, Changes_in_f_v = Changes_in_f_v, Changes_in_f_x = Changes_in_f_x, Changes_in_sigma2e = Changes_in_sigma2e))
  }
}

factor_analysis_with_client_intercepts = function(y,f1,f2,K,X = NULL, initial_V = NULL, initial_v = NULL, initial_beta = NULL, initial_sigma2e = NULL, club_U = FALSE, max_iter = 200, conv_thres = 10^(-5), test_data = FALSE, n_test = NULL , X_test = NULL, y_test = NULL, f1_test = NULL, f2_test = NULL){
  start.time = Sys.time()
  prediction_error = c()
  training_error = c()
  Changes_in_beta = c()
  Changes_in_sigma2e = c()
  Changes_in_V = c()
  Changes_in_f_V = c()
  Changes_in_v = c()
  Changes_in_f_v = c()
  
  N=length(f1)
  R = length(unique(f1))
  C = length(unique(f2))
  n_is=tmat(rep(1,N),f1)
  m_js=tmat(rep(1,N),f2)
  #print(initial_beta==NULL)
  ## initialise V to be from Standard Normal distribution and covariance matrix Psi from chi square distribution, 
  ## beta to be ols.
  if(!is.null(X)){
    q = ncol(X)
    if(!is.null(initial_beta)){
      old_beta = initial_beta 
    }else{
      fit0 = lm(y~X-1)
      old_beta=fit0$coef
    }
    f_x = X%*% old_beta
    Changes_in_f_x = c()
  }else{
    f_x = 0
  }
  
  if(!is.null(initial_V)){
    old_V = initial_V
  }else{
    old_V=matrix(rnorm(C*K),C,K)
  }
  if(!is.null(initial_v)){
    old_v = initial_v
  }else{
    old_v=rnorm(1)
  }
  old_V_tilde = cbind(old_v, old_V)
  if(!is.null(initial_sigma2e)){
    old_sigma2e = initial_sigma2e   
  }else{
    old_sigma2e = 1
  }
  change = 1
  iter = 1
  f_v = 0
  
  while(iter<=max_iter && change>=conv_thres){
    sum_vvT =  sum_z_z_T(old_V_tilde[f2,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f1)
    if(!(is.null(X))){
      r_tilde = y-X%*%old_beta
    }else{
      r_tilde = y
    }
    
    if(club_U){
      project_on_V = function(r){
        
        mean_term_r = (tmat((data.frame(old_V_tilde[f2,,drop = FALSE])*(r)/as.numeric(old_sigma2e)),f1))
        
        project_i = function(i){
          Variance_i = solve(matrix(sum_vvT[i,],K+1,K+1) + diag(K+1))
          return(solve(Variance_i)%*%mean_term_r[i,])    
        }
        
        E_U_r = sapply(1:R,project_i)
        
        return(colSums(t(E_U_r)[f1,]*(old_V_tilde)[f2,,drop = FALSE]))
      }
      
      
      X_residual = X
      
      for(i in 1:q){
        X_residual[,i] = X_residual[,i] - project_on_V(X[,i])
      }
      
      new_beta = solve(t(X)%*% X_residual) %*% t(X_residual) %*% y
      r_tilde = y-X%*%new_beta
      new_f_x = X%*%new_beta
    }
    
    
    negative_log_likelihood = N*(log(old_sigma2e)) + sum((r_tilde^2)/old_sigma2e)
    E_U = matrix(0, R, K+1)
    E_U_U_t = matrix(0, R, (K+1)^2)
    ## E-step
    sum_vvT =  sum_z_z_T(old_V_tilde[f2,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f1)
    mean_term = (tmat((data.frame(old_V_tilde[f2,,drop = FALSE])*(r_tilde)/as.numeric(old_sigma2e)),f1))
    for(i in 1:R){
      variance_u_i = solve(matrix(sum_vvT[i,],K+1,K+1) + diag(K+1))
      E_U[i,] = variance_u_i %*% mean_term[i,]
      E_U_U_t[i,] = as.vector(variance_u_i + E_U[i,]%*% t(E_U[i,]))
      negative_log_likelihood = negative_log_likelihood + log(det(matrix(sum_vvT[i,],K+1,K+1) + diag(K+1))) - t(mean_term[i,])%*%E_U[i,]           
      
    }
    Variance = tmat(E_U_U_t[f1,],f2)
    

    sum_r_tilde_square_i = tmat((r_tilde^2),f2)
    mean_term = tmat((data.frame(E_U[f1,])*r_tilde),f2)
    
    ## V and Psi update
    
    
    new_V = old_V
    c_j = 0
    alpha_j = 0
    #new_Psi = old_Psi
    new_sigma2e = 0
    for(j in 1:C){
      Variance_j = matrix(Variance[j,],K+1,K+1)
      new_V[j,] = solve(Variance_j[2:(K+1),2:(K+1)])%*%(mean_term[j,2:(K+1)] - old_v * Variance_j[2:(K+1),1])
      c_j = c_j + mean_term[j,1] - sum(new_V[j,]* Variance_j[2:(K+1),1])
      alpha_j = alpha_j + Variance_j[1,1]
      new_sigma2e = new_sigma2e + (sum_r_tilde_square_i[j] + t(c(old_v, new_V[j,]))%*%Variance_j%*%c(old_v, new_V[j,]) - 2*t(c(old_v, new_V[j,]))%*%mean_term[j,])
     }   
    new_sigma2e = new_sigma2e/N
    
    new_v = c_j/ alpha_j 
    new_V_tilde = cbind(new_v, new_V)
    
    ## beta update
    new_f_v = rowSums(E_U[f1,,drop = FALSE]*new_V_tilde[f2,,drop = FALSE])
    change_fv = sqrt(sum((new_f_v - f_v)^2))/sqrt(sum((f_v)^2))
    Changes_in_f_v = c(Changes_in_f_v, change_fv)
    f_v = new_f_v
    r_tilde = y - rowSums(E_U[f1,,drop = FALSE]*new_V_tilde[f2,,drop = FALSE])
    change_V = sqrt(sum((new_V_tilde - old_V_tilde)^2))/sqrt(sum(old_V_tilde^2))
    change_sigma2e = abs(new_sigma2e/old_sigma2e - 1)
    
    
    Changes_in_V = c(Changes_in_V, change_V)
    Changes_in_sigma2e = c(Changes_in_sigma2e, change_sigma2e)
    
    
    change = max(change_fv, change_V, change_sigma2e)
    old_sigma2e = as.numeric(new_sigma2e)
    old_V_tilde = new_V_tilde
    old_V = new_V
    old_v = new_v

    ngl = c(ngl, negative_log_likelihood)
    if(!is.null(X)){
      new_beta = lm(r_tilde~X-1)$coef
      new_f_x = X%*% new_beta
      change_f_x = sqrt(sum((new_f_x - f_x)^2))/sqrt(sum((f_x)^2))
      Changes_in_f_x = c(Changes_in_f_x, change_f_x)
      change_beta = sqrt(sum((new_beta - old_beta)^2))/sqrt(sum(old_beta^2))
      change = max(change, change_beta,change_f_x)
      Changes_in_beta = c(Changes_in_beta, change_beta)
      old_beta = new_beta
      f_x = X%*% old_beta
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V_tilde[f2,,drop = FALSE])-X%*%old_beta)^2)
      training_error = c(training_error, train_error)
      residuals = y_test - X_test %*% old_beta
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
          }else{
            U = E_U[f1_test[i],]
          }
          if(f2_test[i]==0){
            V = 0
          }else{
            V = old_V_tilde[f2_test[i],]
          }
          residuals[i] = residuals[i] - sum(U*V)
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "\n") 
      }
    }else{
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V_tilde[f2,,drop = FALSE]))^2)
      training_error = c(training_error, train_error)
      residuals = y_test
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
          }else{
            U = E_U[f1_test[i],]
          }
          if(f2_test[i]==0){
            V = 0
          }else{
            V = old_V_tilde[f2_test[i],]
          }
          residuals[i] = residuals[i] - sum(U*V)
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "ngl = ", negative_log_likelihood, "training_error = ",train_error, "\n") 
      }
    }
    iter = iter + 1    
  }
  end.time = Sys.time()
  time_taken = end.time - start.time
  if(is.null(X)){
    return(list(time_taken = time_taken, V = new_V, v = new_v, sigma2e = new_sigma2e, U = E_U, ngl = ngl, training_error = training_error, prediction_error = prediction_error, Changes_in_f_v = Changes_in_f_v, Changes_in_V = Changes_in_V, Changes_in_sigma2e = Changes_in_sigma2e))
    
  }else{
    return(list(time_taken = time_taken, beta = new_beta, V = new_V,v = new_v, sigma2e = new_sigma2e, U = E_U, ngl = ngl, training_error = training_error, prediction_error = prediction_error, Changes_in_beta = Changes_in_beta, Changes_in_V = Changes_in_V, Changes_in_f_v = Changes_in_f_v, Changes_in_f_x = Changes_in_f_x, Changes_in_sigma2e = Changes_in_sigma2e))
  }
}

factor_analysis_with_client_intercepts_both_sides = function(y,f1,f2,K,X = NULL, initial_V = NULL, initial_v = NULL, initial_Lam = NULL, initial_lam = NULL, initial_beta = NULL, initial_sigma2e = NULL, club_U = FALSE, max_iter = 200, conv_thres = 10^(-5), test_data = FALSE, n_test = NULL , X_test = NULL, y_test = NULL, f1_test = NULL, f2_test = NULL){
  start.time = Sys.time()
  prediction_error = c()
  training_error = c()
  Changes_in_beta = c()
  Changes_in_sigma2e = c()
  Changes_in_V = c()
  Changes_in_f_V = c()
  Changes_in_f_lam = c()
  Changes_in_Lam = c()
  
  N=length(f1)
  R = length(unique(f1))
  C = length(unique(f2))
  n_is=tmat(rep(1,N),f1)
  m_js=tmat(rep(1,N),f2)
  #print(initial_beta==NULL)
  ## initialise V to be from Standard Normal distribution and covariance matrix Psi from chi square distribution, 
  ## beta to be ols.
  if(!is.null(X)){
    q = ncol(X)
    if(!is.null(initial_beta)){
      old_beta = initial_beta 
    }else{
      fit0 = lm(y~X-1)
      old_beta=fit0$coef
    }
    f_x = X%*% old_beta
    Changes_in_f_x = c()
  }else{
    f_x = 0
  }
  
  if(!is.null(initial_V)){
    old_V = initial_V
  }else{
    old_V=matrix(rnorm(C*K),C,K)
  }
  if(!is.null(initial_v)){
    old_v = initial_v
  }else{
    old_v = rnorm(1)
  }
  if(!is.null(initial_Lam)){
    old_Lam = initial_lam
  }else{
    old_Lam = matrix(rnorm(R*K),R,K)
  }
  if(!is.null(initial_lam)){
    old_lam = initial_lam
  }else{
    old_lam = rnorm(1)
  }
  old_V_tilde = cbind(old_v, old_V)
  old_Lam_tilde = cbind(old_lam, old_Lam)
  
  if(!is.null(initial_sigma2e)){
    old_sigma2e = initial_sigma2e   
  }else{
    old_sigma2e = 1
  }
  change = 1
  iter = 1
  f_v = 0
  f_lam = 0
  
  while(iter<=max_iter && change>=conv_thres){
    if(!(is.null(X))){
      r_tilde = y-X%*%old_beta - f_lam
    }else{
      r_tilde = y - f_lam
    }
    E_U = matrix(0, R, K+1)
    E_U_U_t = matrix(0, R, (K+1)^2)
    ## E-step
    sum_vvT =  sum_z_z_T(old_V_tilde[f2,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f1)
    mean_term = (tmat((data.frame(old_V_tilde[f2,,drop = FALSE])*(r_tilde)/as.numeric(old_sigma2e)),f1))
    for(i in 1:R){
      variance_u_i = solve(matrix(sum_vvT[i,],K+1,K+1) + diag(K+1))
      E_U[i,] = variance_u_i %*% mean_term[i,]
      E_U_U_t[i,] = as.vector(variance_u_i + E_U[i,]%*% t(E_U[i,]))
      
    }
    Variance = tmat(E_U_U_t[f1,],f2)
    
    
    mean_term = tmat((data.frame(E_U[f1,])*r_tilde),f2)
    
    ## V and Psi update
    
    
    new_V = old_V
    c_j = 0
    alpha_j = 0
    #new_Psi = old_Psi
    new_sigma2e = 0
    for(j in 1:C){
      Variance_j = matrix(Variance[j,],K+1,K+1)
      new_V[j,] = solve(Variance_j[2:(K+1),2:(K+1)])%*%(mean_term[j,2:(K+1)] - old_v * Variance_j[2:(K+1),1])
      c_j = c_j + mean_term[j,1] - sum(new_V[j,]* Variance_j[2:(K+1),1])
      alpha_j = alpha_j + Variance_j[1,1]
      new_sigma2e = new_sigma2e + t(c(old_v, new_V[j,]))%*%Variance_j%*%c(old_v, new_V[j,]) 
    }   
    new_v = c_j/ alpha_j 
    new_V_tilde = cbind(new_v, new_V)
    
    new_f_v = rowSums(E_U[f1,,drop = FALSE]*new_V_tilde[f2,,drop = FALSE])
    change_fv = sqrt(sum((new_f_v - f_v)^2))/sqrt(sum((f_v)^2))
    Changes_in_f_V = c(Changes_in_f_V, change_fv)
    f_v = new_f_v
    change_V = sqrt(sum((new_V_tilde - old_V_tilde)^2))/sqrt(sum(old_V_tilde^2))
    
    
    Changes_in_V = c(Changes_in_V, change_V)
    
    sum_lam_lamT =  sum_z_z_T(old_Lam_tilde[f1,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f2)
    if(!(is.null(X))){
      r_tilde = y-X%*%old_beta - f_v
    }else{
      r_tilde = y - f_v
    }
    
    E_gam = matrix(0, C, K+1)
    E_gam_gam_t = matrix(0, C, (K+1)^2)
    ## E-step
    mean_term = (tmat((data.frame(old_Lam_tilde[f1,,drop = FALSE])*(r_tilde)/as.numeric(old_sigma2e)),f2))
    for(i in 1:C){
      variance_u_i = solve(matrix(sum_lam_lamT[i,],K+1,K+1) + diag(K+1))
      E_gam[i,] = variance_u_i %*% mean_term[i,]
      E_gam_gam_t[i,] = as.vector(variance_u_i + E_gam[i,]%*% t(E_gam[i,]))
      
    }
    Variance = tmat(E_gam_gam_t[f2,],f1)
    
    
    mean_term = tmat((data.frame(E_gam[f2,])*r_tilde),f1)
    
    new_Lam = old_Lam
    c_j = 0
    alpha_j = 0
    #new_Psi = old_Psi
    for(j in 1:R){
      Variance_j = matrix(Variance[j,],K+1,K+1)
      new_Lam[j,] = solve(Variance_j[2:(K+1),2:(K+1)])%*%(mean_term[j,2:(K+1)] - old_lam * Variance_j[2:(K+1),1])
      c_j = c_j + mean_term[j,1] - sum(new_Lam[j,]* Variance_j[2:(K+1),1])
      alpha_j = alpha_j + Variance_j[1,1]
      new_sigma2e = new_sigma2e + t(c(old_lam, new_Lam[j,]))%*%Variance_j%*%c(old_lam, new_Lam[j,]) 
    }   
    
    new_lam = c_j/ alpha_j 
    new_Lam_tilde = cbind(new_lam, new_Lam)
    
    ## beta update
    new_f_lam = rowSums(E_gam[f2,,drop = FALSE]*new_Lam_tilde[f1,,drop = FALSE])
    change_f_lam = sqrt(sum((new_f_lam - f_lam)^2))/sqrt(sum((f_lam)^2))
    Changes_in_f_lam = c(Changes_in_f_lam, change_f_lam)
    f_lam = new_f_lam   
    change_lam = sqrt(sum((new_Lam_tilde - old_Lam_tilde)^2))/sqrt(sum(old_Lam_tilde^2))
    old_Lam = new_Lam
    old_lam = new_lam
    old_Lam_tilde = new_Lam_tilde
    
    
    Changes_in_Lam = c(Changes_in_Lam, change_lam)
    
    r_tilde = r_tilde + f_v
    new_sigma2e = (new_sigma2e + sum(r_tilde^2) -2* sum(r_tilde*(f_v + f_lam)))/N
    change_sigma2e = abs(new_sigma2e/old_sigma2e - 1)
    Changes_in_sigma2e = c(Changes_in_sigma2e, change_sigma2e)
    
    
    change = max(change_fv, change_V, change_sigma2e, change_f_lam, change_lam)
    old_sigma2e = as.numeric(new_sigma2e)
    old_V_tilde = new_V_tilde
    old_V = new_V
    old_v = new_v
    
    r_tilde = y - f_v - f_lam
    
    if(!is.null(X)){
      #print(old_beta)
      new_beta = lm(r_tilde~X-1)$coef
      new_f_x = X%*% new_beta
      change_f_x = sqrt(sum((new_f_x - f_x)^2))/sqrt(sum((f_x)^2))
      Changes_in_f_x = c(Changes_in_f_x, change_f_x)
      change_beta = sqrt(sum((new_beta - old_beta)^2))/sqrt(sum(old_beta^2))
      change = max(change, change_beta,change_f_x)
      Changes_in_beta = c(Changes_in_beta, change_beta)
      old_beta = new_beta
      f_x = X%*% old_beta
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V_tilde[f2,,drop = FALSE])-X%*%old_beta - f_lam)^2)
      training_error = c(training_error, train_error)
      residuals = y_test - X_test %*% as.matrix(old_beta)
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
            l = 0
          }else{
            U = E_U[f1_test[i],]
            l = old_Lam_tilde[f1_test[i], ]
          }
          if(f2_test[i]==0){
            V = 0
            g = 0
          }else{
            V = old_V_tilde[f2_test[i],]
            g = E_gam[f2_test[i], ]
          }
          residuals[i] = residuals[i] - sum(U*V) - sum(l*g)
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "training_error = ",train_error, "\n") 
      }
    }else{
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V_tilde[f2,,drop = FALSE]) - f_lam)^2)
      training_error = c(training_error, train_error)
      residuals = y_test
      if(test_data){
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
            l = 0
          }else{
            U = E_U[f1_test[i],]
            l = old_Lam_tilde[f1_test[i], ]
          }
          if(f2_test[i]==0){
            V = 0
            g = 0
          }else{
            V = old_V_tilde[f2_test[i],]
            g = E_gam[f2_test[i], ]
          }
          residuals[i] = residuals[i] - sum(U*V) - sum(l*g)
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "training_error = ",train_error, "\n") 
      }
    }
    iter = iter + 1    
  }
  end.time = Sys.time()
  time_taken = end.time - start.time
  if(is.null(X)){
    return(list(time_taken = time_taken, V = new_V, v = new_v, sigma2e = new_sigma2e, U = E_U, Lam = new_Lam_tilde, gam = E_gam, training_error = training_error, prediction_error = prediction_error, Changes_in_f_v = Changes_in_f_v, Changes_in_V = Changes_in_V, Changes_in_sigma2e = Changes_in_sigma2e, Changes_in_f_lam = Changes_in_f_lam, Changes_in_Lam = Changes_in_Lam))
    
  }else{
    return(list(time_taken = time_taken, beta = new_beta, V = new_V,v = new_v, sigma2e = new_sigma2e, U = E_U, Lam = new_Lam_tilde, gam = E_gam, training_error = training_error, prediction_error = prediction_error, Changes_in_beta = Changes_in_beta, Changes_in_V = Changes_in_V, Changes_in_f_v = Changes_in_f_v, Changes_in_f_x = Changes_in_f_x, Changes_in_sigma2e = Changes_in_sigma2e , Changes_in_f_lam = Changes_in_f_lam, Changes_in_Lam = Changes_in_Lam))
  }
}


initial_V = NULL; initial_v = NULL; initial_Lam = NULL; initial_lam = NULL; initial_beta = NULL; initial_sigma2e = NULL; club_U = FALSE; max_iter = 200; conv_thres = 10^(-5); test_data = TRUE

fit = factor_analysis_with_client_intercepts_both_sides(y, f1, f2, 2, X = X, test_data = TRUE, n_test = n_test, X_test =  X_test, y_test = y_test, f1_test = f1_test, f2_test = f2_test)




factor_analysis_with_client_item_intercept = function(y,f1,f2,K,X = NULL, initial_V = NULL, initial_beta = NULL, initial_A = NULL,initial_B = NULL, initial_sigma2e = NULL, initial_sigma2a = NULL, initial_sigma2b = NULL, club_V = FALSE, club_U = FALSE, club_a = FALSE, club_b = FALSE, max_iter = 200, conv_thres = 10^(-5), test_data = FALSE, n_test = NULL , X_test = NULL, y_test = NULL, f1_test = NULL, f2_test = NULL){
  start.time = Sys.time()
  prediction_error = c()
  training_error = c()
  Changes_in_beta = c()
  Changes_in_sigma2e = c()
  Changes_in_V = c()
  Changes_in_f_v = c()
  Changes_in_sigma2a = c()
  Changes_in_sigma2b = c()
  Changes_in_f_a = c()
  Changes_in_f_b = c()
  
  
  N=length(f1)
  q = ncol(X)
  R = length(unique(f1))
  C = length(unique(f2))
  n_is=tmat(rep(1,N),f1)
  m_js=tmat(rep(1,N),f2)
  #print(initial_beta==NULL)
  ## initialise V to be from Standard Normal distribution and covariance matrix Psi from chi square distribution, 
  ## beta to be ols.
  if(!is.null(X)){
    if(!(is.null(initial_beta))){
      old_beta = initial_beta 
      f_x = X%*%old_beta
    }else{
      fit0 = lm(y~X-1)
      old_beta=fit0$coef
      f_x = X%*%old_beta
      new_f_x = f_x
      Changes_in_f_x = c()
    } 
  }else{
    f_x = 0
    new_f_x = f_x
  }
  
  if(!is.null(initial_V)){
    old_V = initial_V
  }else{
    old_V=matrix(rnorm(C*K),C,K)
  }
  if(!is.null(initial_sigma2e)){
    old_sigma2e = as.numeric(initial_sigma2e)   
  }else{
    old_sigma2e = 1
  }
  if(!(is.null(initial_sigma2a))){
    old_sigma2a = as.numeric(initial_sigma2a)   
  }else{
    old_sigma2a = 1
  }
  if(!is.null(initial_A)){
    old_A = initial_A  
  }else{
    old_A = rnorm(R)
  }
  if(!is.null(initial_sigma2b)){
    old_sigma2b = as.numeric(initial_sigma2b)  
  }else{
    old_sigma2b = 1
  }
  if(!is.null(initial_B)){
    old_B = initial_B  
  }else{
    old_B = rnorm(C)
  }
  change = 1
  iter = 1
  f_v = rep(0,N)
  f_a = old_A[f1]
  f_b = old_B[f2]
  
  while(iter<=max_iter && change>=conv_thres){
    
    var_a_given_y = 1/(n_is/(old_sigma2e) + 1/ (old_sigma2a))
    var_b_given_y = 1/(m_js/(old_sigma2e) + 1/ (old_sigma2b))
    
    
    
    if(club_a){
      fit_a = function(r){
        return(((tmat(r,f1))/(n_is + (old_sigma2e)/(old_sigma2a)))[f1])
      }
      X_residual = X
      for(i in 1:q){
        X_residual[,i] = X[,i] - fit_a(X[,i])   
      }
      new_beta = solve(t(X)%*% X_residual) %*% t(X_residual) %*% (y - f_v - f_b)
      
      new_f_x = X%*% new_beta
      
    }
    ## a update
    
    r_a = y - new_f_x - f_v - f_b
    
    E_a = (tmat(r_a,f1))/(n_is + (old_sigma2e)/(old_sigma2a))
    
    new_f_a = E_a[f1]
    
    new_sigma2a = mean((E_a)^2) + mean(var_a_given_y)
    
    if(club_b){
      fit_b = function(r){
        return(((tmat(r,f2))/(m_js + (old_sigma2e)/(old_sigma2b)))[f2])
      }
      X_residual = X
      for(i in 1:q){
        X_residual[,i] = X[,i] - fit_b(X[,i])   
      }
      new_beta = solve(t(X)%*% X_residual) %*% t(X_residual) %*% (y - f_v - f_a)
      
      new_f_x = X%*% new_beta
      
    }
    ## b update
    
    r_b = y - new_f_x - f_v - new_f_a
    
    E_b = (tmat(r_b,f2))/(m_js + (old_sigma2e)/(old_sigma2b))
    
    new_f_b = E_b[f2]
    
    new_sigma2b = mean((E_b)^2) + mean(var_b_given_y)  
    
    sum_vvT =  sum_z_z_T(old_V[f2,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f1)
    
    
    if(club_U){
      project_on_V = function(r){
        
        mean_term_r = (tmat((data.frame(old_V[f2,, drop = FALSE])*(r)/as.numeric(old_sigma2e)),f1))
        
        project_i = function(i){
          Variance_i = solve(matrix(sum_vvT[i,],K,K) + diag(K))
          return(solve(Variance_i)%*%mean_term_r[i,])    
        }
        
        E_U_r = sapply(1:R,project_i)
        
        return(colSums(t(E_U_r)[f1,]*(old_V)[f2,]))
      }
      
      
      X_residual = X
      
      for(i in 1:q){
        X_residual[,i] = X_residual[,i] - project_on_V(X[,i])
      }
      
      new_beta = solve(t(X)%*% X_residual) %*% t(X_residual) %*% (y - f_a - f_b)
      
      new_f_x = X%*%new_beta
    }
    r_tilde = y-new_f_x - new_f_a - new_f_b
    E_U = matrix(0, R, K)
    E_U_U_t = matrix(0, R, K^2)
    ## E-step
    sum_vvT =  sum_z_z_T(old_V[f2,,drop = FALSE]/sqrt(as.numeric(old_sigma2e)),f1)
    mean_term = (tmat((data.frame(old_V[f2,,drop = FALSE])*(r_tilde)/as.numeric(old_sigma2e)),f1))
    for(i in 1:R){
      variance_u_i = solve(matrix(sum_vvT[i,],K,K) + diag(K))
      E_U[i,] = variance_u_i %*% mean_term[i,]
      E_U_U_t[i,] = as.vector(variance_u_i + E_U[i,]%*% t(E_U[i,]))
      
    }
    sum_r_tilde_square_i = tmat((r_tilde^2),f2)
    Variance = tmat(E_U_U_t[f1,],f2)
    
    if(club_V){
      project_on_U = function(r){
        
        mean_term_r = tmat((data.frame(E_U[f1,])*r),f2)
        
        project_j = function(j){
          Variance_j = matrix(Variance[j,],K,K)
          return(solve(Variance_j)%*%mean_term_r[j,])    
        }
        
        V_r = sapply(1:C,project_j)
        return(colSums(E_U[f1,]*t(V_r)[f2,]))
        
        
      }
      
      
      q = ncol(X)
      X_residual = X
      
      for(i in 1:q){
        X_residual[,i] = X_residual[,i] - project_on_U(X[,i])
      }
      
      new_beta = solve(t(X)%*% X_residual) %*% t(X_residual) %*% (y - f_a - f_b)
      new_f_x = X%*% new_beta
      r_tilde = y - X%*% new_beta - new_f_a - new_f_b
      
    }
    
    
    mean_term = tmat((data.frame(E_U[f1,,drop = FALSE])*r_tilde),f2)
    
    ## V and Psi update
    
    
    new_V = old_V
    #new_Psi = old_Psi
    new_sigma2e = sum(var_a_given_y[f1] + var_b_given_y[f2])
    for(j in 1:C){
      Variance_j = matrix(Variance[j,],K,K)
      new_V[j,] = solve(Variance_j)%*%mean_term[j,]
      new_sigma2e = new_sigma2e + sum_r_tilde_square_i[j] - t(new_V[j,])%*%Variance_j%*%new_V[j,]
      #new_Psi[j] = sum_r_tilde_square_i[j]/m_js[j] - t(new_V[j,])%*%Variance_j%*%new_V[j,]/m_js[j]
    }
    new_sigma2e = new_sigma2e/N
    ## beta update
    new_f_v = rowSums(E_U[f1,,drop = FALSE]*new_V[f2,,drop = FALSE])
    change_V = sqrt(sum((new_V - old_V)^2))/sqrt(sum(old_V^2))
    change_fv = sqrt(sum((new_f_v - f_v)^2))/sqrt(sum(f_v^2))
    f_v = new_f_v
    
    Changes_in_f_v = c(Changes_in_f_v, change_fv)
    change_fa = sqrt(sum((new_f_a - f_a)^2))/sqrt(sum(f_a^2))
    Changes_in_f_a = c(Changes_in_f_a, change_fa)
    change_fb = sqrt(sum((new_f_b - f_b)^2))/sqrt(sum(f_b^2))
    Changes_in_f_b = c(Changes_in_f_b, change_fb)
    f_a = new_f_a
    f_b = new_f_b
    change_sigma2e = abs(new_sigma2e/old_sigma2e - 1)
    change_sigma2a = abs(new_sigma2a/old_sigma2a - 1)
    change_sigma2b = abs(new_sigma2b/old_sigma2b - 1)
    change = max(change_V, change_fa,change_fv, change_sigma2e, change_sigma2a, change_fb, change_sigma2b)
    
    Changes_in_V = c(Changes_in_V, change_V)
    Changes_in_sigma2e = c(Changes_in_sigma2e, change_sigma2e)
    Changes_in_sigma2a = c(Changes_in_sigma2a, change_sigma2a)
    Changes_in_sigma2b = c(Changes_in_sigma2b, change_sigma2b)
    
    old_sigma2e = as.numeric(new_sigma2e)
    old_sigma2a = as.numeric(new_sigma2a)
    old_sigma2b = as.numeric(new_sigma2b)
    old_V = new_V
    
    if(!is.null(X)){
      r_tilde = y - new_f_v - new_f_a - new_f_b
      new_beta = lm(r_tilde~X-1)$coef
      change_beta = sqrt(sum((new_beta - old_beta)^2))/sqrt(sum(old_beta^2))
      Changes_in_beta = c(Changes_in_beta, change_beta)
      old_beta = new_beta
      new_f_x = X%*% new_beta
      change_f_x = sqrt(sum((new_f_x - f_x)^2))/sqrt(sum(f_x^2))
      f_x = new_f_x
      Changes_in_f_x = c(Changes_in_f_x, change_f_x)
      change = max(change, change_beta,change_f_x)
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V[f2,,drop = FALSE])-X%*%old_beta - E_a[f1] - E_b[f2])^2)
      training_error = c(training_error, train_error)
      if(test_data){
        residuals = y_test - X_test %*% old_beta
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
            a = 0
          }else{
            U = E_U[f1_test[i],]
            a = E_a[f1_test[i]]
          }
          if(f2_test[i]==0){
            V = 0
            b = 0
          }else{
            V = old_V[f2_test[i],]
            b = E_b[f2_test[i]]
          }
          residuals[i] = residuals[i] - sum(U*V)- a - b
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "chng sigma2a", change_sigma2a, "chng fa", change_fa,"chng sigma2b", change_sigma2b, "chng fb", change_fb,  "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng beta = ", change_beta, "chng V = ",change_V, "chng sigma2e",change_sigma2e,"chng fx = ", change_f_x, "chng fv = ",change_fv, "chng sigma2a", change_sigma2a, "chng fa", change_fa,"chng sigma2b", change_sigma2b, "chng fb", change_fb, "training_error = ",train_error, "\n") 
      }
    }else{
      train_error = mean((y-rowSums(E_U[f1,,drop = FALSE]*old_V[f2,,drop = FALSE]) - E_a[f1] - E_b[f2] )^2)
      training_error = c(training_error, train_error)
      if(test_data){
        residuals = y_test 
        for(i in 1:n_test){
          if(f1_test[i]==0){
            U = 0 
            a = 0
          }else{
            U = E_U[f1_test[i],]
            a = E_a[f1_test[i]]
          }
          if(f2_test[i]==0){
            V = 0
            b = 0
          }else{
            V = old_V[f2_test[i],]
            b = E_b[f2_test[i]]
          }
          residuals[i] = residuals[i] - sum(U*V)- a - b
        }
        pred_error = mean((residuals)^2)
        
        prediction_error = c(prediction_error, pred_error)
        
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "chng sigma2a", change_sigma2a, "chng fa", change_fa, "chng sigma2b", change_sigma2b, "chng fb", change_fb,  "training_error = ",train_error, "prediction error = ",pred_error, "\n") 
      }else{
        cat(iter, "chng V = ",change_V,"chng fv = ",change_fv, "chng sigma2e",change_sigma2e, "chng sigma2a", change_sigma2a, "chng fa", change_fa, "chng sigma2b", change_sigma2b, "chng fb", change_fb,  "training_error = ",train_error, "\n") 
      }
    }
    iter = iter + 1    
  }
  end.time = Sys.time()
  time_taken = end.time - start.time
  if(is.null(X)){
    return(list(time_taken = time_taken, V = new_V, A = E_a, B = E_b, sigma2e = new_sigma2e,sigma2a = new_sigma2a, Changes_in_f_a = Changes_in_f_a, Changes_in_sigma2a = Changes_in_sigma2a, Changes_in_f_b = Changes_in_f_b, Changes_in_sigma2b = Changes_in_sigma2b, U = E_U, training_error = training_error, prediction_error = prediction_error, Changes_in_f_v = Changes_in_f_v, Changes_in_V = Changes_in_V, Changes_in_sigma2e = Changes_in_sigma2e))
    
  }else{
    return(list(time_taken = time_taken, Changes_in_f_a = Changes_in_f_a, A = E_a, B = E_b,  Changes_in_sigma2a = Changes_in_sigma2a,Changes_in_f_b = Changes_in_f_b, Changes_in_sigma2b = Changes_in_sigma2b,beta = new_beta, V = new_V, sigma2e = new_sigma2e, U = E_U, training_error = training_error, prediction_error = prediction_error, Changes_in_beta = Changes_in_beta, Changes_in_V = Changes_in_V, Changes_in_f_v = Changes_in_f_v, Changes_in_f_x = Changes_in_f_x, Changes_in_sigma2e = Changes_in_sigma2e))
  }
}



factor_analysis = function(X,y,f1,f2,K, add_client_intercept = 0, all_item_intercept = 1, club_V = 0, club_U = 0, club_a = 0, club_b = 0, max_iter = 200, conv_thres = 10^(-5), initial_V = -1,initial_beta = -1, initial_A = -1,initial_B = -1, initial_sigma2e = -1, initial_sigma2a = -1, initial_sigma2b = -1){
  if(add_client_intercept && add_item_intercept){
    fit = factor_analysis_with_client_item_intercept(X,y,f1,f2,K, initial_V, initial_beta, initial_A,initial_B, initial_sigma2e, initial_sigma2a, initial_sigma2b, club_V, club_U, club_a, club_b)
  }else if(add_client_intercept){
    fit = factor_analysis_with_client_intercept(X,y,f1,f2,K, initial_V, initial_beta, initial_A , initial_sigma2e, initial_sigma2a , club_V, club_U , club_a)
  }else{
    fit = factor_analysis_without_intercepts(X,y,f1,f2,K, initial_V, initial_beta , initial_sigma2e , club_V , club_U)
  }
  
  return(fit)
}


