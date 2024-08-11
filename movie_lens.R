source("factor_analysis.R")
library(softImpute)
user_info = read.table("movie_lens_100k_dataset/u.user", header = FALSE, sep = "|", quote = "", stringsAsFactors = FALSE)
user_info = data.frame(user_info)
colnames(user_info) = c("user_id","age","gender","occupation","zip code")
occupations = unique(user_info$occupation)
o = length(unique(user_info$occupation))
for(i in 1:o){
  if(!(i==2)){
    name = paste("occupation",occupations[i])
    user_info[[name]] = ifelse(user_info$occupation==occupations[i],1,0)
    
  }
}
movie_info = read.table("movie_lens_100k_dataset/u.item", header = FALSE, sep = "|", quote = "", stringsAsFactors = FALSE)
movie_info = movie_info[,-4]
movie_info = data.frame(movie_info)
colnames(movie_info) = c("movie_id","movie_title","video_release_date","IMDb_URL","unknown","Action","Adventure","Animation","Children's","Comedy","Crime","Documentary","Drama","Fantasy","Film-Noir","Horror","Musical","Mystery","Romance","Sci-Fi","Thriller","War","Western")




results = function(t, K, occupation,  max_iter = 200,c = 6){
    
    train_data = read.table(sprintf("movie_lens_100k_dataset/u%d.base",t), header = FALSE)
    test_data = read.table(sprintf("movie_lens_100k_dataset/u%d.test",t), header = FALSE)
    train_data = data.frame(train_data)
    colnames(train_data) = c("user_id","item_id","rating","timestamp")
    test_data = data.frame(test_data)
    colnames(test_data) = c("user_id","item_id","rating","timestamp")
    
    y = train_data$rating
    train_movie = movie_info[as.numeric(train_data$item_id), ]
    train_movie = train_movie[,6:23]
    if(occupation==1){
      train_user = user_info[train_data$user_id,c(2:3,6:25)]
      test_user = user_info[test_data$user_id,c(2:3, 6:25)]
      
    }else{
      train_user = user_info[train_data$user_id,2:3]
      test_user = user_info[test_data$user_id,2:3]
    }
    X = cbind(train_movie, train_user)
    X$gender = ifelse(X$gender == "M", 1, 0)
    #X$occupation = as.factor(X$occupation)
    X = as.matrix(X)
    N = nrow(X)
    X = cbind(1,X)
    
    
    y_test = test_data$rating
    test_movie = movie_info[as.numeric(test_data$item_id), ]
    test_movie = test_movie[,6:23]
    
    X_test = cbind(test_movie, test_user)
    X_test$gender = ifelse(X_test$gender == "M", 1, 0)
    #X$occupation = as.factor(X$occupation)
    X_test = as.matrix(X_test)
    X_test = cbind(1,X_test)
    n_test = nrow(X_test)
    f1 = train_data$user_id
    f2 = train_data$item_id
    
    f1_test = test_data$user_id
    f2_test = test_data$item_id
    
    
    range_f1_train=sort(unique(f1))
    range_f2_train=sort(unique(f2))
    
    train_R=length(range_f1_train)
    train_C=length(range_f2_train)
    
    
    bij_f1=list()
    bij_f2=list()
    
    for(i in 1:(train_R)){
    bij_f1[[range_f1_train[i]]]=i
    
    }
    for(i in 1:(train_C)){
    bij_f2[[range_f2_train[i]]]=i
    
    }
    
    
    
    for(i in 1:N){
    f1[i]=bij_f1[[f1[i]]]
    }
    
    for(i in 1:N){
    f2[i]=bij_f2[[f2[i]]]
    }
    
    R = max(f1)
    C = max(f2)
    
    
    for(i in 1:n_test){
    if(f1_test[i] <= length(bij_f1)){
      if(!(is.null(bij_f1[[f1_test[i]]]))){
        f1_test[i]=bij_f1[[f1_test[i]]]
      }else{
        f1_test[i] = 0
      }
    }else{
      f1_test[i] = 0}
    }
    for(i in 1:n_test){
    if(f2_test[i] <= length(bij_f2)){
      if(!(is.null(bij_f2[[f2_test[i]]]))){
        f2_test[i] = bij_f2[[f2_test[i]]]
      }else{
        f2_test[i] = 0
      }
    }else{
      f2_test[i] = 0}
    }
    
    sel = which(f1_test>0 & f2_test>0)
    X_test = X_test[sel,]
    y_test = y_test[sel]
    f1_test = f1_test[sel]
    f2_test = f2_test[sel]
    X_0 = matrix(rep(1,N), ncol = 1)
    n_test = length(y_test)
    X_test_0 = matrix(rep(1,n_test), ncol = 1)
    
    
    beta = lm(y~X-1)$coef
  
    Y = y - X%*% beta
    Y_test = y_test - X_test%*% beta
    Y = as.numeric(Y)
    incomplete_matrix = Incomplete(i=f1,j=f2,x=Y)
    set.seed(1)
    xsc=biScale(incomplete_matrix,col.scale=FALSE,row.scale=FALSE, maxit = 500)
    
    
    set.seed(1)
    fit=softImpute(xsc,rank.max= K, maxit = 5, lambda = 0)
    results = impute(fit, f1_test, f2_test)
    prediction_error = c(mean((results - Y_test)^2))
    
    max_iter = 400
    u = as.integer(max_iter/5)
    
    for(i in 1:u){
      fit = softImpute(xsc, rank.max = K, warm.start = fit, maxit = 5, lambda = 0)
      results = impute(fit, f1_test, f2_test)
      prediction_error = c(prediction_error, mean((results - Y_test)^2))
    }
    saveRDS(prediction_error, sprintf("softImpute_pred_error_lambda_0_K_%d_t_%d_occupation_%d.RDS",K,t, occupation))
    
    
    set.seed(1)
    fit=softImpute(xsc,rank.max= K, maxit = 5)
    results = impute(fit, f1_test, f2_test)
    prediction_error = c(mean((results - Y_test)^2))
    
    max_iter = 400
    u = as.integer(max_iter/5)
    
    for(i in 1:u){
      fit = softImpute(xsc, rank.max = K, warm.start = fit, maxit = 5)
      results = impute(fit, f1_test, f2_test)
      prediction_error = c(prediction_error, mean((results - Y_test)^2))
    }
    
    saveRDS(prediction_error, sprintf("softImpute_pred_error_K_%d_t_%d_occupation_%d.RDS",K,t, occupation))
    
    set.seed(1)
    fit_1 = factor_analysis_without_intercepts(y,f1,f2,K,X = X, initial_V = NULL, initial_beta = NULL, initial_sigma2e = NULL, club_V = FALSE, club_U = FALSE, max_iter = 200, conv_thres = 10^(-5), test_data = TRUE, n_test = n_test, X_test = X_test, y_test = y_test, f1_test = f1_test, f2_test = f2_test)
    saveRDS(fit_1, sprintf("factor_analysis_without_intercepts_K_%d_t_%d_occupation_%d",K, t, occupation))
    
    
    set.seed(1)
    fit_1 = factor_analysis_with_client_intercepts(y,f1,f2,K,X = X, initial_V = NULL, initial_beta = NULL, initial_sigma2e = NULL, club_U = FALSE, max_iter = 200, conv_thres = 10^(-5), test_data = TRUE, n_test = n_test, X_test = X_test, y_test = y_test, f1_test = f1_test, f2_test = f2_test)
    saveRDS(fit_1, sprintf("factor_analysis_with_client_intercepts_K_%d_t_%d_occupation_%d",K, t, occupation))
    
    
    set.seed(1)
    fit_1 = factor_analysis_with_client_item_intercept(y,f1,f2,K,X = X, initial_V = NULL, initial_beta = NULL, initial_sigma2e = NULL, club_U = FALSE, max_iter = 200, conv_thres = 10^(-5), test_data = TRUE, n_test = n_test, X_test = X_test, y_test = y_test, f1_test = f1_test, f2_test = f2_test)
    saveRDS(fit_1, sprintf("factor_analysis_with_client_item_intercept_K_%d_t_%d_occupation_%d",K, t, occupation))


    if(K==1){
    set.seed(1)
    fit = factor_analysis_with_just_client_intercept(y,f1,f2,X = X, initial_v = NULL, initial_beta = NULL, initial_sigma2e = NULL,  max_iter = 200, conv_thres = 10^(-5), test_data = TRUE, n_test = n_test , X_test = X_test, y_test = y_test, f1_test = f1_test, f2_test = f2_test)
    saveRDS(fit, sprintf("just_client_intercept_occupation_%d_t_%d",occupation, t))

    set.seed(1)
    fit = factor_analysis_with_just_client_item_intercepts(y,f1,f2,X = X, initial_beta = NULL, initial_sigma2e = NULL,  max_iter = 200, conv_thres = 10^(-5), test_data = TRUE, n_test = n_test , X_test = X_test, y_test = y_test, f1_test = f1_test, f2_test = f2_test) 
    saveRDS(fit, sprintf("just_client_item_intercept_occupation_%d_t_%d",occupation, t))
    
    }
        
    if(occupation==1){
    
    incomplete_matrix = Incomplete(i=f1,j=f2,x=y)
    xsc=biScale(incomplete_matrix,col.scale=FALSE,row.scale=FALSE, maxit = 500)
    set.seed(1)
    fit=softImpute(xsc,rank.max= K, maxit = 5)
    results = impute(fit, f1_test, f2_test)
    prediction_error = c(mean((results - y_test)^2))
    
    max_iter = 400
    u = as.integer(max_iter/5)
    
    for(i in 1:u){
      fit = softImpute(xsc, rank.max = K, warm.start = fit, maxit = 5)
      results = impute(fit, f1_test, f2_test)
      prediction_error = c(prediction_error, mean((results - y_test)^2))
    }
    
    saveRDS(prediction_error, sprintf("softImpute_pred_error_K_%d_t_%d.RDS",K,t))
   
    set.seed(1)
    fit=softImpute(xsc,rank.max= K, maxit = 5, lambda = 0)
    results = impute(fit, f1_test, f2_test)
    prediction_error = c(mean((results - y_test)^2))
    
    max_iter = 400
    u = as.integer(max_iter/5)
    
    for(i in 1:u){
      fit = softImpute(xsc, rank.max = K, warm.start = fit, maxit = 5, lambda = 0)
      results = impute(fit, f1_test, f2_test)
      prediction_error = c(prediction_error, mean((results - y_test)^2))
    }
    
    
    saveRDS(prediction_error, sprintf("softImpute_pred_error_lam_0_K_%d_t_%d.RDS",K,t))    
    

    
    
     if(K==1){
        set.seed(1)
        fit = factor_analysis_with_just_client_intercept(y,f1,f2,X = X_0, initial_v = NULL, initial_beta = NULL, initial_sigma2e = NULL,  max_iter = 200, conv_thres = 10^(-5), test_data = TRUE, n_test = n_test , X_test = X_test_0, y_test = y_test, f1_test = f1_test, f2_test = f2_test)
        saveRDS(fit, sprintf("just_client_intercept_occupation_X0_%d_t_%d",occupation, t))
        
        
        set.seed(1)
        fit = factor_analysis_with_just_client_item_intercepts(y,f1,f2,X = X_0, initial_beta = NULL, initial_sigma2e = NULL,  max_iter = 200, conv_thres = 10^(-5), test_data = TRUE, n_test = n_test , X_test = X_test_0, y_test = y_test, f1_test = f1_test, f2_test = f2_test)
        saveRDS(fit, sprintf("just_client_item_intercept_occupation_X0_%d_t_%d",occupation, t))
        
    }
  
  
  }
  

  
  

}
