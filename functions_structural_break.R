##########################################
# functions of statistical test 
##########################################
index.outliers = function(data.xx)
{
  c = 1.5
  Q1 = quantile(data.xx, 0.25,type=5)
  Q3 = quantile(data.xx, 0.75, type=5)
  IQD = Q3 - Q1
  lower = Q1 - c*IQD
  upper = Q3 + c*IQD
  index = which(data.xx<lower|data.xx>upper)
}


structural_break_test = function(tt = c(1:21), wt, ko, rm.outliers = FALSE, method = 'chow.test', 
                                 nBoot = 10000)
{
  # tt = aa[, 1]; wt = aa[, c(2:11)]; ko = aa[, 12:21]; rm.outliers = FALSE; nBoot = 10000;
  tt  = as.numeric(tt);
  wt = as.matrix(wt)
  ko = as.matrix(ko)
  
  data.wt = c()
  tt.wt = c()
  for(n in 1:ncol(wt))
  {
    tt.wt = c(tt.wt, tt[which(!is.na(wt[, n]))])
    data.wt = c(data.wt, as.numeric(wt[which(!is.na(wt[,n])) ,n]))
  }
  cat(length(tt.wt), ' tt -- ', length(data.wt), ' wt data points \n')
  
  if(rm.outliers){
    ii_out = index.outliers(data.wt)
    if(length(ii_out)>0){
      tt.wt = tt.wt[-ii_out]
      data.wt = data.wt[-ii_out]
    }
    cat(length(tt.wt), ' tt -- ', length(data.wt), ' wt data points after outlier removal \n')
  }
  
  data.ko = c()
  tt.ko = c()
  for(n in 1:ncol(ko))
  {
    tt.ko = c(tt.ko, tt[which(!is.na(ko[, n]))])
    data.ko = c(data.ko, as.numeric(ko[which(!is.na(ko[,n])) ,n]))
  }
  cat(length(tt.ko), ' tt -- ', length(data.ko), ' ko data points \n')
  
  if(rm.outliers){
    ii_out = index.outliers(data.ko)
    if(length(ii_out)>0){
      tt.ko = tt.ko[-ii_out]
      data.ko = data.ko[-ii_out]
    }
    cat(length(tt.ko), ' tt -- ', length(data.ko), ' ko data points after outlier removal \n')
  }
  
  wt.ko = c(data.wt, data.ko)
  tt.wt.ko = c(tt.wt, tt.ko)
  cat(length(tt.wt.ko), ' tt -- ', length(wt.ko), ' pool of wt and ko \n' )
  
  if(method == 'chow_test'){
    ##########################################
    # the chow test here is modified based on the 
    # https://github.com/lengfei5/GachonProt/blob/main/scripts/proteomics_functions.R
    # and the paper 
    # https://www.researchgate.net/publication/336476565_Performance_of_Methods_Determining_Structural_Break_in_Linear_Regression_Models
    ##########################################
    cat('---------------------------------------\n')
    cat('chow_test for coefficients is performed \n ')
    
    k = 2 # nb of parameteres
    n1 = length(data.wt)
    n2 = length(data.ko)
    n = length(wt.ko)
    
    ## fitting the pool of wt and ko
    fit = lm(wt.ko ~ tt.wt.ko)
    ee = sum(fit$residuals^2)
    
    ## fitting wt
    fit1 = lm(data.wt ~ tt.wt)
    ee1 = sum(fit1$residuals^2)
    
    # fit ko 
    fit2 = lm(data.ko ~ tt.ko)
    ee2 = sum(fit2$residuals^2)
  
    cat('F_CPT in chow test is used \n')
    
    cat('---------------------------------------\n')
    
    F1 = (ee-ee1)/n2/(ee1/(n1-k)) # F_CPT 
    pval.F1 = pf(F1, n2, (n1-k), lower.tail = FALSE, log.p = FALSE)
    
    #cat('F_DBT is used \n')
    F2 = (ee-ee1-ee2)/k/((ee1+ee2)/(n - 2*k)) # F_CBT
    pval.F2 = pf(F2, k, (n - 2*k), lower.tail = FALSE, log.p = FALSE)
    
    plot(tt.wt, data.wt, ylim = range(wt.ko), col = 'darkblue')
    abline(fit1$coefficients, col = 'darkblue', lwd = 2.0)
    points(tt.ko, data.ko, col = 'red')
    abline(fit$coefficients, col = 'red', lwd = 2.0)
    #print(pval.wt.ko)
    #return(pval.F1)
    pval = pval.F2
    
  }
  
  if(method == 'hpt_test'){
    
    cat('---------------------------------------\n')
    cat('hpt_test for variance is performed \n ')
    
    k = 2 # nb of parameteres
    n1 = length(data.wt)
    n2 = length(data.ko)
    n = length(wt.ko)
    
    ## fitting wt
    fit1 = lm(data.wt ~ tt.wt)
    ee1 = sum(fit1$residuals^2)
    ss1 = sqrt(ee1/n1)
    
    pred2 = predict(fit1, data.frame(tt.wt = tt.ko))
    
    ee2 = sum((pred2 - data.ko)^2)
    ss2 = sqrt(ee2/n2)
    
    barplot(c(ss1, ss2), names.arg = c('wt', 'ko'), col = c('darkblue', 'darkred'), ylab = 'estimated sigma')
    
    hpt = ee2/(ee1/(n1-k))
    pval = pchisq(hpt, n2, lower.tail = FALSE, log.p = FALSE)
    
    cat('---------------------------------------\n')
    
  }
  
  ##########################################
  # mz test from the same paper mentioned above
  # some original code from 
  # https://github.com/myaseen208/SupMZ/blob/master/R/supmz.R
  ##########################################
  if(method == 'mz_test')
  {
    cat('---------------------------------------\n')
    cat('mz_test for coefficients and variance is performed \n ')
    
    k = 2 # nb of parameteres
    n1 = length(data.wt)
    n2 = length(data.ko)
    n = length(wt.ko)
    
    ## fitting the pool of wt and ko
    fit = lm(wt.ko ~ tt.wt.ko)
    ss =  sigma(fit)^2
    
    
    ## fitting wt
    fit1 = lm(data.wt ~ tt.wt)
    ss1 = sigma(fit1)^2
    
    # fit ko 
    fit2 = lm(data.ko ~ tt.ko)
    ss2 = sigma(fit2)^2
    mz = (n - k) * log(ss) - ((n1 - k) * log(ss1) + (n2 - k)*log(ss2))
    
    ## simulate the mz distribution under null hypothesis
    fitted <- as.numeric(fitted(fit))
    cat('--- start simulation for null hypothesis -- \n')
    # nBoot = 10000
    mzs = rep(NA, nBoot)
    set.seed(seed = 2022)
    for(i in 1:nBoot)
    {
      if(i%%1000 == 0) cat(i, '\n')
      data.sim = fitted + rnorm(n = length(wt.ko), mean = 0, sd = sqrt(ss))
      data1 = data.sim[1:n1]
      data2 = data.sim[c((n1+1):length(wt.ko))]
      
      fm0 = lm(data.sim ~ tt.wt.ko)
      ss00 = sigma(fm0)^2
      fm1 = lm(data1 ~ tt.wt)
      ss01 = sigma(fm1)^2
      fm2 = lm(data2 ~ tt.ko)
      ss02 = sigma(fm2)^2
      mzs[i] = (n-k)*log(ss00) - ((n1-k)*log(ss01) + (n2 -k)*log(ss02))
      
    }
    
    Fn = ecdf(mzs)
    pval = 1 - Fn(mz)
    
    hist(mzs, xlim = range(mzs, mz), xlab = 'MZ statistics')
    abline(v = mz, lwd = 2.0, col = 'red')
    cat('---------------------------------------\n')
    
  }
  
  cat('pval --', pval, '\n')
  return(pval)
  
}


