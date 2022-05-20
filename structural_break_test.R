##########################################################################
##########################################################################
# Project:
# Script purpose: test structual break in vertebrae data of axolotl in wt and mutant
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri May 20 09:39:58 2022
##########################################################################
##########################################################################


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

chow_test = function(tt = c(1:21), wt, ko, rm.outliers = TRUE)
{
  # tt = aa[, 1]; wt = aa[, c(2:11)]; ko = aa[, 12:21]
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
  
  ##########################################
  # the chow test here is modified based on the 
  # https://github.com/lengfei5/GachonProt/blob/main/scripts/proteomics_functions.R
  # and the paper 
  # https://www.researchgate.net/publication/336476565_Performance_of_Methods_Determining_Structural_Break_in_Linear_Regression_Models
  ##########################################
  ## fitting the pool of wt and ko
  fit = lm(wt.ko ~ tt.wt.ko)
  ee = sum(fit$residuals^2)
  
  ## fitting wt
  fit1 = lm(data.wt ~ tt.wt)
  ee1 = sum(fit1$residuals^2)
  
  # fit ko 
  fit2 = lm(data.ko ~ tt.ko)
  ee2 = sum(fit2$residuals^2)
  
  k = 2 # nb of parameteres
  n1 = length(data.wt)
  n2 = length(data.ko)
  n = length(wt.ko)
  
  F1 = (ee-ee1)/n2/(ee1/(n1-k)) # F_CPT 
  pval.F1 = pf(F1, n2, (n1-k), lower.tail = FALSE, log.p = FALSE)
  
  F2 = (ee-ee1-ee2)/k/((ee1+ee2)/(n - 2*k)) # F_CBT
  pval.F2 = pf(F2, k, (n - 2*k), lower.tail = FALSE, log.p = FALSE)
  
 
  
  plot(tt.wt, data.wt, ylim = range(wt.ko), col = 'darkblue')
  abline(fit1$coefficients, col = 'darkblue', lwd = 2.0)
  points(tt.ko, data.ko, col = 'red')
  abline(fit$coefficients, col = 'red', lwd = 2.0)
  
  #print(pval.wt.ko)
  return(pval.F1)
  
}


##########################################
# main 
##########################################
require(openxlsx)
library(tidyr)
library(dplyr)
require(ggplot2)
library(stringr)

aa = openxlsx::read.xlsx('vertebrae_length.xlsx', sheet = 1)
colnames(aa)[1] = 'vnb'
colnames(aa)[c(2:11)] = paste0('wt_ax', c(1:10))
colnames(aa)[c(12:21)] = paste0('ko_ax', c(1:10))

as_tibble(aa) %>% gather(ax, lengths, 2:21) %>%
  separate(ax, c('condition', 'ax')) %>%
  ggplot(aes(x = vnb, y = lengths,  group=condition, color = condition)) + 
  geom_point() +
  #geom_boxplot()
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 10)) +
  scale_fill_manual(values=c('#999999','#E69F00')) +
  labs( x = 'vertebrae lengths', y = 'lengths')
  
aa = as.matrix(aa)
chow_test(tt = aa[, 1], wt = aa[, c(2:11)], ko = aa[, 12:21])


