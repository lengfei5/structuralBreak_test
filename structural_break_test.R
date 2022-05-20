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
# main 
##########################################
rm(list =ls())

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

source('functions_structural_break.R')

structural_break_test(tt = aa[, 1], wt = aa[, c(2:11)], ko = aa[, 12:21], rm.outliers = TRUE, method = 'chow_test')

structural_break_test(tt = aa[, 1], wt = aa[, c(2:11)], ko = aa[, 12:21], rm.outliers = FALSE, method = 'hpt_test')

structural_break_test(tt = aa[, 1], wt = aa[, c(2:11)], ko = aa[, 12:21], rm.outliers = FALSE, method = 'mz_test')
