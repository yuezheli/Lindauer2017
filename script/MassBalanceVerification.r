# this script aims to verify mass balance of mAb in the mouse model

rm(list = ls())
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

W0 = 170 # initial tumor size (mm^3)
MW = 149000 # antibody molecular weight

# prepare a fixed model for simulation
modo <- mread("../model/Lindauer_mus") %>% init(A10 = 7.45e-6) %>% init(A6 = 49800) %>% init(A11 = W0)  %>% zero_re()  

# dosing regiment
sunev <- function(doseperweight = 5)
{
  MW = 149000
  mouseweight = 20 # unit: g
  dose = doseperweight * mouseweight
  
  dosingevent <- c(
#    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 0 * 24),
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 5 * 24),
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 9 * 24),
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 13 * 24)
  )
  
  return(dosingevent)
}


modo %>% 
  ev( sunev( 5 ) )  %>%
  mrgsim(delta = 1, end = 1080) %>%
  as_tibble() %>% dplyr::filter(row_number() != 1) %>%
  as.data.frame() %>% 
  plot_sims(.f = totalmAb ~ time) %>% print()
