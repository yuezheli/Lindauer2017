# local sensitivity analysis for Lindauer model (mouse model)

rm(list = ls())
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(sensitivity)

# read in the model
mod <- mread("../model/Lindauer_mus") 

# fix the randomness in the model
mod  = mod %>% zero_re()

# set up dosing event
sunev <- function(doseperweight = 5, mouseweight = 20)
{
  MW = 149000
  # mouseweight = 20 # unit: g
  dose = doseperweight * mouseweight
  
  dosingevent <- c(
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 0 * 24),
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 5 * 24),
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 9 * 24),
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 13 * 24)
  )
  
  return(dosingevent)
}


# set up function to track changes in tumor volume

SimulationTrial <- function(paranumber, foldofchange = 5, event= sunev(doseperweight, mouseweight), tmpmod = mod)
{
  allparameters = param(tmpmod)
  new_values = allparameters[paranumber]
  new_values[[1]] = unlist(allparameters[paranumber])[[1]] * foldofchange
  
  modB <- param(tmpmod, new_values)
  
  out <- tryCatch(
    {
      modB %>% ev(event)  %>%
        mrgsim(delta = 1, end = 480) %>%
        as_tibble() 
    },
    error = function(cond)
    {
      print( names(allparameters[paranumber]) )
      return(NA)
    },
    warning=function(cond) {
      
      return(NA) 
      
    },
    finally={ 
      sim = modB %>% ev(event)  %>%
        mrgsim(delta = 1, end = 20 * 24) %>%
        as_tibble() %>% data.frame()
      return(  log10( tail(sim$TV, n=1)/head(sim$TV, n=1) )  ) 
    }
  )
}


M = length(param(mod)) # number of parameters in the model

outcome = rep(NA, 2*M) %>% matrix(ncol = M) 
varnames = rep(NA, M)

# calculate the baseline change
baselinechange = SimulationTrial(1, foldofchange = 1, event= sunev(doseperweight = 0.1, mouseweight = 20))

for( i in 1:30)
{
  varnames[i] = names(param(mod)[i])
  outcome[1,i] = SimulationTrial(i, foldofchange = 0.2, event= sunev(doseperweight = 0.1, mouseweight = 20))
  outcome[2,i] = SimulationTrial(i, foldofchange = 5, event= sunev(doseperweight = 0.1, mouseweight = 20))
}


rownames(outcome) <- c('20%','500%')  
colnames(outcome) <- varnames  

# sort the change in tumor size based on 5* 
outcome = outcome[, order(outcome['20%', ], na.last = NA)]

outcome5 = outcome[,1:6]

# display the normalized result
barplot(outcome5 - baselinechange, horiz = T, xaxt='n', ylab = ' ',
        xlab = 'log of tumor size change (normalized)',
        beside=T, col=c('springgreen','indianred2'), 
        legend = TRUE, args.legend = list(x = "topright"), 
        xlim = c(-1,1), 
        las = 1, cex.names = 0.7)
x <- seq(-1,1, length=5)     
axis(1, at=pretty(x),  lab=paste0(pretty(x) ), las=TRUE)
# print(colnames(outcome))
#options(scipen=999)
#print(t(outcome5-baselinechange), digits = 3)