# The aim of this global analysis document is to focus on PKPD/ OD parameters. 
# PK parameters should either have no variation, or should have a much smaller variation compared to the PD parameters. 
# global sensitivity analysis code adapted from https://github.com/metrumresearchgroup/pbpk-qsp-mrgsolve/blob/master/docs/global_sensitivity_analysis.md

rm(list = ls())
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(sensitivity)
library(PKPDmisc)
library(mrgsim.parallel)

# set up parallel computation
number_of_cores = 144
options(mc.cores = number_of_cores)

# set up simulation parameters for sobol2007 function
simulationboot = 1e3 # number of simulation time in sobol2007
sampleperparam = 150 # sampling size per parameter, for sobol2007
set.seed(88771) # seed adopted from the GitHub page mentioned at the beginning

##----------------------- prepare model for sensitivity analysis -----------------------##

W0 = 170 # initial tumor size (mm^3)
MW = 149000 # antibody molecular weight

# read in the model
mod <- mread("../model/Lindauer_mus") 

mod  = mod %>% update(end = 480, delta = 1) %>% 
  zero_re()# fix the randomness in the model


# dosing regiment
sunev <- function(doseperweight = 5)
{
  MW = 149000 # molecular weight of the antibody
  mouseweight = 20 # unit: g
  dose = doseperweight * mouseweight
  
  dosingevent <- c(
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 0 * 24),
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 5 * 24),
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 9 * 24),
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 13 * 24)
  )
  
  return(dosingevent)
}

##----------------------- Global sensitivity analysis -----------------------##

setpk <- c('TVW0', 'TVL0', 'TVCL')
setpd <- c('TVPLQ', 'TVCLup','TVKdeg' ,'TVFR', 'V_blood', 'Tmulti', 'EMAXTP', 'EC50TP', 'KdegPD1', 'TVL1')
setkil <- c('TVSLtg')
# binding affinity of mAb to FcRn and PD-1 are not analyzed; the assumption is that once an antibody is given, these 2 parameters should be fixed; 
# it is also assumed that blood T cells concentration is a constant

# create sampling method
gen_samples <- function(n, l, which = names(l), 
                        factor = c(0.01,100), N = NULL) { # NEW ARGUMENT: N, the absolute sampling number per parameter
  
  vars <- tidyselect::vars_select(names(l), !!(enquo(which)))
  
  l <- as.list(l)[vars]
  
  l <- lapply(l, function(x) x*factor )
  
  if(is.numeric(N)) { # NEW
    n <- N*2  
  } else {
    n <- length(l)*n*2
  }
  
  df <- as.data.frame(l)
  
  len <- length(df)
  
  X <- matrix(ncol=len, nrow=n)
  
  colnames(X) <- names(df)
  
  Y <- X
  
  for(i in seq(len)){
    r <- exp(runif(n, log(df[1,i]), log(df[2,i])))
    X[,i] <- r
    r <- exp(runif(n, log(df[1,i]), log(df[2,i])))
    Y[,i] <- r
  }
  
  return(list(x1 = as.data.frame(X), x2 = as.data.frame(Y)))
}

N <- sampleperparam * length(c(setpk,setpd, setkil))

# set different variation in each group
sets <- list(
  list(which = setpk, factor = c(0.9, 1.1), N = N, l = param(mod)), 
  list(which = setpd, factor = c(0.01,100), N = N, l = param(mod)), 
  list(which = setkil, factor = c(0.99,1.01), N = N, l = param(mod))
)

# The following lines of code are commented out for efficiency in reproducing
# the manuscript figures. Instead of performing the Sobol sensitivity analysis
# we load the results of a previous analysis.
#samps <- map(sets, do.call, what = gen_samples)
#samp <- map(c(1,2), ~ map_dfc(samps, .x))

##----------------------- Global sensitivity analysis using tumor volume change as the readout -----------------------##

# tumor volume change
startendvar <- function(y)
{
  x = tail(y, n=1)/head(y, n=1)
  
  if(is.na(x))
  {
    return(0)
  }
  else{
    return(x)
  }
}

batch_tumorvol <- function(x) {
  mod %>% 
    idata_set(x) %>%
    ev(sunev()) %>%
    mrgsim(obsonly = TRUE) %>% 
    group_by(ID) %>% 
    summarise(tvvar = startendvar(TV)) %>% 
    pull(tvvar)
}

# # warning: the following line takes several minutes to run
# tvglobal = sobol2007(batch_tumorvol, X1=samp[[1]], X2=samp[[2]], nboot=simulationboot)
# plot(tvglobal)
# lines(c(0,31), c(0.05, 0.05), type="l", lty = 2)
# # save a minimal set of the Sobol analysis data
# # that is sufficient for visualization
#sobolx <- list(S = tvglobal$S , T = tvglobal$T)
#saveRDS(sobolx, file = "../data/SobolData.rds")
# load previously saved Sobol sensitivity analysis data (minimal dataset)
x <- readRDS("../data/SobolData.rds")
### plot
globSens <- tibble(Parameter = c(setpk,setpd, setkil),
                   main = x$S$original,
                   main_lo = x$S$'min. c.i.',
                   main_hi = x$S$'max. c.i.',
                   total = x$T$original,
                   total_lo = x$T$'min. c.i.',
                   total_hi = x$T$'max. c.i.') %>%
  gather(effect, Index, -Parameter, -main_lo, -main_hi, -total_lo, -total_hi) %>%
  mutate(lo = ifelse(effect == "main", main_lo, total_lo),
         hi = ifelse(effect == "main", main_hi, total_hi),
         Effect = factor(effect))

fig_sens_glob <- ggplot(data=globSens, aes(x=Parameter, y=Index, group=Effect, col=Effect)) +
  geom_point(position = position_dodge(width=0.3)) +
  geom_errorbar(aes(ymax=hi, ymin=lo), position=position_dodge(width=0.3), width=0) +
  theme_bw() +
  geom_hline(aes(yintercept=0.05), lty=2) +
  theme(legend.title = element_text(size = 15), legend.text=element_text(size=12)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12))

png('../img/globalsens_tv.png', width = 800, height = 300)
print(fig_sens_glob)
dev.off()


##----------------------- Global sensitivity analysis using endosomal mAb exposure as the readout -----------------------##
# please note this section of code can be extremely computational expensive, thus the codes are not run in the current situation

# exposure of endosomal antibody concentration
batch_endo_pc <- function(x) {
  x <- mutate(x, ID = row_number())
  mc_mrgsim_ei(mod, sunev(), x, nchunk = number_of_cores, .parallel = TRUE) %>% 
    na.omit() %>%
    group_by(ID) %>% 
    summarise(mAb = auc_partial(time, ces)) %>% 
    pull(mAb)
}

Sys.time()
#png('../img/pd_globalsens_endo.png', width = 600, height = 300)
#sobol2007(batch_endo_pc, X1=samp[[1]], X2=samp[[2]], nboot=simulationboot) %>% plot()
#lines(c(0,31), c(0.05, 0.05), type="l", lty = 2)
#dev.off()
Sys.time()


##----------------------- Global sensitivity analysis using tumor receptor occupancy as the readout -----------------------##

batch_ROt_pc <- function(x) {
  x <- mutate(x, ID = row_number())
  mc_mrgsim_ei(mod, sunev(), x, nchunk = number_of_cores, .parallel = TRUE) %>% 
    group_by(ID) %>% 
    summarise(rot = mean(ROt, na.rm = TRUE)) %>% 
    pull(rot)
}

batch_ROt <- function(x) {
  mrgsim_ei(mod, sunev(), x) %>% 
    group_by(ID) %>% 
    summarise(rot = mean(ROt, na.rm = TRUE)) %>% 
    pull(rot)
}

# note the code is not run here, because it requires the sample (variable samp) generated at the beginning
# without that line of code, this cannot work

#png('../img/pd_globalsens_ROt.png', width = 600, height = 300)
Sys.time()
#sobol2007(batch_ROt_pc, X1=samp[[1]], X2=samp[[2]], nboot=simulationboot) %>% plot()
#lines(c(0,31), c(0.05, 0.05), type="l", lty = 2)
Sys.time()
#dev.off()

