rm(list = ls())

set.seed(0)

# load libraries
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
##------------------- set up dosing events -------------------##
# human dosing regiment: dose every 3 weeks, dose for 6 month
humandosing <- function(Dose)
{
  humanweight = 71 * 1000 # assume human bodyweight = 71kg
  MW	= 149000	# g/mol ; Molecular weight of antibody
  dose = Dose * humanweight
  e  <- c(
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 0 * 24), 
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 3 * 7 * 24), 
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 6 * 7 * 24), 
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 9 * 7 * 24), 
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 12 * 7 * 24), 
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 15 * 7 * 24), 
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 18 * 7 * 24), 
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 21 * 7 * 24), 
    ev(amt = dose*1000/(MW) , cmt = 'A1', time = 24 * 7 * 24)
  )
  return(e)
}

##-------------------------- prepare the model --------------------------##
mod <- mread("../model/Lindauer_homo") 


DiameterReduction <- function(Dose = 0, tvl0 = 0.0017, tvsltg = 2.886e-7, N = 300, modhomo = mod)
{
  modA <- modhomo %>% param(TVL0 = tvl0) %>% param(TVSLtg = tvsltg) # change growth parameter
  
  # at dosing event
  e <- humandosing(Dose)
  
  # simulation
  sim <- modA %>% ev(e) %>%
    mrgsim(end = 26 * 7 * 24, nid = N) %>%
    as_tibble() %>% select(time, radiuschange) %>% group_by(time) %>% 
    summarize(lo=quantile(radiuschange, 0.1), hi=quantile(radiuschange, 0.9), med = quantile(radiuschange, 0.5) )
  
  ratio_high = 1-tail(sim$lo, n=1)/head(sim$lo, n=1)
  ratio_low = 1-tail(sim$hi, n=1)/head(sim$hi, n=1)
  ratio_med = 1-tail(sim$med, n=1)/head(sim$med, n=1)
  
  ratio = c(ratio_low, ratio_high, ratio_med)
  
  return(ratio) # return changes in tumor diameter
  
}

doses = c(0, 0.1, 0.2, 0.3, 0.75, 1, 2, 3, 5, 10)

diameterratio_slowgrowth = data.frame(matrix(rep(NA, length(doses) * 3), ncol = 3))
diameterratio_slowallo = data.frame(matrix(rep(NA, length(doses) * 3), ncol = 3))
diameterratio_medgrowth = data.frame(matrix(rep(NA, length(doses) * 3), ncol = 3))
diameterratio_medallo = data.frame(matrix(rep(NA, length(doses) * 3), ncol = 3))
diameterratio_fastgrowth = data.frame(matrix(rep(NA, length(doses) * 3), ncol = 3))
diameterratio_fastallo = data.frame(matrix(rep(NA, length(doses) * 3), ncol = 3))

for(i in 1:length(doses))
{
  diameterratio_fastallo[i, ]  =  DiameterReduction(Dose = doses[i], tvl0 = 0.0088, tvsltg = 2.575e-6)   # fast growth & allometric scaling
  diameterratio_fastgrowth[i, ]  =  DiameterReduction(Dose = doses[i], tvl0 = 0.0088, tvsltg = 1.542e-6) # fast growth & growth-proportional scaling
  diameterratio_medallo[i, ]  =  DiameterReduction(Dose = doses[i], tvl0 = 0.0036, tvsltg = 2.575e-6)    # medium growth & allometric scaling
  diameterratio_medgrowth[i, ]  =  DiameterReduction(Dose = doses[i], tvl0 = 0.0036, tvsltg = 6.247e-7)  # medium growth & growth-proportional scaling
  diameterratio_slowallo[i, ]  =  DiameterReduction(Dose = doses[i], tvl0 = 0.0017, tvsltg = 2.575e-6)   # slow growth & allometric scaling
  diameterratio_slowgrowth[i, ]  =  DiameterReduction(Dose = doses[i]) # slow growth & growth scaling
}

# combine all the data into 1 data frame

dt1 <- data.frame(doses = c(doses,doses, doses, doses, doses, doses),
                  ymin = c(diameterratio_slowgrowth$X1,
                           diameterratio_slowallo$X1,
                           diameterratio_medgrowth$X1,
                           diameterratio_medallo$X1,
                           diameterratio_fastgrowth$X1,
                           diameterratio_fastallo$X1),
                  ymax = c(diameterratio_slowgrowth$X2,
                           diameterratio_slowallo$X2,
                           diameterratio_medgrowth$X2,
                           diameterratio_medallo$X2,
                           diameterratio_fastgrowth$X2,
                           diameterratio_fastallo$X2),
                  grp =  rep(c('Growth: slow | Scaling: growth','Growth: slow | Scaling: allo', 'Growth: med | Scaling: growth', 
                               'Growth: med | Scaling: allo', 'Growth: fast | Scaling: growth', 'Growth: fast | Scaling: allo'),each = length(doses)))

dt2 = data.frame(doses = c(doses,doses, doses, doses, doses, doses),
                 ymed = c(diameterratio_slowgrowth$X3,
                          diameterratio_slowallo$X3,
                          diameterratio_medgrowth$X3,
                          diameterratio_medallo$X3,
                          diameterratio_fastgrowth$X3,
                          diameterratio_fastallo$X3),
                 grp = rep(c('Growth: slow | Scaling: growth','Growth: slow | Scaling: allo', 'Growth: med | Scaling: growth', 
                             'Growth: med | Scaling: allo', 'Growth: fast | Scaling: growth', 'Growth: fast | Scaling: allo'), each = length(doses))
)

disp <- ggplot() + 
  geom_ribbon(data = dt1,aes(x = doses,ymin = ymin, ymax = ymax,fill = grp), alpha = 0.7) + 
  geom_line(data = dt2, aes(x = doses, y = ymed, linetype = grp), show.legend=T) + 
  labs(y = 'tumor diameter reduction', x = 'doses (mg/kg)') +
  scale_y_continuous(breaks = seq(-0.2, 1, by = 0.2), limits = c(-0.2, 0.85)) + 
  scale_x_continuous(breaks = seq(0, 10, by = 2) ) + theme_bw() + 
  theme(panel.grid.minor = element_blank(), legend.title = element_blank(), 
        legend.text=element_text(size=20), text = element_text(size=20)) 

png('../img/Fig3a.png', width = 1000, height = 500)
print(disp)
dev.off()