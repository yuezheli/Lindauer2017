rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)

# read in the model
mod <- mread("../model/Lindauer_mus") 

# fix all random effect in the model
modmus = mod %>% zero_re() %>%
        init(A6 = 49800) %>%
        init(A10 = 7e-6) %>%
        init(A11 = 170)

MW = 149000 # antibody molecular weight

mouseweight = 20 # mouse weight = 20g

dose = 10 * mouseweight

e_10 <- c(
  ev(amt = dose*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose*1000/(MW) , cmt = 'A1', time = 7 * 24), 
  ev(amt = dose*1000/(MW) , cmt = 'A1', time = 14 * 24)
)


samplingtime_mouse = c(1,   3,   6,   16,  24,  72,  120, 
                       168, 169, 171, 174, 184, 192, 240, 288,
                       336, 337, 339, 342, 352, 360, 408, 456)

sim_10 <- modmus %>% ev(e_10) %>%
  mrgsim(delta = 1, end = 480) %>%
  as_tibble() %>% filter(time %in% samplingtime_mouse) %>% 
  select(time, CP, TV) 

# see data; this code is not necessary for the purpose of this script
par(mfrow = c(1,2))
plot(sim_10$time, sim_10$CP)
plot(sim_10$time, sim_10$TV)

# save result for comparison
write.table(sim_10, '../data/mrgresult.csv', row.names = FALSE)
