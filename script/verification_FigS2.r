rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

# read in the model
mod <- mread("../model/Lindauer_mus") 

W0 = 170 # initial tumor size (mm^3)
MW = 149000 # antibody molecular weight

mouseweight = 20 # mouse weight = 20g


# read in the observed PK data
data_1 <- read.csv('../data/FigS2/dose_1.csv', header = TRUE)
data_point1 <- read.csv('../data/FigS2/dose_point1.csv', header = TRUE)
data_10 <- read.csv('../data/FigS2/dose_10.csv', header = TRUE)

data_point4 <- read.csv('../data/FigS2/dose_point4.csv', header = TRUE)
data_1point4 <- read.csv('../data/FigS2/dose_1point4.csv', header = TRUE)
data_5 <- read.csv('../data/FigS2/dose_5.csv', header = TRUE)

# 6 different dosing
dose_1 = 1 * mouseweight
dose_point1 = 0.1 * mouseweight
dose_10 = 10 * mouseweight

dose_point4 = 0.4 * mouseweight
dose_1point4 = 1.4 * mouseweight
dose_5 = 5 * mouseweight
# 6 different doses and dosing regiments
e10 <- c(
  ev(amt = dose_10*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose_10*1000/(MW) , cmt = 'A1', time = 7 * 24), 
  ev(amt = dose_10*1000/(MW) , cmt = 'A1', time = 14 * 24)
)

e1 <- c(
  ev(amt = dose_1*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose_1*1000/(MW) , cmt = 'A1', time = 7 * 24), 
  ev(amt = dose_1*1000/(MW) , cmt = 'A1', time = 14 * 24)
)


e_point1 <- c(
  ev(amt = dose_point1*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose_point1*1000/(MW) , cmt = 'A1', time = 7 * 24), 
  ev(amt = dose_point1*1000/(MW) , cmt = 'A1', time = 14 * 24)
)

e_point4 <- c(
  ev(amt = dose_point4*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose_point4*1000/(MW) , cmt = 'A1', time = 4 * 24)
)

e_1point4 <- c(
  ev(amt = dose_1point4*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose_1point4*1000/(MW) , cmt = 'A1', time = 4 * 24)
)

e_5 <- c(
  ev(amt = dose_5*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose_5*1000/(MW) , cmt = 'A1', time = 4 * 24)
)


N = 100 # population size for simulation

samplingtime_mouse = c(0,   1,   3,   6,   16,  24,  72,  120, 
                       168, 169, 171, 174, 184, 192, 240, 288,
                       336, 337, 339, 342, 352, 360, 408, 456)

samplingtime_rat = c(0, 24, 96, 
                     120, 192, 288,
                     456)

## ------------- simulation with data filtered based on samling frequency in Table S1 ------------ ##
sim10 <- mod %>% ev(e10) %>%
  mrgsim(delta = 1, end = 20 * 24, nid  = N) %>% 
  as_tibble() %>% filter(time %in% samplingtime_mouse) %>% 
  select(time, CP) %>% group_by(time)  %>%
  summarize(lo=quantile(CP, 0.1), hi=quantile(CP, 0.9), med = quantile(CP, 0.5)) 

sim1 <- mod %>% ev(e1) %>%
  mrgsim(delta = 1, end = 20 * 24, nid  = N) %>%
  as_tibble() %>% filter(time %in% samplingtime_mouse) %>% 
  select(time, CP) %>% group_by(time) %>% 
  summarize(lo=quantile(CP, 0.1), hi=quantile(CP, 0.9), med = quantile(CP, 0.5))

sim_point1 <- mod %>% ev(e_point1) %>%
  mrgsim(delta = 1, end = 20 * 24, nid  = N) %>%
  as_tibble() %>% filter(time %in% samplingtime_mouse) %>% 
  select(time, CP) %>% group_by(time) %>% 
  summarize(lo=quantile(CP, 0.1), hi=quantile(CP, 0.9), med = quantile(CP, 0.5))

sim_point4 <- mod %>% ev(e_point4) %>%
  mrgsim(delta = 1, end = 20 * 24, nid  = N) %>%
  as_tibble() %>% filter(time %in% samplingtime_rat) %>% 
  select(time, CP) %>% group_by(time) %>% 
  summarize(lo=quantile(CP, 0.1), hi=quantile(CP, 0.9), med = quantile(CP, 0.5))

sim_1point4 <- mod %>% ev(e_1point4) %>%
  mrgsim(delta = 1, end = 20 * 24, nid  = N) %>%
  as_tibble() %>% filter(time %in% samplingtime_rat) %>% 
  select(time, CP) %>% group_by(time) %>% 
  summarize(lo=quantile(CP, 0.1), hi=quantile(CP, 0.9), med = quantile(CP, 0.5))

sim_5 <- mod %>% ev(e_5) %>%
  mrgsim(delta = 1, end = 20 * 24, nid  = N) %>%
  as_tibble() %>% filter(time %in% samplingtime_rat) %>% 
  select(time, CP) %>% group_by(time) %>% 
  summarize(lo=quantile(CP, 0.1), hi=quantile(CP, 0.9), med = quantile(CP, 0.5))


## ------------- plot raw simulation data against the observed data ------------ ##

plasmaconc_10 <- ggplot() + 
  #  geom_line(data = sim10, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim10, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.5) + 
  geom_point(data = data_10, aes(x = time, y = plasma_conc, color = 'observed'), size = 1) + 
  labs(y = 'plasma conc (ug/ml)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 10mg/kg, 79AEA') + xlim(0.01, 20)


plasmaconc_1 <- ggplot() + 
  #  geom_line(data = sim1, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim1, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.5) + 
  geom_point(data = data_1, aes(x = time, y = plasma_conc, color = 'observed'), size = 1) + 
  labs(y = 'plasma conc (ug/ml)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 1mg/kg, 79AEA') + xlim(0.01, 20)

plasmaconc_point1 <- ggplot() + 
  #  geom_line(data = sim_point1, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim_point1, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.5) + 
  geom_point(data = data_point1, aes(x = time, y = plasma_conc, color = 'observed'), size = 1) + 
  labs(y = 'plasma conc (ug/ml)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 0.1mg/kg, 79AEA') + xlim(0.01, 20)


plasmaconc_point4 <- ggplot() + 
  #  geom_line(data = sim_point4, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim_point4, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.5) + 
  geom_point(data = data_point4, aes(x = time, y = plasma_conc, color = 'observed'), size = 1) + 
  labs(y = 'plasma conc (ug/ml)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 0.4mg/kg, 80ADW') + xlim(0.01, 20)

plasmaconc_1point4 <- ggplot() + 
  #  geom_line(data = sim_1point4, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim_1point4, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.5) + 
  geom_point(data = data_1point4, aes(x = time, y = plasma_conc, color = 'observed'), size = 1) + 
  labs(y = 'plasma conc (ug/ml)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 1.4mg/kg, 80ADW') + xlim(0.01, 20)


plasmaconc_5 <- ggplot() + 
  #  geom_line(data = sim_5, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim_5, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.5) + 
  geom_point(data = data_5, aes(x = time, y = plasma_conc, color = 'observed'), size = 1) + 
  labs(y = 'plasma conc (ug/ml)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 5mg/kg, 80ADW') + xlim(0.01, 20)

figs2 = grid.arrange(plasmaconc_5, plasmaconc_10, plasmaconc_1, plasmaconc_1point4, 
             plasmaconc_point1, plasmaconc_point4, ncol=2)

# save the plot
ggsave(
  '../img/FigS2.png',
  plot = figs2,
  scale = 1.2,
  dpi = 300,
  limitsize = FALSE,
)
