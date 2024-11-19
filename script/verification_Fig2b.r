rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

# read in the model
mod <- mread("../model/Lindauer_mus") 


MW = 149000 # antibody molecular weight
mouseweight = 20 # mouse weight = 20g
# dosing amount
dose_point1 = 0.1 * mouseweight
dose_point4 = 0.4 * mouseweight
dose_1point4 = 1.4 * mouseweight
dose_5 = 5 * mouseweight

e_point1 <- c(
  ev(amt = dose_point1*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose_point1*1000/(MW) , cmt = 'A1', time = 5 * 24), 
  ev(amt = dose_point1*1000/(MW) , cmt = 'A1', time = 9 * 24), 
  ev(amt = dose_point1*1000/(MW) , cmt = 'A1', time = 13 * 24)
)

e_point4 <- c(
  ev(amt = dose_point4*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose_point4*1000/(MW) , cmt = 'A1', time = 5 * 24),
  ev(amt = dose_point4*1000/(MW) , cmt = 'A1', time = 9 * 24),
  ev(amt = dose_point4*1000/(MW) , cmt = 'A1', time = 13 * 24)
)

e_1point4 <- c(
  ev(amt = dose_1point4*1000/(MW) , cmt = 'A1', time = 0 * 24),
  ev(amt = dose_1point4*1000/(MW) , cmt = 'A1', time = 5 * 24),
  ev(amt = dose_1point4*1000/(MW) , cmt = 'A1', time = 9 * 24),
  ev(amt = dose_1point4*1000/(MW) , cmt = 'A1', time = 13 * 24)
)

e_5 <- c(
  ev(amt = dose_5*1000/(MW) , cmt = 'A1', time = 0 * 24),
  ev(amt = dose_5*1000/(MW) , cmt = 'A1', time = 5 * 24),
  ev(amt = dose_5*1000/(MW) , cmt = 'A1', time = 9 * 24),
  ev(amt = dose_5*1000/(MW) , cmt = 'A1', time = 13 * 24)
)

# simulation
N = 100 # pop size

sim_point1 <- mod %>% ev(e_point1) %>% mrgsim(delta = 10, end = 20 * 24, nid = N) %>% as_tibble() %>% select(CP, ROt) %>% filter(CP > 0.001, ROt > 0.1) %>% na.omit()

sim_point4 <- mod %>% ev(e_point4) %>% mrgsim(delta = 10, end = 20 * 24, nid = N) %>% as_tibble() %>% select(CP, ROt) %>% filter(CP > 0.001, ROt > 0.1) %>% na.omit()

sim_1point4 <- mod %>% ev(e_1point4) %>% mrgsim(delta = 10, end = 20 * 24, nid = N) %>% as_tibble() %>% select(CP, ROt) %>% filter(CP > 0.001, ROt > 0.1) %>% na.omit()

sim_5 <- mod %>% ev(e_5) %>% mrgsim(delta = 10, end = 20 * 24, nid = N) %>% as_tibble() %>% select(CP, ROt) %>% filter(CP > 0.001, ROt > 1) %>% na.omit()

simuldata <- rbind(sim_point1, sim_point4, sim_1point4, sim_5)

# observed data pulled from Lindauer et al.,2017
data <- read.csv('../data/Fig2b.csv', header = TRUE)

tumorreceptoroccupancy <- ggplot() + 
  geom_point(data = simuldata, aes(x = CP, y = ROt, colour = 'simulated'), alpha = 0.7, stroke = 0, size = 2) + 
  geom_point(data = data, aes(x = plasma_conc, y = ROt, colour = 'observed'), size = 3) + 
  labs(x = 'plasma mAb concentration (mg/L)', y = 'tumor receptor occupancy (%)', colour = '') +  
    scale_x_continuous(trans='log10', limits = c(0.1, 150)) +  theme(legend.position="bottom", legend.text=element_text(size=14))

# print(tumorreceptoroccupancy)

# save the plot
ggsave(
  '../img/Fig2b.png',
  plot = tumorreceptoroccupancy,
  scale = 1,
  dpi = 300,
  limitsize = FALSE,
)