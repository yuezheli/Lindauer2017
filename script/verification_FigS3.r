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

# FigS3 dose regime: 
# antibody: Rat (80ADW)
# Administration route: intravenous
# sex/ strain: female, C56BL/6
# regumen: multidose, day 0 and day 4
# sampling time: 24 hours after day 4 after first & second doses

## ------------- observed data from Fig S3 ------------- ##
data_0 <- read.csv('../data/FigS3/dose_0.csv', header = TRUE)
data_point1 <- read.csv('../data/FigS3/dose_point1.csv', header = TRUE)
data_point4 <- read.csv('../data/FigS3/dose_point4.csv', header = TRUE)
data_1point4 <- read.csv('../data/FigS3/dose_1point4.csv', header = TRUE)
data_5 <- read.csv('../data/FigS3/dose_5.csv', header = TRUE)


## ------------- pop simulation ------------- ##
N = 100 # pop size

# dosing amount
dose_point1 = 0.1 * mouseweight
dose_point4 = 0.4 * mouseweight
dose_1point4 = 1.4 * mouseweight
dose_5 = 5 * mouseweight

# dosing regiment
e_point1 <- c(
  ev(amt = dose_point1*1000/(MW) , cmt = 'A1', time = 0 * 24), 
  ev(amt = dose_point1*1000/(MW) , cmt = 'A1', time = 4 * 24)
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


# simulation 
sim0 <- mod %>% 
  mrgsim(delta = 1, end = 20 * 24, nid = N) %>%
  as_tibble() %>% select(time, TV) %>% group_by(time) %>% 
  summarize(lo=quantile(TV, 0.1), hi=quantile(TV, 0.9), med = quantile(TV, 0.5))


sim_point1 <- mod %>% ev(e_point1) %>%
  mrgsim(delta = 1, end = 20 * 24, nid = N) %>%
  as_tibble() %>% select(time, TV) %>% group_by(time) %>% 
  summarize(lo=quantile(TV, 0.1), hi=quantile(TV, 0.9), med = quantile(TV, 0.5))

sim_point4 <- mod %>% ev(e_point4) %>%
  mrgsim(delta = 1, end = 20 * 24, nid = N) %>%
  as_tibble() %>% select(time, TV) %>% group_by(time) %>% 
  summarize(lo=quantile(TV, 0.1), hi=quantile(TV, 0.9), med = quantile(TV, 0.5))

sim_1point4 <- mod %>% ev(e_1point4) %>%
  mrgsim(delta = 1, end = 20 * 24, nid = N) %>%
  as_tibble() %>% select(time, TV) %>% group_by(time) %>% 
  summarize(lo=quantile(TV, 0.1), hi=quantile(TV, 0.9), med = quantile(TV, 0.5))

sim_5 <- mod %>% ev(e_5) %>%
  mrgsim(delta = 1, end = 20 * 24, nid = N) %>%
  as_tibble() %>% select(time, TV) %>% group_by(time) %>% 
  summarize(lo=quantile(TV, 0.1), hi=quantile(TV, 0.9), med = quantile(TV, 0.5))

## ------------- Plot simulation result and compare against observed data ------------- ##

TumorVolume_0 <- ggplot() + 
  #  geom_line(data = sim0, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim0, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.3) + 
  geom_line(data = data_0, aes(x = time, y = tumor_volume, color = 'observed'), size = 1) + 
  labs(y = 'tumor volume (uL)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 0mg/kg') + xlim(0.01, 18) + 
  scale_y_continuous(trans='log10', limits = c(1,3000)) + theme(legend.position = "bottom", legend.title = element_blank()) 

TumorVolume_point1 <- ggplot() + 
  #  geom_line(data = sim_point1, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim_point1, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.3) + 
  geom_line(data = data_point1, aes(x = time, y = tumor_volume, color = 'observed'), size = 1) + 
  labs(y = 'tumor volume (uL)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 0.1mg/kg') + xlim(0.01, 18) + 
  scale_y_continuous(trans='log10', limits = c(1,3000)) + theme(legend.position = "bottom", legend.title = element_blank()) 

TumorVolume_point4 <- ggplot() + 
  #  geom_line(data = sim_point4, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim_point4, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.3) +   
  geom_line(data = data_point4, aes(x = time, y = tumor_volume, color = 'observed'), size = 1) + 
  labs(y = 'tumor volume (uL)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 0.4mg/kg') + xlim(0.01, 18) + 
  scale_y_continuous(trans='log10', limits = c(1,3000)) + theme(legend.position = "bottom", legend.title = element_blank())

TumorVolume_1point4 <- ggplot() +
  #  geom_line(data = sim_1point4, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim_1point4, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.3) + 
  geom_line(data = data_1point4, aes(x = time, y = tumor_volume, color = 'observed'), size = 1) + 
  labs(y = 'tumor volume (uL)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 1.4mg/kg') + xlim(0.01, 18) + 
  scale_y_continuous(trans='log10', limits = c(1,3000))  + theme(legend.position = "bottom", legend.title = element_blank())

TumorVolume_5 <- ggplot() + 
  #  geom_line(data = sim_5, aes(x = time/24, y = med, color = 'median, predicted'), size = 0.5) + 
  geom_ribbon(data = sim_5, aes(x = time/24, ymin=lo,ymax=hi, color = 'predicted'), fill="gray", alpha = 0.3) + 
  geom_line(data = data_5, aes(x = time, y = tumor_volume, color = 'observed'), size = 1) + 
  labs(y = 'tumor volume (uL)', x = 'time (days)', 
       color = 'data source') + ggtitle('dose = 5mg/kg') + xlim(0.01, 18) + 
  scale_y_continuous(trans='log10', limits = c(1,3000)) + theme(legend.position = "bottom", legend.title = element_blank()) 

blankgraph <- textGrob(" ")

figs3 = grid.arrange(TumorVolume_1point4, TumorVolume_5, blankgraph, 
             TumorVolume_0, TumorVolume_point1, TumorVolume_point4, ncol=3)

# save the plot
ggsave(
  '../img/FigS3.png',
  plot = figs3,
  scale = 1.2,
  dpi = 300,
  limitsize = FALSE,
)