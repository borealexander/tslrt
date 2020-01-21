## ----setup, include=FALSE------------------------------------------------

knitr::opts_chunk$set(echo = TRUE)


## ----packages, echo = FALSE, message = FALSE, warning = FALSE------------

library(data.table)
library(survival)
library(ggplot2)
library(dplyr)
library(tslrt)
library(knitr)
library(kableExtra)



## ----parameter_setup_delay-----------------------------------------------

# number of patients in each arm
n_c <- 100
n_e <- 100

# median in control before switch
median_c_1 <- 10
# median in control after switch
median_c_2 <- 15
# median in experimental
median_e <- 15

# delay time where switching occurs
delay <- 3

# number of events to stop (event-driven study)
end_event <- 150

# recruitment period
rec_period <- 12
# recruitment parameter (assuming recruitment follows a powermodel)
rec_power <- 1

# different proportion of patients that switches
p <- c(0, 0.25, 0.5, 0.75, 1)

# significance level (one-sided)
alpha <- 0.025

# number of simulations
M <- 100


## ----MC_delay------------------------------------------------------------

res_delay <- MC_delay_crossover(n_c = n_c,
                                n_e = n_e,
                                median_c_1 = median_c_1,
                                median_c_2 = median_c_2,
                                median_e = median_e,
                                delay = delay,
                                end_event = end_event,
                                rec_period = rec_period,
                                rec_power = rec_power,
                                p = p, 
                                alpha = alpha,
                                M = M)


## ----MC_delay_table------------------------------------------------------

res_delay$result


## ----MC_delay_graph, fig.width = 7, fig.height = 5-----------------------

res_delay$plot


## ----parameter_setup_prog------------------------------------------------

# number of patients in each arm
n_c <- 100
n_e <- 100

# median in control before switch
median_c <- 10
# median in experimental
median_e <- 15
# median for progressions
median_prog <- 3

# number of events to stop (event-driven study)
end_event <- 150

# recruitment period
rec_period <- 12
# recruitment parameter (assuming recruitment follows a powermodel)
rec_power <- 1

# different proportion of patients that switches
p <- c(0, 0.25, 0.5, 0.75, 1)

# significance level (one-sided)
alpha <- 0.025

# number of simulations
M <- 100


## ----MC_prog-------------------------------------------------------------

res_prog <- MC_exp_prog_crossover(n_c = n_c,
                                  n_e = n_e,
                                  median_c = median_c,
                                  median_e = median_e,
                                  median_prog = median_prog,
                                  end_event = end_event,
                                  rec_period = rec_period,
                                  rec_power = rec_power,
                                  p = p, 
                                  alpha = alpha,
                                  M = M)


## ----MC_prog_table-------------------------------------------------------

res_prog$result


## ----MC_prog_graph, fig.width = 7, fig.height = 5------------------------

res_prog$plot


