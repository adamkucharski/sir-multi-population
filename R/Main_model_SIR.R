# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# R Code for SIR epidemic in 2 coupled populations
# Author: AJ Kucharski (2017-)
#
# Main code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Set up source functions
library(deSolve)
library(magrittr)
library(tidyverse)

setwd("~/Documents/GitHub/sir-two-population/")

source("R/Functions_SIR.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define model parameters

popsize1 <- 10000 # population 1 size
popsize2 <- 8000 # population 2 size
dt <- 1  # time increment for simulation output (days)
time.vals <- seq(0,300,7) # time intervals

r0 <- 1.5 # basic reproduction number
  
# Set parameter vector
theta <- c(beta = NA,
          gamma = 1/7, # mean duration of infectiousness
          alpha = 0.005, # relative contact rate between the two populations
          npop1 = popsize1, # size of population 1
          npop2 = popsize2, # size of population 2
          rep = 0.25 # proportion of cases reported
          )

theta[["beta"]] <- r0*theta[["gamma"]] # define beta as function of R0 and gamma
  
# Set initial conditions - assume everyone susceptible
initial_inf <- 1 # initially infectious in population 1
theta_init <- c(s_init=NA,i1_init=initial_inf,r_init=0,s2_init=NA,i2_init=0,r2_init=0)
theta_init[["s_init"]] <- theta[["npop1"]]-theta_init[["i1_init"]]
theta_init[["s2_init"]] <- theta[["npop2"]]-theta_init[["i2_init"]]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run simulation
output_sim <- run_simulation(dt,theta,theta_init,time.vals)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot new reported cases over time
par(mfrow=c(2,1),mar = c(3.2,3.5,1,1),mgp=c(2.2,0.6,0),las=1)

# Total cases
plot(time.vals,output_sim$C1_trace,type="l",col='white',xlab="days",ylab="total cases",ylim=c(0,1000))
lines(time.vals, output_sim$C1_trace ,col='blue') # population 1
lines(time.vals, output_sim$C2_trace,col='red') # population 2
lines(time.vals, output_sim$C1_trace + output_sim$C2_trace,col='black')  # total cases

# Reported cases
plot(time.vals, output_sim$C1_trace,type="l",col='white',xlab="days",ylab="reported cases",ylim=c(0,300))
lines(time.vals, cases_reported(output_sim$C1_trace,theta) ,col='blue') # population 1
lines(time.vals, cases_reported(output_sim$C2_trace,theta),col='red') # population 2
lines(time.vals, cases_reported( output_sim$C1_trace + output_sim$C2_trace,theta) ,col='black')  # total cases


# Save plots
dev.copy(pdf,paste("plots/SIR_model.pdf",sep=""),width=8,height=6)
dev.off()