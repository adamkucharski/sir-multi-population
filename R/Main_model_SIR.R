# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# R Code for SIR epidemic in 2 coupled populations
# Author: AJ Kucharski (2017-)
#
# Main code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Set up source functions
library(deSolve)

setwd("~/Dropbox/LSHTM/Code_files/Code_examples/basic_SIR_meta_model/")
source("Functions_SIR.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define model parameters

popsize1 <- 10000 # population 1 size
popsize2 <- 10000 # population 2 size
dt <- 1  # time increment for simulation output (days)
time.vals <- seq(0,300,1) # time range
r0 <- 2 # basic reproduciton number
  
# Set parameter vector
theta <- c(beta=NA,
          gamma=1/14, # mean duration of infectiousness
          alpha=0.005, # relative contact rate between the two populations
          npop1=popsize1,
          npop2=popsize2,
          rep=0.5 # proportion of cases reported
          )

theta[["beta"]] <- r0*theta[["gamma"]] # define beta as function of R0 and gamma
  
# Set initial conditions - assume everyone susceptible
initial_inf <- 1 # initially infectious in population 1
theta_init <- c(s_init=NA,i1_init=initial_inf,r_init=0,s2_init=NA,i2_init=0,r2_init=0)
theta_init[["s_init"]] <- theta[["npop1"]]-theta_init[["i1_init"]]
theta_init[["s2_init"]] <- theta[["npop2"]]-theta_init[["i2_init"]]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run simulation
output_sim <- Run_simulation(dt,theta,theta_init,time.vals)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot new reported cases over time
par(mfrow=c(1,1),mar = c(3,3,1,1),mgp=c(1.8,0.6,0),las=1)

plot(time.vals,output_sim$C1_trace,type="l",col='white',xlab="days",ylab="cases",ylim=c(0,150))
lines(time.vals,output_sim$C1_trace,col='blue') # population 1
lines(time.vals,output_sim$C2_trace,col='red') # population 2
lines(time.vals,output_sim$C1_trace+ output_sim$C2_trace,col='black')  # total cases

# Save plots
dev.copy(pdf,paste("plots/SIR_model.pdf",sep=""),width=8,height=6)
dev.off()