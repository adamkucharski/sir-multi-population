# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# R Code for SIR epidemic in 2 coupled populations
# Author: AJ Kucharski (2017-)
#
# Simulation functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define simulation model ODEs
simulate_deterministic <- function(theta, init.state, times) {
  SIR_ode <- function(time, state, theta) {
    
    ## define variables
    S <- state[["s_init"]]
    I <- state[["i1_init"]]
    R <- state[["r_init"]]
    C <-state[["c_init"]]
    
    S2 <- state[["s2_init"]]
    I2 <- state[["i2_init"]]
    R2 <- state[["r2_init"]]
    C2 <-state[["c2_init"]]
    
    S3 <- state[["s3_init"]]
    I3 <- state[["i3_init"]]
    R3 <- state[["r3_init"]]
    C3 <- state[["c3_init"]]
    
    N1 <- theta[["npop1"]]
    N2 <- theta[["npop2"]]
    N3 <- theta[["npop3"]]
    
    ## extract parameters
    alpha1 <- theta[["alpha1"]]
    alpha2 <- theta[["alpha2"]]
    alpha3 <- theta[["alpha3"]]
    beta <- theta[["beta"]]
    gamma <- theta[["gamma"]]
    repR <- theta[["rep"]]
    
    foi1 <- beta*(I/N1+alpha1*I2/N2+alpha3*I3/N1)
    foi2 <- beta*(I2/N2+alpha1*I/N1+alpha2*I3/N2)
    foi3 <- beta*(I3/N3+alpha2*I2/N2+alpha3*I/N1)
    
    # Define transition equations
    # Population 1
    dS <- -S*foi1
    dI <- S*foi1-gamma*I
    dR <- gamma*I
    dC <- S*foi1
    
    # Population 2
    dS2 <- -S2*foi2
    dI2 <- S2*foi2-gamma*I2
    dR2 <- gamma*I2
    dC2 <- S2*foi2
    
    # Population 3
    dS3 <- -S3*foi3
    dI3 <- S3*foi3-gamma*I3
    dR3 <- gamma*I3
    dC3 <- S3*foi3
    
    return(list(c(dS,dI,dR,dC,dS2,dI2,dR2,dC2,dS3,dI3,dR3,dC3)))
  }
  
  traj <- as.data.frame(ode(init.state, times, SIR_ode, theta, method = "ode45"))
  return(traj)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Output simulated deterministic epidemic
run_simulation <- function(dt, theta, theta_init,time.vals){

  # Read in initial conditions
  init1 <- c(s_init=theta_init[["s_init"]],
          i1_init=theta_init[["i1_init"]],
          r_init=theta_init[["r_init"]],
          c_init=0,
          s2_init=theta_init[["s2_init"]],
          i2_init=theta_init[["i2_init"]],
          r2_init=theta_init[["r2_init"]],
          c2_init=0,
          s3_init=theta_init[["s3_init"]],
          i3_init=theta_init[["i3_init"]],
          r3_init=theta_init[["r3_init"]],
          c3_init=0)
  
  # Output simulation data
  output <- simulate_deterministic(theta,init1,seq(0,max(time.vals),0.1))
  S_traj <- output[match(time.vals,output$time),"s_init"]
  I_traj <- output[match(time.vals,output$time),"i1_init"]
  
  # Calculate incidence in two populations
  cases1 <- output[match(time.vals,output$time),"c_init"]
  casecount <- cases1-c(0,cases1[1:(length(time.vals)-1)])
  
  cases2 <- output[match(time.vals,output$time),"c2_init"]
  casecount2 <- cases2-c(0,cases2[1:(length(time.vals)-1)])
  
  cases3 <- output[match(time.vals,output$time),"c3_init"]
  casecount3 <- cases3-c(0,cases3[1:(length(time.vals)-1)])
  
  # Return selected outputs
  return(list(C1_trace=casecount,C2_trace=casecount2,C3_trace=casecount3,C2_out=cases2,S_trace=S_traj,I_trace=I_traj))
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Poisson observation model

cases_reported <- function(case_series,theta) {
  
  sapply(case_series,function(x){rpois(1,lambda = theta[["rep"]]*x)})
  
}





