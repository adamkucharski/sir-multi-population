# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# R Code for SIR epidemic in 2 coupled populations
# Author: AJ Kucharski (2017-)
#
# Simulation functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Output simulated deterministic epidemic
Run_simulation<-function(dt, theta, theta_init,time.vals){
  
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
      
      N1 <- theta[["npop1"]]
      N2 <- theta[["npop2"]]
      
      ## extract parameters
      alpha <- theta[["alpha"]]
      beta <- theta[["beta"]]
      gamma <- theta[["gamma"]]
      repR <- theta[["rep"]]
      
      foi1 <- beta*(I+alpha*I2)/N1
      foi2 <- beta*(I2+alpha*I)/N2
      
      # Define transition equations
      dS <- -S*foi1
      dI <- S*foi1-gamma*I
      dR <- gamma*I
      dC <- repR*S*foi1
      
      dS2 <- -S2*foi2
      dI2 <- S2*foi2-gamma*I2
      dR2 <- gamma*I2
      dC2 <- repR*S2*foi2
      
      return(list(c(dS,dI,dR,dC,dS2,dI2,dR2,dC2)))
    }
    
    traj <- as.data.frame(ode(init.state, times, SIR_ode, theta, method = "ode45"))
    return(traj)
  }
  
  # Read in initial conditions
  init1 <- c(s_init=theta_init[["s_init"]],
          i1_init=theta_init[["i1_init"]],
          r_init=theta_init[["r_init"]],
          c_init=0,
          s2_init=theta_init[["s2_init"]],
          i2_init=theta_init[["i2_init"]],
          r2_init=theta_init[["r2_init"]],
          c2_init=0)
  
  # Output simulation data
  output <- simulate_deterministic(theta,init1,seq(0,max(time.vals),0.1))
  S_traj <- output[match(time.vals,output$time),"s_init"]
  I_traj <- output[match(time.vals,output$time),"i1_init"]
  
  # Calculate incidence in two populations
  cases1 <- output[match(time.vals,output$time),"c_init"]
  casecount <- cases1-c(0,cases1[1:(length(time.vals)-1)])
  
  cases2 <- output[match(time.vals,output$time),"c2_init"]
  casecount2 <- cases2-c(0,cases2[1:(length(time.vals)-1)])
  
  # Return selected outputs
  return(list(C1_trace=casecount,C2_trace=casecount2,C2_out=cases2,S_trace=S_traj,I_trace=I_traj))
  
}

