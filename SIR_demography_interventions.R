###############################################################################################

#Adds to the basic SIR model (SIR with interventions)

###############################################################################################
library(deSolve)
library(ggplot2)

## SIER model with vaccination ##

#1. Define model function
#1-1. SIR model with vaccination (if vaccination provides perfect protection from infection)
SIR_V1<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + I + R
    
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + birth*N - death*S - mu*S
    dI <- beta*S*I/N - death*I - gamma*I
    dR <- gamma*I + mu*S - death*R 
    
    # return the rates of change as a list
    list(c(dS, dI, dR))
  })
}

#1-2. SIR model with vaccination ( if vaccination provides imperfect protection from infection)
SIR_V2<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V
    
    #compartments without vaccination
    dS_NotV <- -beta*S_NotV*I_NotV/N + birth*NotN - death*S_NotV - mu*S_NotV
    dI_NotV <- beta*S_NotV*I_NotV/N - death*I_NotV - gamma*I_NotV
    dR_NotV <- gamma*I_NotV - death*R_NotV 
    #compartments with vaccination
    dS_V <- -beta*S_V*I_V/N + birth*N - death*S_V + mu*S_NotV
    dI_V <- beta*S_V*I_V/N - death*I_V - gamma*I_V
    dR_V <- gamma*I_V - death*R_V 
    
    # return the rates of change as a list
    list(c(dS, dI, dR))
  })
}

#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                birth = 0.03, #birth rate (per capita)
                death = 0.03, #all-cause mortality rate
                omega = 0, # waning immunity
                mu = 0.6 # vaccination rate
)

state <- c(S = 99999, #population of 100,000, 1 person starts of infected
           I = 1, 
           R = 0
)


T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step 

#Run the base-case
output <- ode(y = state, times = times, func = SEIR, parms = parameters)



