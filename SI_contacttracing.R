###############################################################################################

#Adds to the basic SI model (SI with contact tracing)

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)
library(tidyverse)

# 1. Define model function
SI_CC<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + I + Tx
    
    #SI w/o demography equations from lecture
    dS <- -beta*S*I/N + rt*Tx
    dI <- beta*S*I/N - di*I - c*b*Tx
    dTx <- di*I + c*b*Tx - rt*Tx
    # Cumulative incidence 
    dC <- beta*S*I/N
    # return the rates of change as a list
    list(c(dS, dI, dTx, dC))
  })
}
#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                rt = 0.1, # 1/duration of treatment and tracing
                di = 0.1, # rate of seeking treatment
                c = 0.05, # rate of contact tracing
                b = 0.05# probability that a traced individual is infected
)
# Initial state
state <- c(S = 99999, #population of 100,000, 1 person starts of infected
             I = 1, 
             Tx = 0,
             C = 0)
T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step 

#Run the base-case
output <- ode(y = state, times = times, func = SI_CC, parms = parameters)

C <- output[,'C']

output[,-c('C')]
# Plot result from base-case
output_t <- melt(as.data.frame(output),id.vars='time')
ggplot(output_t)+
  geom_point(aes(time,value,color=variable))+
  theme_bw()

# Sensitivity analysis on contact tracing rate
c_list <- c(0,0.05,0.1) # a vector of contact tracing rate
output_dt <- data.frame() # empty data to save outcomes

for (this_c in c_list){
  temp <- paramters
  # Replace latent period with the next value in the last
  temp['c'] <- this_c
  # Run ode solver
  this_output <- data.frame(ode(y = state, times = times, func = SI_CC, parms = parameters))
  # Record current value of t_lat
  this_output$c = as.character(c)
  # Stack the result 
  output_dt <- rbind(output_dt, this_output)
  
}

# Plot the results with varying latent period
ggplot(output_dt)+
  geom_line(aes(x=time, y=I, color=c, group=c))+
  ylab("Infected")+
  xlab("Time")+
  theme_bw()
