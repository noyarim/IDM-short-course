###############################################################################################

#Adds to the basic SIR model (SIRS)

###############################################################################################
library(deSolve)
library(ggplot2)

## SIRS model (a model with waning immunity) ##

#1. Define model function
OpenSIRS<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + I + R
    
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + birth*N - death*S + omega*R
    dI <- beta*S*I/N - death*I - gamma*I
    dR <- gamma*I - death*R - omega*R
    
    # return the rates of change as a list
    list(c(dS, dI, dR))
  })
}

#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                birth = 0.03, #birth rate (per capita)
                death = 0.03, #all-cause mortality rate
                omega = 0.3 # waning immunity
)

state <- c(S = 99999, #population of 100,000, 1 person starts of infected
           I = 1, 
           R = 0
)


T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step 

#Run the base-case
output <- ode(y = state, times = times, func = OpenSIRS, parms = parameters)

#Plot the result from base-case
plot(output)

#Run one-way sensitivity analysis on the rate of waning immunity
omega_list <- seq(0,0.3,by=0.1)
osa_list <- data.frame()

for (this_omega in omega_list){
  
  parameters$omega = this_omega
  this_output <- as.data.frame(ode(y = state, times=times, func=OpenSIRS, parms = parameters))
  this_output$omega = as.character(this_omega)
  osa_list <- rbind(osa_list, this_output)
  
}

# Plot the results with varying rate of waning immunity
ggplot(osa_list)+
    geom_line(aes(x=time, y=I, color=omega, group=omega))+
    ylab("Infected")+
    xlab("Time")+
    theme_bw()

#Run two-way sensitivity analysis on omega & beta
beta_list <- seq(0.4,0.6,by=0.1)
twsa_list <- data.frame()

for (this_beta in beta_list){
  parameters$beta = this_beta
  for (this_omega in omega_list){
    
    parameters$omega = this_omega
    this_output <- as.data.frame(ode(y = state, times=times, func=OpenSIRS, parms = parameters))
    this_output$omega = as.character(this_omega)
    this_output$beta = as.character(this_beta)
    twsa_list <- rbind(twsa_list, this_output)
  }
}
# Plot the results with varying gamma and beta
ggplot(twsa_list)+
  geom_line(aes(x=time, y=I))+
  facet_grid(beta~omega)+
  ylab("Infected")+
  xlab("Time")+
  theme_bw()
