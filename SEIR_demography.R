###############################################################################################

#Adds to the basic SIR model (SEIR)

###############################################################################################
library(deSolve)
library(ggplot2)

## SIRS model (a model with waning immunity) ##

#1. Define model function
SEIR<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + E + I + R
    
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + birth*N - death*S + omega*R
    dE <- beta*S*I/N - sigma*E - death*E
    dI <- sigma*E - death*I - gamma*I
    dR <- gamma*I - death*R - omega*R
    
    # return the rates of change as a list
    list(c(dS, dE, dI, dR))
  })
}

#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                birth = 0.03, #birth rate (per capita)
                death = 0.03, #all-cause mortality rate
                omega = 0.3, # waning immunity
                sigma = 0.1 # 1/incubation period
)

state <- c(S = 99999, #population of 100,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           R = 0
)


T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step 

#Run the base-case
output <- ode(y = state, times = times, func = SEIR, parms = parameters)

#Run the sensitivity analysis on the rate of waning immunity
sigma_list <- seq(0,0.5,by=0.1)
output_list <- data.frame()

for (this_sigma in sigma_list){
  
  parameters['sigma'] = this_sigma
  this_output <- as.data.frame(ode(y = state, times=times, func=SEIR, parms = parameters))
  this_output$sigma = as.character(this_sigma)
  output_list <- rbind(output_list, this_output)
  
}

# Plot the results with varying rate of waning immunity
ggplot(output_list)+
  geom_line(aes(x=time, y=I+E, color=sigma, group=sigma))+
  ylab("Infected")+
  xlab("Time")+
  theme_bw()

# With death and birth turn on and off
# Turn off birth and death and introduce E first to show how it affects the timing of epidemics
# Then turn on birth and death and introduce E to show how it interacts to produce different prevalence in equibrilium
