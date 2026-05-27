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
omega_list <- seq(0,0.3,by=0.1) # list of omegas
osa_dt <- data.frame() # empty data to save outputs

for (this_omega in omega_list){
  # Update omega
  parameters['omega']= this_omega
  # Solve ODE
  this_output <- as.data.frame(ode(y = state, times=times, func=OpenSIRS, parms = parameters))
  # Save omega value in the outcome table
  this_output$omega = as.character(this_omega)
  # Stack the ode result
  osa_dt <- rbind(osa_dt, this_output)
  
}

# Plot the results with varying rate of waning immunity
ggplot(osa_dt)+
    geom_line(aes(x=time, y=I, color=omega, group=omega))+
    ylab("Infected")+
    xlab("Time")+
    theme_bw()

#Run two-way sensitivity analysis on omega & beta
omega_list <- seq(0,0.3,by=0.1) # list of omegas
beta_list <- seq(0.4,0.6,by=0.1) # list of beta values
twsa_dt <- data.frame() # data to save the ode outcomes

for (this_beta in beta_list){
  # Update beta
  parameters['beta'] = this_beta
  for (this_omega in omega_list){
    # update omega
    parameters['omega'] = this_omega
    # Run ode solver
    this_output <- as.data.frame(ode(y = state, times=times, func=OpenSIRS, parms = parameters))
    # Save omega
    this_output$omega = as.character(this_omega)
    # Save beta
    this_output$beta = as.character(this_beta)
    # Stack the ode result
    twsa_dt <- rbind(twsa_dt, this_output)
  }
}
# Plot the results with varying omega and beta
# Change the label for omega and beta
omega_lb <- sapply(omega_list,function(x) paste0("omega=",x))
twsa_dt$omega <- factor(twsa_dt$omega, labels = omega_lb)
beta_lb <- sapply(beta_list,function(x) paste0("beta=",x))
twsa_dt$beta <- factor(twsa_dt$beta, labels = beta_lb)
# Check if label is created correctly
head(twsa_dt)
# Plot the two-way sensitivity analysis result
ggplot(twsa_dt)+
  geom_line(aes(x=time, y=I))+
  facet_grid(beta~omega)+
  ylab("Infected")+
  xlab("Time")+
  theme_bw()
