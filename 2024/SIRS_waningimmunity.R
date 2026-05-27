###############################################################################################

#SIRS model (closed population) - adds waning immunity to the basic SIR model

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)

## SIRS model (a model with waning immunity) ##
# Backward flow from R to S indicates those who lost infection-induced immunity over time.

#1. Define model function
SIRS<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S + I + R

    #SIRS equations from lecture
    dS <- -beta*S*I/N  + omega*R
    dI <- beta*S*I/N - gamma*I
    dR <- gamma*I - omega*R

    # return the rates of change as a list
    list(c(dS, dI, dR))
  })
}

#2. Define parameters and starting compartment sizes
# Omega is 1/duration of immunity, assuming a constant rate of waning immunity.
params <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                omega = 0.3 # waning immunity
)

state <- c(S = 99999, #population of 100,000, 1 person starts of infected
           I = 1,
           R = 0
)


#3. Run SIRS model
T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step
output <- ode(y = state, times = times, func = SIRS, parms = params)

#Plot the number of individuals in S, I, R compartments over time
output_t <- melt(as.data.frame(output), id.vars='time')
ggplot(output_t)+
  geom_line(aes(time,value,color=variable))+
  theme_bw()

#4. One-way sensitivity analysis on the rate of waning immunity
omega_list <- seq(0,0.3,by=0.1) # a vector of omegas
osa_dt <- data.frame() # empty data to save outputs

for (this_omega in omega_list){
  # Update omega
  params['omega']= this_omega
  # Solve ODE
  this_output <- as.data.frame(ode(y = state, times=times, func=SIRS, parms = params))
  # Save omega value in the outcome table
  this_output$omega = as.character(this_omega)
  # Stack the ode result
  osa_dt <- rbind(osa_dt, this_output)

}

# Plot the results with varying rate of waning immunity
# How does changing waning immunity affect the disease dynamic?
ggplot(osa_dt)+
    geom_line(aes(x=time, y=I, color=omega, group=omega))+
    ylab("Infected")+
    xlab("Time")+
    theme_bw()

#5. Two-way sensitivity analysis on omega & beta
# Effect of waning immunity on epidemic growth can depend on how fast the disease spreads (beta).
omega_list <- seq(0,0.3,by=0.1) # a vector of omegas
beta_list <- seq(0.4,0.6,by=0.1) # a vector of beta values
twsa_dt <- data.frame() # data to save the ode outcomes

for (this_beta in beta_list){
  # Update beta
  params['beta'] = this_beta
  for (this_omega in omega_list){
    # update omega
    params['omega'] = this_omega
    # Run ode solver
    this_output <- as.data.frame(ode(y = state, times=times, func=SIRS, parms = params))
    # Save omega
    this_output$omega = as.character(this_omega)
    # Save beta
    this_output$beta = as.character(this_beta)
    # Stack the ode result
    twsa_dt <- rbind(twsa_dt, this_output)
  }
}

# Change the label for omega and beta
omega_lb <- sapply(omega_list,function(x) paste0("omega=",x))
omega_lb
twsa_dt$omega <- factor(twsa_dt$omega, labels = omega_lb)
beta_lb <- sapply(beta_list,function(x) paste0("beta=",x))
beta_lb
twsa_dt$beta <- factor(twsa_dt$beta, labels = beta_lb)
# Check if label is created correctly
head(twsa_dt)

# Plot the two-way sensitivity analysis result.
# Has the effect of waning immunity on epidemic changed with different betas?
ggplot(twsa_dt)+
  geom_line(aes(x=time, y=I))+
  facet_grid(beta~omega)+
  ylab("Infected")+
  xlab("Time")+
  theme_bw()
