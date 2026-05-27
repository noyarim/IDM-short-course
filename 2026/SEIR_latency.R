###############################################################################################

#SEIR model with latent state (adds latency to the basic SIR model)

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)

## SEIR model (a model with latent state) ##
# E compartment indicates those who were exposed to the infectious agent but have not developed infectiousness.
# Rate of progressing from E to I is internally calculated based on the latent period (sigma = 1/t_lat).
# Note that E does not infect others (new infections = beta*S*I/N, not beta*S*(E+I)/N).

#1. Define model function
# SEIR model with births and deaths
OpenSEIR<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S + E + I + R
    sigma = 1/t_lat # 1/latent period

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
# t_lat is the duration of latency before getting infectious.
# In this example, we assumed no waning immunity (omega=0) but this assumption can be easily relaxed.
params <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                birth = 0.03, #birth rate (per capita)
                death = 0.03, #all-cause mortality rate
                omega = 0, # waning immunity
                t_lat = 3 # latent period from E

)

state <- c(S = 99999, #population of 100,000, 1 person starts of infected
           E = 0,
           I = 1,
           R = 0
)


T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step

#Run the base-case
out_seir_demo <- data.frame(ode(y = state, times = times, func = OpenSEIR, parms = params))

#Plot the number of I compartment over time
ggplot(data=as.data.frame(out_seir_demo))+
  geom_point(aes(time,I))+
  ylab("Infected")+
  theme_bw()

#3. One-way sensitivity analysis on the latent period - see how infections change with varying duration of latency
tlat_list <- c(0.01,seq(1,5,by=1)) # a vector of latent period
output_dt <- data.frame() # an empty dataset to save outcomes

for (this_tlat in tlat_list){
  # Replace latent period with the next value in the last
  params['t_lat'] = this_tlat
  # Run ode solver
  this_output <- data.frame(ode(y = state, times=times, func=OpenSEIR, parms = params))
  # Record current value of t_lat
  this_output$t_lat = as.character(this_tlat)
  # Stack the result
  output_dt <- rbind(output_dt, this_output)

}

# Plot the number of infections (size of I compartment) with varying latent period.
# Here we plot only I to assess the size of population who are infectious;
# if you want 'true' number of infections, compute sum of I and E compartments.
# With longer disease latency, how does the epidemic curve change?
ggplot(output_dt)+
  geom_line(aes(x=time, y=I, color=t_lat, group=t_lat))+
  ylab("Infected")+
  xlab("Time")+
  theme_bw()

#4. Continuous sensitivity analysis showing the change in the epidemic peak
# A function to find the epidemic peak given a time frame
find_peak <- function(output){
  infected <- output$I # output$E + output$I if you want to include exposed compartment
  peak <- which.max(infected)
  return(peak)
}

tlat_list <- seq(0.03,10,by=0.1) # a vector of t_lat
output_dt <- data.frame() # An empty dataset to save ode outcomes
for (this_tlat in tlat_list){
  # Replace t_lat with the next value
  params['t_lat'] = this_tlat
  # Run ode solver
  this_output <- data.frame(ode(y = state, times=times, func=OpenSEIR, parms = params))
  # Find the epidemic peak
  this_peak <- find_peak(this_output)
  # Create a dataset with peak and t_lat in this cycle
  this_output <- data.frame(peak=this_peak, tlat=this_tlat)
  # Stack the dataset
  output_dt <- rbind(output_dt, this_output)

}
# Plot the relationship between epidemic peak and latent period
ggplot(output_dt)+
  geom_point(aes(x=tlat, y=peak))+
  ylab("Epi peak")+
  xlab("Latent period")+
  theme_bw()


## Additional code for SEIR model without demography
# SEIR model without births and deaths
SEIR_nodemo<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S + E + I + R
    sigma = 1/t_lat # 1/latent period

    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + omega*R
    dE <- beta*S*I/N - sigma*E
    dI <- sigma*E - gamma*I
    dR <- gamma*I - omega*R

    # return the rates of change as a list
    list(c(dS, dE, dI, dR))
  })
}
