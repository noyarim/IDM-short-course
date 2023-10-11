###############################################################################################

#Adds to the basic SIR model (SEIR)

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)

## SEIR model (a model with latent state) ##

#1. Define model function
# SEIR model without births and deaths
SEIR_nodemo<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
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

# SEIR model with births and deaths
OpenSEIR<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
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
parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
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
out_seir_demo <- data.frame(ode(y = state, times = times, func = OpenSEIR, parms = parameters))

#Plot the result of base-case
ggplot(data=as.data.frame(out_seir_demo))+
  geom_point(aes(time,I))
out_seir_demo_t <- melt(out_seir_demo, id.vars="time")
out_seir_demo_t$demo <- "Births and deaths"
out_seir_comb <- rbind(out_seir_nodemo_t, out_seir_demo_t)

ggplot(data=out_seir_comb)+
  geom_line(aes(x=time,y=value, color=variable))+
  facet_wrap(.~demo)+
  theme_bw()

#Run the sensitivity analysis on the latent period (SEIR_demo)
tlat_list <- seq(0,100,by=0.1)
output_list <- data.frame()

for (this_tlat in tlat_list){
  
  parameters['t_lat'] = this_tlat
  this_output <- data.frame(ode(y = state, times=times, func=OpenSEIR, parms = parameters))
  this_output$t_lat = as.character(this_tlat)
  output_list <- rbind(output_list, this_output)
  
}

# Plot the results with varying latent period
ggplot(output_list)+
  geom_line(aes(x=time, y=I+E, color=t_lat, group=t_lat))+
  ylab("Infected")+
  xlab("Time")+
  theme_bw()+
  theme(legend.position='none')

# Discrete sensitivity analysis on latent period to show change in dynamics
tlat_list <- c(0.03,3,10)
output_list <- data.frame()
for (this_tlat in tlat_list){
  
  parameters['t_lat'] = this_tlat
  this_output <- data.frame(ode(y = state, times=times, func=OpenSEIR, parms = parameters))
  this_output$t_lat = as.character(this_tlat)
  output_list <- rbind(output_list, this_output)
  
}
ggplot(output_list)+
  geom_line(aes(x=time, y=I+E, color=t_lat, group=t_lat))+
  ylab("Infected")+
  xlab("Time")+
  theme_bw()

# Continuous sensitivity analysis showing the change in the epidemic peak
find_peak <- function(output){
  infected <- output$E + output$I
  peak <- which.max(infected)
  return(peak)
}

tlat_list <- seq(0.03,10,by=0.1)
output_list <- data.frame()
for (this_tlat in tlat_list){
  
  parameters['t_lat'] = this_tlat
  this_output <- data.frame(ode(y = state, times=times, func=OpenSEIR, parms = parameters))
  this_peak <- find_peak(this_output)
  this_output <- data.frame(peak=this_peak, tlat=this_tlat)
  output_list <- rbind(output_list, this_output)
  
}
ggplot(output_list)+
  geom_point(aes(x=tlat, y=peak))+
  ylab("epi peak")+
  xlab("tlat")+
  theme_bw()
