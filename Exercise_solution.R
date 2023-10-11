###############################################################################################

#Exercise solution

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)

## SEIRS model (a model with latent state and waning immunity) ##

#1. Define model function
# SEIRS model with births and deaths
OpenSEIRS<-function(t, state, parameters) {
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
parameters <- c(beta = 0.9, #effective contact rate (aka transmission rate)
                gamma = 1/3, #recovery rate (1/duration infection)
                birth = 0.01, #birth rate (per capita)
                death = 0.01, #all-cause mortality rate
                omega = 1/(365*3), # waning immunity
                t_lat = 5 # latent period from E
                
)

state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           R = 0
)


T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 

#Run the base-case with beta = 0.5
parameters['beta'] <- 0.5
out_lowbeta <- data.frame(ode(y = state, times = times, func = OpenSEIRS, parms = parameters))
out_lowbeta$totI <- out_lowbeta$E + out_lowbeta$I
#Run the base-case with beta = 0.5
parameters['beta'] <- 0.9
out_highbeta <- data.frame(ode(y = state, times = times, func = OpenSEIRS, parms = parameters))
out_highbeta$totI <- out_highbeta$E + out_highbeta$I

# Plot S,E,I,R over time
out_lowbeta_t <- melt(out_lowbeta, id.vars='time')
ggplot(out_lowbeta_t)+
  geom_line(aes(time,value,color=variable))+
  ggtitle("Beta = 0.5")+
  theme_bw()

out_highbeta_t <- melt(out_highbeta, id.vars='time')
ggplot(out_highbeta_t)+
  geom_line(aes(time,value,color=variable))+
  ggtitle("Beta = 0.9")+
  theme_bw()

# Plot E + I over time & peak of epidemics
totI <- data.frame(time = out_lowbeta$time, lowbeta = out_lowbeta$totI, highbeta = out_highbeta$totI)
ggplot(totI)+
  geom_line(aes(time,lowbeta, color='low beta'))+
  geom_line(aes(time,highbeta, color='high beta'))+  
  ylab("E + I")+
  annotate(geom="vline", 
           x = c(which.max(totI$lowbeta),which.max(totI$highbeta)),
           xintercept = c(which.max(totI$lowbeta),which.max(totI$highbeta)),
           linetype = c("dashed","dashed"),
           color = c("#0072B2","#D55E00"))+
  annotate(geom="text",
           label = c(as.character(which.max(totI$lowbeta)),as.character(which.max(totI$highbeta))),
           x = c(which.max(totI$lowbeta),which.max(totI$highbeta)),
           y = c(max(totI$lowbeta), max(totI$highbeta)),
           color = c("#0072B2","#D55E00"),
           hjust = -0.5)+
  ggtitle("Total number of infections (E+I)")+
  theme_bw()+
  scale_color_manual(values = c("#D55E00","#0072B2"))
  
  


