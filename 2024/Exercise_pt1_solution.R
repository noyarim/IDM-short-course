#Set up and load packages
library(deSolve)
library(ggplot2)
library(dplyr)

#################
#Question 1     #
#################

#start by defining the basic SIR model function
BasicSIR <- function(t, state, params) {
  with(as.list(c(state, params)),{ 
    N = S + I + R #define N (total population size)
    
    #SIR model equations from slides - rates of change in and out of each compartment 
    dS <- -beta*S*I/N
    dI <- beta*S*I/N - gamma*I
    dR <- gamma*I
    
    #add a 4th state to track cumulative cases over time
    dC <- beta*S*I/N
    
    #return the rates of change as a list - IMPORTANT: must be same order as "state"
    list(c(dS, dI, dR, dC)) 
  })
}
# Define parameters and starting compartment sizes
# from the prompt, gamma = 1/2 weeks - but we need to estimate beta
gamma <- 1/14
R0_bounds <- c(3, 10)
beta_bounds <- gamma*R0_bounds
print(beta_bounds)

#################
#Question 2     #
#################

#run this model using best and worst case beta values
state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           I = 1, 
           R = 0,
           C = 0 #track cumulative number of infections
)
T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 
params <- c(gamma=gamma)

#Run with lower beta (R0=3)
params['beta'] <- beta_bounds[1]
out_best <- data.frame(ode(y = state, times = times, func = BasicSIR, parms = params))
#Run with higher beta (R0=10)
params['beta'] <- beta_bounds[2]
out_worst <- data.frame(ode(y = state, times = times, func = BasicSIR, parms = params))

#ANSWER TO 2A
print(paste("Best case:", round(out_best[T_end+1, 'C']), "cumulative cases"))
print(paste("Worst case:", round(out_worst[T_end+1, 'C']), "cumulative cases"))

#ANSWER TO 2B
print(paste("Best case epidemic peak occurs at", out_best[which.max(out_best[,"I"]),"time"], "days, with",
            round(max(out_best[,"I"])), "cases."))
print(paste("Worst case epidemic peak occurs at", out_worst[which.max(out_worst[,"I"]),"time"], "days, with",
            round(max(out_worst[,"I"])), "cases."))

#we can also plot trajectories over time
ggplot() + 
  geom_line(data=out_best, aes(x=time, y=I, color="Best Case")) +
  geom_line(data=out_worst, aes(x=time, y=I, color="Worst Case")) +
  labs(x="Day", y="Number of people infected", color="") +
  ggtitle("Cases over time")

