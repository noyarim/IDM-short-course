#Set up and load packages
library(deSolve)
library(ggplot2)
library(dplyr)

#################
#Question 1     #
#################

## Answers to the question 1 and 2 ##
# Update the SIRS to SEIRS model that incorporates latency
SEIRS<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S + E + I + R
    
    #SEIRS equations from lecture
    dS <- -beta*S*I/N + omega*R
    dE <- beta*S*I/N - sigma*E
    dI <- sigma*E - gamma*I
    dR <- gamma*I - omega*R
    
    #add a 4th state to track cumulative cases over time
    dC <- beta*S*I/N
    
    # return the rates of change as a list
    list(c(dS, dE, dI, dR, dC))
  })
}







#################
#Question 2     #
#################
#define model parameters (gamma, beta, omega same as in part 2)
gamma <- 1/14
R0_bounds <- c(3, 10)
beta_bounds <- gamma*R0_bounds
omega <- 1/(3*30.5) #1 over time to immunological waning - in days, since gamma is in days
# Add sigma (rate of transition from E to I = 1 over duration of latent period)
sigma <- 1/5 

params <- c(beta = beta_bounds[1], # place holder
            gamma = gamma,
            omega = omega,
            sigma = sigma)

#Also need to update "state" vector - to add E
state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           R = 0,
           C = 0 #track cumulative number of infections
)


#################
#Question 3a     #
#################
T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 

#Run with lower beta (R0=3) ,or "Best case scenario"
params['beta'] <- beta_bounds[1]
out_best3 <- data.frame(ode(y = state, times = times, func = SEIRS, parms = params))
#Run with higher beta (R0=10), or "Worst case scenario"
params['beta'] <- beta_bounds[2]
out_worst3 <- data.frame(ode(y = state, times = times, func = SEIRS, parms = params))
#Cumulative cases at Day 300
print(paste("Best case with latency:", round(out_best3[T_end+1, 'C']), "cumulative cases"))
print(paste("Worst case with latency:", round(out_worst3[T_end+1, 'C']), "cumulative cases"))




#################
#Question 3b     #
#################

# The likely size and timing of the epidemic peak
# When consider only I
# Best case scenario
print(paste("Including latent infections, best case epidemic peak occurs at",
            out_best3[which.max(out_best3[,"I"]),"time"], "days, with",
            round(max(out_best3[,"I"])), "cases."))
# Worst case scenario
print(paste("Including latent infections, worst case epidemic peak occurs at",
            out_worst3[which.max(out_worst3[,"I"]),"time"], "days, with",
            round(max(out_worst3[,"I"])), "cases."))

# Additional note: if we consider both E and I
out_best3 <- cbind(out_best3, Inf_all=out_best3$E+out_best3$I)
out_worst3 <- cbind(out_worst3, Inf_all=out_worst3$E+out_worst3$I)
# Best case scenario
print(paste("Including latent infections, best case epidemic peak occurs at",
            out_best3[which.max(out_best3[,"Inf_all"]),"time"], "days, with",
            round(max(out_best3[,"Inf_all"])), "infections."))
# Worst case scenario
print(paste("Including latent infections, worst case epidemic peak occurs at",
            out_worst3[which.max(out_worst3[,"Inf_all"]),"time"], "days, with",
            round(max(out_worst3[,"Inf_all"])), "infections."))



#################
#Question 3c     #
#################
# ANSWER: Including latency pushes out the timing of the epidemic peak, and also flattens it (reducing cumulative cases - because they are delayed)


# NEED TO RUN PART 1 and 2 solutions to compare them with the result in Part 3
# We can also plot trajectories over time, to compare across models (Part 2,3 in best and worst scenario)
ggplot() + 
  geom_line(data=out_best, aes(x=time, y=I, color="Best Case, no waning immunity or latency")) +
  geom_line(data=out_worst, aes(x=time, y=I, color="Worst Case, no waning immunity or latency")) +
  geom_line(data=out_best2, aes(x=time, y=I, color="Best Case, with waning immunity only")) +
  geom_line(data=out_worst2, aes(x=time, y=I, color="Worst Case, with waning immunity only")) +
  geom_line(data=out_best3, aes(x=time, y=I, color="Best Case, with latency & waning immunity")) +
  geom_line(data=out_worst3, aes(x=time, y=I, color="Worst Case, with latency & waning immmunity")) +
  labs(x="Day", y="Number of people infected", color="") +
  ggtitle("Cases over time")


