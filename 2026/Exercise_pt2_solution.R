#Set up and load packages
library(deSolve)
library(ggplot2)
library(dplyr)

#################
#Question 1     #
#################

#Now we need an SIRS model - to include waning immunity
#define new model function
SIRS <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + I + R
    
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + omega*R
    dI <- beta*S*I/N - gamma*I
    dR <- gamma*I - omega*R
    
    #add a 4th state to track cumulative cases over time
    dC <- beta*S*I/N
    
    # return the rates of change as a list
    list(c(dS, dI, dR, dC))
  })
}

#################
#Question 2     #
#################

#define model parameters (gamma and beta same as in part 1)
gamma <- 1/14
R0_bounds <- c(3, 10)
beta_bounds <- gamma*R0_bounds
#add omega to parameters
omega <- 1/(3*30.5) #1 over time to immunological waning - in days, since gamma is in days
params <- c(gamma=gamma,
            omega=omega)

state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           I = 1, 
           R = 0,
           C = 0 #track cumulative number of infections
)
T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 

#################
#Question 3a    #
#################

#Run with lower beta (R0=3)
params['beta'] <- beta_bounds[1]
out_best2 <- data.frame(ode(y = state, times = times, func = SIRS, parms = params))
#Rerun with higher beta (R0=10)
params['beta'] <- beta_bounds[2]
out_worst2 <- data.frame(ode(y = state, times = times, func = SIRS, parms = params))

print(paste("Best case with immunological waning:", round(out_best2[T_end+1, 'C']), "cumulative cases"))
print(paste("Worst case with immunological waning:", round(out_worst2[T_end+1, 'C']), "cumulative cases"))
#we find that cumulative cases can exceed population size because each person can get infected more than once.

#################
#Question 3b    #
#################

print(paste("With immunological waning, best case epidemic peak occurs at",
            out_best2[which.max(out_best2[,"I"]),"time"], "days, with",
            round(max(out_best2[,"I"])), "cases."))

print(paste("With immunological waning, worst case epidemic peak occurs at",
            out_worst2[which.max(out_worst2[,"I"]),"time"], "days, with",
            round(max(out_worst2[,"I"])), "cases."))

ggplot() + 
  geom_line(data=out_best2, aes(x=time, y=I, color="Best Case, with waning immunity")) +
  geom_line(data=out_worst2, aes(x=time, y=I, color="Worst Case, with waning immunity")) +
  labs(x="Day", y="Number of people infected", color="") +
  ggtitle("Cases over time")
#Note - In the above figure, it looks like a steady state has been reached in the worst case version, when waning immunity is included.

#################
#Question 4     #
#################

#R0 doesn't change, but the effective reproductive number does
#You weren't asked to calculate the effective reproductive number
#But if you wanted to:
Rt_best2 <- R0_bounds[1]*out_best2$S/sum(state)
Rt_worst2 <- R0_bounds[2]*out_worst2$S/sum(state)

Rt_all <- rbind(data.frame("Rt"=Rt_best2, "version"="Best case, with waning immunity"),
                data.frame("Rt"=Rt_worst2, "version"="Worst case, with waning immunity"))
Rt_all <- cbind("time"=out_best2$time, Rt_all)

ggplot(Rt_all, aes(x=time, y=Rt, color=version)) +
  geom_line() +
  labs(x="Day", y="Effective reproductive number", color="") +
  ggtitle("Effective reproductive number over time")
#Notice how the effective reproductive number stabilizes at 1 when we incorporate waning immunity

#################
#Question 5     #
#################

#Intuition is that incorporating demography won't matter much over the short (300 day) horizon
#Here, we show this with an updated (open) model 

OpenSIRS<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + I + R
    
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + birth*N - death*S + omega*R
    dI <- beta*S*I/N - death*I - gamma*I
    dR <- gamma*I - death*R - omega*R
    
    #add a 4th state to track cumulative cases over time
    dC <- beta*S*I/N
    
    # return the rates of change as a list
    list(c(dS, dI, dR, dC))
  })
}
params <- c(gamma=gamma,
            omega=omega, #1 over time to immunological waning
            birth=(120/10000)*(1/365), #births per population per day
            death=(120/10000)*(1/365)) #births per population per day

#Run with lower beta (R0=3)
params['beta'] <- beta_bounds[1]
out_best2_demo <- data.frame(ode(y = state, times = times, func = OpenSIRS, parms = params))
#Run with higher beta (R0=10)
params['beta'] <- beta_bounds[2]
out_worst2_demo <- data.frame(ode(y = state, times = times, func = OpenSIRS, parms = params))
ggplot() + 
  geom_line(data=out_best2, aes(x=time, y=I, color="Best Case, without demography")) +
  geom_line(data=out_worst2, aes(x=time, y=I, color="Worst Case, without demography")) +
  geom_line(data=out_best2_demo, aes(x=time, y=I, color="Best Case, with demography")) +
  geom_line(data=out_worst2_demo, aes(x=time, y=I, color="Worst Case, with demography")) +
  labs(x="Day", y="Number of people infected", color="")
#with and without demography lines completely overlap


