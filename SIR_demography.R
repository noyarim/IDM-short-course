###############################################################################################

#SIR model with births and deaths

###############################################################################################

library(deSolve) #differential equation solver
library(tidyverse) 
library(PRIMsrc)

#1. Define model function
OpenSIR<-function(t, state, parameters) {
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
                 omega = 0.0 # waning immunity
)

state <- c(S = 99999, #population of 100,000, 1 person starts of infected
           I = 1, 
           R = 0
)

T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step 

#3. Calculate R0 and equilibrium distributions
R0 <- parameters[["beta"]]/(parameters[["gamma"]] + parameters[["death"]])
print(R0)

S_star <- 1/R0
I_star <- (parameters[["death"]]/parameters[["beta"]])*(R0-1)
R_star <- 1-(S_star + I_star)
print(c("S_star"=S_star, "I_star"=I_star, "R_star"=R_star)) #prevalences at equilibrium

#4. Run model
output <- ode(y = state, times = times, func = OpenSIR, parms = parameters)

output_long <- pivot_longer(as.data.frame(output), cols=2:ncol(output), names_to="state", values_to="size")
plot_trace(output_long)

output_df <- data.frame(output)
output_df$N <- output_df$S + output_df$I + output_df$R
print(output_df$I[T_end+1]/output_df$N[T_end+1])

