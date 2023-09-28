###############################################################################################

#SIR model with births and deaths

###############################################################################################

library(deSolve) #differential equation solver
library(tidyverse)

plot_trace <-function(out) {
  fig <- ggplot(out, aes(x=time, y=size/N, color=state)) + 
    geom_line(linewidth=1.25) +
    labs(x='Time', y='Compartment size', color='') +
    theme_bw() + theme(panel.grid=element_blank())
  return(fig)
}

#1. Define model function
OpenSIR<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + I + R
    
    #SIR w/ demography equations from the slides
    dS <- -beta*S*I/N + birth*N - death*S
    dI <- beta*S*I/N - death*I - gamma*I
    dR <- gamma*I - death*R
    
    # return the rates of change as a list
    list(c(dS, dI, dR))
  })
}

#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.5, #effective contact rate
                 gamma = 0.3, #recovery rate (1/duration infection)
                 birth = 0.02, #birth rate (per capita)
                 death = 0.02 #all-cause mortality rate
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

#4. Run model
output <- ode(y = state, times = times, func = OpenSIR, parms = parameters)
output <- as.data.frame(output) %>% mutate(N=S+I+R)
output_long <- pivot_longer(output, cols=c("S","I","R"), names_to="state", values_to="size")
plot_trace(output_long)

print(output$I[T_end+1]/output$N[T_end+1])
print(c("S_star"=S_star, "I_star"=I_star, "R_star"=R_star)) #prevalences at equilibrium
#see how this matches the equilibrium distribution

#5. Calculate Rt
output <- output %>% mutate(Rt=R0*S/N)

#changes in Rt corerlate w/ the oscillatory pattern in the trace graph
#changes in Rt are fueled by rises and falls in the susceptible population (from births and new infections, respectively)
plot_trace(output_long) +
  geom_line(data=output, aes(x=time, y=Rt), color="black") +
  geom_hline(yintercept=1, linetype="dashed", color="red")

#6. How do dynamics change as birth/death rate changes?
output_all <- list()
for(rates in c(0, 0.01, 0.03, 0.1)) {
    parameters[["birth"]] <- parameters[["death"]] <- rates
    output <- ode(y = state, times = times, func = OpenSIR, parms = parameters)
    output <- as.data.frame(output) %>%
      mutate(N=S+I+R,
             demo_rate=rates,
             R0=parameters[["beta"]]/(parameters[["gamma"]] + rates),
             Rt=R0*S/N,
             demo_lab=paste0("Birth/death rate: ", demo_rate),
             R0_lab=paste0("R0: ", round(R0, 2)),
             I_star_lab=paste0("% I at t=500: ", round(100*I[time==500]/N[time==500], 1)))
    output_long <- pivot_longer(output, cols=c("S","I","R"), names_to="state", values_to="size")
    output_all <- c(output_all, list(output_long))
}

output_all <- bind_rows(output_all)

plot_trace(output_all) + facet_wrap(~demo_lab+R0_lab+I_star_lab)
# lower birth/death rate:
# takes more time to reach a steady state (when Rt=1)
# higher R0 (less likely for someone to die before infecting someone else)
# lower % infected at steady state (slower replenishment of susceptibles)
