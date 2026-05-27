###############################################################################################

#SIR model with births and deaths

###############################################################################################

library(deSolve) #differential equation solver
library(tidyverse)

#Load plotting function - same as SIR_basic
plot_trace <-function(out) {
  fig <- ggplot(out, aes(x=time, y=size/N, color=state)) +
    geom_line(linewidth=1.25) +
    labs(x='Time', y='Compartment proportions', color='') +
    theme_bw() + theme(panel.grid=element_blank())
  return(fig)
}

#1. Define the main model function
# It looks the same as BasicSIR from the last R lab, but with births and deaths added
OpenSIR<-function(t, state, params) {
  with(as.list(c(state, params)),{
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
params <- list(beta = 0.5, #effective contact rate
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

#3. Run model and view output
output <- ode(y = state, times = times, func = OpenSIR, parms = params)
output <- as.data.frame(output) %>% mutate(N=S+I+R)
output_long <- pivot_longer(output, cols=c("S","I","R"), names_to="state", values_to="size")
plot_trace(output_long)
output[501, c("S", "I", "R")]/100000

#4. We didn't go over this in the slides, but in this model there are equations for the steady state distribution
# Calculate R0 and then the steady state/equilibrium distributions
R0 <- params[["beta"]]/(params[["gamma"]] + params[["death"]])
print(paste0("R0: ", R0))

# if R0 <= 1, there is no endemic equilibrium
# if R0 > 1, this is the endemic equilibrium:
S_star <- 1/R0
I_star <- (params[["death"]]/params[["beta"]])*(R0-1)
R_star <- 1-(S_star + I_star)

#see how this matches the proportions in the graph/from the model
print(paste0("Estimated prevalence at equilibrium: ", round(100*output$I[[T_end+1]]/output$N[[T_end+1]], 2), "%"))

print("Simulated prevalences at equilibrium:")
print(c("S_star"=S_star*100, "I_star"=I_star*100, "R_star"=R_star*100)) #prevalences at equilibrium

#5. Calculate Rt and assess its relationship to infection dynamics
output <- output %>% mutate(Rt=R0*S/N)

#changes in Rt correlate w/ the oscillatory pattern in the trace graph
#changes in Rt are fueled by rises and falls in the susceptible population (from births and new infections, respectively)
plot_trace(output_long) +
  geom_line(data=output, aes(x=time, y=Rt), color="black") +
  geom_hline(yintercept=1, linetype="dashed", color="red")

#6. Use 1-way sensitivity analysis to assess how dynamics change as birth/death rates rise/fall
output_all <- list()
for(rates in c(0, 0.01, 0.03, 0.1)) {
    params[["birth"]] <- params[["death"]] <- rates
    output <- ode(y = state, times = times, func = OpenSIR, parms = params)
    output <- as.data.frame(output) %>%
      mutate(N=S+I+R,
             demo_rate=rates,
             R0=params[["beta"]]/(params[["gamma"]] + rates),
             Rt=R0*S/N,
             demo_lab=paste0("Birth/death rate: ", demo_rate),
             R0_lab=paste0("R0: ", round(R0, 2)),
             I_star_lab=paste0("Prevalence at t=500: ", round(100*I[time==500]/N[time==500], 1)))
    output_long <- pivot_longer(output, cols=c("S","I","R"), names_to="state", values_to="size")
    output_all <- c(output_all, list(output_long))
}

output_all <- bind_rows(output_all)

plot_trace(output_all) + facet_wrap(~demo_lab+R0_lab+I_star_lab)

#7. How would you modify the code or parameters to simulate a chronic disease (with demography)?
params <- list(beta = 0.5, #effective contact rate
                 gamma = 0, #recovery rate (1/duration infection)
                 birth = 0.02, #birth rate (per capita)
                 death = 0.02 #all-cause mortality rate
)

output <- ode(y = state, times = times, func = OpenSIR, parms = params)
output <- as.data.frame(output) %>% mutate(N=S+I+R)
output_long <- pivot_longer(output, cols=c("S","I","R"), names_to="state", values_to="size")
plot_trace(output_long)
