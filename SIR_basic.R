###############################################################################################

#Core/Basic SIR model code

###############################################################################################

#1. Load packages
#if you are installing packages for the first time (only run once):
#install.packages("deSolve")
#install.packages("tidyverse")

#if you've already installed the packages:
library(deSolve) #differential equation solver
library(tidyverse)  #for manipulating dataframes and graphing in R
library(cowplot) #for combining multiple graphs on the same plot

#2. Start by defining:
# A list (or vector) of model parameters (e.g. effective contact rate, recovery rate)
# A vector of initial compartment sizes (often 1 infected person and the rest susceptible)
# A vector of time steps corresponding to how long we want to run the model.

parameters <- list(beta = 0.5, #effective contact rate
                gamma = 0.3 #recovery rate (1/duration of infection)
)

state <- c(S = 99999, #population of 100,000, 1 person starts of infected
           I = 1, 
           R = 0
) #these could also be proportions (0.99 and 0.01, for instance)

T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step 

#3. Define function for a basic SIR model without demography
# This will be used with the deSolve package to simulate how your population moves between compartments over time

BasicSIR<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{ 
    N = S + I + R #define N (total population size)
    
    #SIR model equations from slides - rates of change in and out of each compartment 
    dS <- -beta*S*I/N
    dI <- beta*S*I/N - gamma*I
    dR <- gamma*I
    
    #return the rates of change as a list - IMPORTANT: must be same order as "state"
    list(c(dS, dI, dR)) 
  })
}

#4. Run the SIR model using the ode() function that is part of the deSolve package.
output <- ode(y = state, times = times, func = BasicSIR, parms = parameters)

#Note: ode() is a function in the deSolve package. 
# You must fill it in using the syntax: 
#    "ode(y= [vector of initial compartment sizes], 
#         times = [vector of time steps], 
#         func = [name of SIR model function], 
#         parms = [vector or list of parameter values])"  


#View and analyze model output
print(head(output, 10)) #I's are rapidly increasing
print(tail(output, 10)) #steady state - state sizes aren't changing

# 5. Plot results. We define a function that will plot the epidemic curve of our model output from step #4. 
output <- as.data.frame(output) %>% mutate(N=S+I+R)
output_long <- pivot_longer(output, cols=c("S","I","R"), names_to="state", values_to="size")

plot_trace <-function(out) {
  fig <- ggplot(out, aes(x=time, y=size/N, color=state)) + 
    geom_line(linewidth=1.25) +
    labs(x='Time', y='Compartment proportions', color='') +
    theme_bw() + theme(panel.grid=element_blank())
  return(fig)
}

plot_trace(output_long %>% filter(time <= 200))

#6. calculate basic reproductive number R0 (doesn't vary over time)
R0 <- parameters[["beta"]]/parameters[["gamma"]]
print(R0)
#calculate effective reproductive number Rt (declines over time as % population susceptible declines)
output <- output %>% mutate(Rt=R0*S/N)

#7.  Let's further examine how infection dynamics and Rt are related
fig1 <- ggplot(output %>% filter(time<=200)) + 
  geom_line(aes(x=time, y=I/N)) + 
  labs(x="Time", y="% infected") +
  theme_bw() + theme(panel.grid=element_blank())
fig2 <- ggplot(output %>% filter(time<=200)) + 
  geom_line(aes(x=time, y=Rt)) +
  geom_hline(yintercept=1, linetype="dashed", color="red") +
  labs(x="Time", y="R(t)") +
  theme_bw() + theme(panel.grid=element_blank())
fig <- plot_grid(fig1, fig2, nrow=2)
fig

#8. Finally, let's see how changing model parameters affects infection dynamics.
# We'll do this w/ a simple 2-way sensitivity analysis, changing beta (effective contact rate) and gamma (recovery rate)
output_all <- list()
for(betas in c(0.2, 0.5, 1)) {
  for(gammas in c(0.1, 0.3, 0.5)) {
    parameters <- c("beta"=betas, 
                    "gamma"=gammas)
    output <- ode(y = state, times = times, func = BasicSIR, parms = parameters)
    output_long <- pivot_longer(as.data.frame(output) %>% mutate(N=S+I+R), 
                                cols=c("S","I","R"), names_to="state", values_to="size")
    output_all <- c(output_all, list(output_long %>% mutate(beta=paste0("Beta: ", betas), 
                                                            gamma=paste0("Gamma: ", gammas),
                                                            R0=paste0("R0: ", round(betas/gammas, 2)))))
  }
}

output_all <- bind_rows(output_all)
output_all <- output_all %>% mutate(label=paste(beta, gamma, R0, sep="; "))

plot_trace(output_all) + facet_wrap(~label)

