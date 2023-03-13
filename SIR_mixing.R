###############################################################################################

#Core/Basic SIR model code

###############################################################################################

#1. Load packages
#if you are installing packages for the first time (only run once):
#install.packages("deSolve")
#install.packages("ggplot2")

#if you've already installed the packages:
library(deSolve) #differential equation solver
library(tidyverse) 

#2. Start by defining:
# A vector or list of model parameters (e.g. effective contact rate, recovery rate)
# A vector of initial compartment sizes (often 1 infected person and the rest susceptible)
# A vector of time steps corresponding to how long we want to run the model.

parameters <- list(contact_matrix=matrix(data=c(5, 0.5, 0.5, 1)*0.5, #0.5 is original beta
                                         nrow=2, ncol=2, byrow=T), #WAIFW: beta_hh, beta_lh, beta_hl, beta_ll
                   gamma = 0.3 #recovery rate (1/duration infection)
)

state <- c(S_High = 24999, #high (h) and low (l) risk groups
           S_Low = 75000,
           I_High = 1, 
           I_Low = 0,
           R_High = 0,
           R_Low = 0
)

T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step 

#3. Define function for a basic SIR model without demography
# This will be used with the deSolve package to simulate how your population moves between compartments over time

MixingSIR<-function(t, state, parameters) {

    state <- matrix(data=state, nrow=ncol(parameters$contact_matrix))
    colnames(state) <- c("S", "I", "R")
    rownames(state) <- c("h", "l")
    
    with(parameters, {
      dS <- -1*state[, "S"]*state[,"I"]%*%contact_matrix/rowSums(state)
      dI <- state[, "S"]*state[,"I"]%*%contact_matrix/rowSums(state) -
        gamma*state[,"I"]
      dR <- gamma*state[,"I"]
      
      #return the rates of change as a list
      list(c(dS, dI, dR)) 
    })
    
}

#4. Run the SIR model using the ode() function that is part of the deSolve package.
output <- ode(y = state, times = times, func = MixingSIR, parms = parameters)

#Note: ode() is a function in the deSolve package. 
# You must fill it in using the syntax: 
#    "ode(y= [vector of initial compartment sizes], 
#         times = [vector of time steps], 
#         func = [name of SIR model function], 
#         parms = [vector or list of parameter values])"  


#5. View and analyze model output
print(head(output)) #I's are rapidly increasing
print(tail(output)) #steady state - state sizes aren't changing

# We define a function that will plot the epidemic curve of our model output from step #4. 
output_long <- pivot_longer(as.data.frame(output), cols=2:ncol(output), names_to="state", values_to="size")
output_long <- output_long %>% separate(state, c("state", "risk"), "_")
  
plot_trace <-function(out) {
  fig <- ggplot(out, aes(x=time, y=size, color=state, linetype=risk)) + 
    geom_line(linewidth=1.25) +
    labs(x='Time', y='Compartment size', color='', linetype='') +
    theme_bw() + theme(panel.grid=element_blank())
  return(fig)
}

plot_trace(output_long %>% filter(time<100))

#HAVEN'T UPDATED PAST HERE YET
#6. calculate R0
R0 <- parameters[["beta"]]/parameters[["gamma"]]
print(R0)


#7. Changing mixing patterns
output_all <- list()
for(betas in c(0.2, 0.5, 1)) {
  for(gammas in c(0.1, 0.3, 0.5)) {
    parameters <- c("beta"=betas, 
                    "gamma"=gammas)
    output <- ode(y = state, times = times, func = BasicSIR, parms = parameters)
    output_long <- pivot_longer(as.data.frame(output), cols=2:ncol(output), names_to="state", values_to="size")
    output_all <- c(output_all, list(output_long %>% mutate(beta=paste0("Beta: ", betas), 
                                                            gamma=paste0("Gamma: ", gammas),
                                                            R0=paste0("R0: ", round(betas/gammas, 2)))))
  }
}

output_all <- bind_rows(output_all)

plot_trace(output_all) + facet_wrap(~beta+gamma+R0)

