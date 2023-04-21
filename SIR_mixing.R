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

parameters <- list(contact_matrix=matrix(data=c(29, 1, 1, 5)*0.1, 
                                         nrow=2, ncol=2, byrow=T), 
                   # before: 5 contacts with 10% probability of infection per infected contact
                   # now: still 10% probability of infection per infected contact, but # contacts varies by risk group (still avg of 5 contacts per person)
                   # WAIFW: effective contact rate to h from h, to l from h, to h from l, to l from l
                   # here, the high-risk population is high-risk because they have more contacts with each other - which will increase their prevalence/incidence
                   # could represent something like sexual partners and STIs
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
      dS <- -1*state[, "S"]*state[,"I"]%*%contact_matrix/sum(state) 
      dI <- state[, "S"]*state[,"I"]%*%contact_matrix/sum(state) -
        gamma*state[,"I"]
      dR <- gamma*state[,"I"]
      
      #return the rates of change as a list
      list(c(dS, dI, dR)) 
    })
    
}

#4. Run the SIR model using the ode() function that is part of the deSolve package.
output <- ode(y = state, times = times, func = MixingSIR, parms = parameters)
output <- as.data.frame(output) 

#5. View and analyze model output
output_long <- pivot_longer(output, cols=2:ncol(output), names_to=c("state", "risk"), names_sep="_", values_to="size")
output_long <- output_long %>% group_by(time, risk) %>% mutate(N=sum(size))

plot_trace_risk <-function(out) {
  fig <- ggplot(out, aes(x=time, y=size/N, color=state, linetype=risk)) + 
    geom_line(linewidth=1.25) +
    labs(x='Time', y='Compartment size', color='', linetype='') +
    theme_bw() + theme(panel.grid=element_blank())
  return(fig)
}

plot_trace_risk(output_long %>% filter(time<100))

#in this example, almost all the high-risk group get infected but many of the low-risk group do not

#6. Compare total infections to a version without risk stratification

#rerun model from SIR_basic (at the beginning of the code demos)
parameters_basic <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                      gamma = 0.3 #recovery rate (1/duration of infection)
)

state_basic <- c(S = 99999, #population of 100,000, 1 person starts of infected
                 I = 1, 
                 R = 0
)

BasicSIR<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{ 
    
    N = S + I + R #define N (total population size)
    
    #SIR model equations from lecture - rates of change in and out of each compartment 
    dS <- -beta*S*I/N
    dI <- beta*S*I/N - gamma*I
    dR <- gamma*I
    
    #return the rates of change as a list
    list(c(dS, dI, dR)) 
  })
}

plot_trace <-function(out) {
  fig <- ggplot(out, aes(x=time, y=size/N, color=state)) + 
    geom_line(linewidth=1.25) +
    labs(x='Time', y='Compartment size', color='') +
    theme_bw() + theme(panel.grid=element_blank())
  return(fig)
}

output_basic <- ode(y = state_basic, times = times, func = BasicSIR, parms = parameters_basic)
output_basic <- as.data.frame(output_basic) %>% mutate(N=S+I+R)
output_long_basic <- pivot_longer(output_basic, cols=c("S","I","R"), names_to="state", values_to="size")
plot_trace(output_long_basic %>% filter(time<100))

#sum across risk groups
output_long_sum <- output_long %>% group_by(time, state) %>% 
  summarise(size=sum(size), N=sum(N))
plot_trace(output_long_sum %>% filter(time<100))

# compare these 2 figures
# without stratified risk groups, epidemic moves more slowly (lower R0) and cumulatively lower % are infected 

#7. Let's test this out by calculating R0 and Rt
R_matrix <- parameters$contact_matrix*(rowSums(matrix(data=state, nrow=ncol(parameters$contact_matrix)))/sum(state))/
  (parameters$gamma)
R0 <- max(eigen(R_matrix)$values) #this is the R0 of the stratified version

R0_basic <- parameters_basic[["beta"]]/parameters_basic[["gamma"]]
print(paste0("R0 basic: ", round(R0_basic, 2), "; RO stratified: ", round(R0, 2)))

#with Rt
if(FALSE) {
  output_basic <- output_basic %>% mutate(Rt=R0_basic*S/N)
  output_long_sum <- output_long_sum %>% mutate(Rt=R0*size[state=="S"]/N)
  
  ggplot() +
    geom_line(data=output_long_sum %>% filter(time<100), aes(x=time, y=Rt, color="Stratified")) +
    geom_line(data=output_basic %>% filter(time<100), aes(x=time, y=Rt, color="Basic")) +
    labs(x="Time", y="R_t", color="Model version") + theme_bw()
}


#8. Changing mixing patterns - but keeping total # contacts the same
contact_matrices <- list(matrix(data=c(35, 0, 0, 5)*0.1,
                             nrow=2, ncol=2, byrow=T),
                      matrix(data=c(29, 1, 1, 5)*0.1, 
                             nrow=2, ncol=2, byrow=T),
                      matrix(data=c(17, 3, 3, 5)*0.1,
                             nrow=2, ncol=2, byrow=T),
                      matrix(data=c(5,5,5,5)*0.1,
                             nrow=2, ncol=2, byrow=T)
                      )

names(contact_matrices) <- c("No mixing", "More assortative mixing", 
                             "Less assortative mixing", "Homogeneous mixing")

output_all <- list()
for(i in names(contact_matrices)) {
    parameters[["contact_matrix"]] <- contact_matrices[[i]]
    output <- ode(y = state, times = times, func = MixingSIR, parms = parameters)
    output_long <- pivot_longer(as.data.frame(output), cols=2:ncol(output), 
                                names_to=c("state", "risk"), names_sep="_", values_to="size")
    output_long <- output_long %>% group_by(time, risk) %>% mutate(N=sum(size))
    output_all <- c(output_all, list(output_long %>% mutate(contact_pattern=i)))
}

output_all <- bind_rows(output_all) %>%
  mutate(contact_pattern=factor(contact_pattern, levels=names(contact_matrices)))

plot_trace_risk(output_all %>% filter(time<=100)) + facet_wrap(~contact_pattern)
#more homogeneous mixing -> more of the low-risk group gets infected, slower dynamics

output_all_sum <- output_all %>% group_by(time, state, contact_pattern) %>% 
  summarise(size=sum(size), N=sum(N))
plot_trace(output_all_sum %>% filter(time<100)) + facet_wrap(~contact_pattern)

#check that homogeneous mixing version is same as basic model without risk stratification
plot_trace(output_all_sum %>% filter(time<200 & contact_pattern=="Homogeneous mixing"))
plot_trace(output_long_basic %>% filter(time<200))

