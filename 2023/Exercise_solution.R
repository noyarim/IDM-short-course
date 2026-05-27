###############################################################################################

#Exercise solution

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

## Answers to the question 1 and 2 ##
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
    
    dC <- beta*S*I/N
    
    # return the rates of change as a list
    list(c(dS, dE, dI, dR, dC))
  })
}

#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.2, # place holder for effective contact rate (aka transmission rate)
                gamma = 1/14, #recovery rate (1/duration infection)
                birth = 12/1000/365,#0.0000027, #birth rate (per capita)
                death = 12/1000/365, #all-cause mortality rate
                omega = 1/(30.5*3),#0.01, # waning immunity
                t_lat = 5 # latent period from E
                
)

beta_low <- 3*(parameters[["death"]] + parameters[["gamma"]])*(parameters[["death"]] + 1/parameters[["t_lat"]])*parameters[["t_lat"]]
beta_high <- 10*(parameters[["death"]] + parameters[["gamma"]])*(parameters[["death"]] + 1/parameters[["t_lat"]])*parameters[["t_lat"]]

state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           R = 0,
           C = 0 #track cumulative number of infections
)


T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 

#Run the base-case with lower beta (R0=3)
parameters['beta'] <- beta_low
out_lowbeta <- data.frame(ode(y = state, times = times, func = OpenSEIRS, parms = parameters))
out_lowbeta$totI <- out_lowbeta$E + out_lowbeta$I
#Run the base-case with higher beta (R0=10)
parameters['beta'] <- beta_high
out_highbeta <- data.frame(ode(y = state, times = times, func = OpenSEIRS, parms = parameters))
out_highbeta$totI <- out_highbeta$E + out_highbeta$I

# Plot S,E,I,R over time
out_lowbeta_t <- melt(out_lowbeta, id.vars='time')
ggplot(out_lowbeta_t)+
  geom_line(aes(time,value,color=variable))+
  ggtitle("R0 = 3")+
  theme_bw()
max(out_lowbeta$C) #cumulative infections

out_highbeta_t <- melt(out_highbeta, id.vars='time')
ggplot(out_highbeta_t)+
  geom_line(aes(time,value,color=variable))+
  ggtitle("R0 = 10")+
  theme_bw()
max(out_highbeta$C) #cumulative infections

# Plot E + I over time & peak of epidemics
totI <- data.frame(time = out_lowbeta$time, lowbeta = out_lowbeta$totI, highbeta = out_highbeta$totI)
ggplot(totI)+
  geom_line(aes(time,lowbeta, color='low beta'))+
  geom_line(aes(time,highbeta, color='high beta'))+  
  ylab("E + I")+
  geom_vline(xintercept = c(which.max(totI$lowbeta),which.max(totI$highbeta)),
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
# Plot only I
onlyI <- data.frame(time = out_lowbeta$time, lowbeta = out_lowbeta$I, highbeta = out_highbeta$I)
ggplot(onlyI)+
  geom_line(aes(time,lowbeta, color='low beta'))+
  geom_line(aes(time,highbeta, color='high beta'))+  
  ylab("Infectious population")+
  geom_vline(xintercept = c(which.max(onlyI$lowbeta),which.max(onlyI$highbeta)),
             linetype = c("dashed","dashed"),
             color = c("#0072B2","#D55E00"))+
  annotate(geom="text",
           label = c(as.character(which.max(onlyI$lowbeta)),as.character(which.max(onlyI$highbeta))),
           x = c(which.max(onlyI$lowbeta),which.max(onlyI$highbeta)),
           y = c(max(onlyI$lowbeta), max(onlyI$highbeta)),
           color = c("#0072B2","#D55E00"),
           hjust = -0.5)+
  ggtitle("Total number of infectious individuals (I)")+
  theme_bw()+
  scale_color_manual(values = c("#D55E00","#0072B2"))


# Answers for the question #3 and #4 

# Model quarantine
#1. Define model function
# SEIRS model with quarantine
OpenSEIRS_qrtn<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + E + I + R + Q
    sigma = 1/t_lat # 1/latent period
    q_t = ifelse(t<t_start_q, 0, q) # before implementing quarantine, q(t) = 0
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + birth*N - death*S + omega*R
    dE <- beta*S*I/N - sigma*E - death*E
    dI <- sigma*E - death*I - gamma*I - q_t*I
    dQ <- q_t*I - death*Q - gamma*Q
    dR <- gamma*I + gamma*Q - death*R - omega*R
    # Cumulative number of infection
    dC <- beta*S*I/N
    # return the rates of change as a list
    list(c(dS, dE, dI, dQ, dR, dC))
  })
}

#2. Define parameters and starting compartment sizes
beta_low <- 3*(parameters[["death"]] + parameters[["gamma"]])*(parameters[["death"]] + 1/parameters[["t_lat"]])*parameters[["t_lat"]]
beta_high <- 10*(parameters[["death"]] + parameters[["gamma"]])*(parameters[["death"]] + 1/parameters[["t_lat"]])*parameters[["t_lat"]]

parameters <- c(beta = beta_low, #effective contact rate (aka transmission rate)
                gamma = 1/14, #recovery rate (1/duration infection)
                birth = 12/1000/365,#0.0000027, birth rate (per capita)
                death = 12/1000/365, #all-cause mortality rate
                omega = 1/(30.5*3),#0.01, # waning immunity
                t_lat = 5, # latent period from E to I
                q = 0.75, # quarantine rate after implementing quarantine
                t_start_q = 50 # time to implement quarantine
                
)

state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           Q = 0,
           R = 0,
           C = 0
)

T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 

# Check if quarantine is working correctly:
test_param1 <- parameters
test_param1['t_start_q'] <- 7
test1 <- ode(y = state, times=times, func=OpenSEIRS_qrtn, parms = test_param1)
test_param2 <- parameters
test_param2['t_start_q'] <- 21
test2 <- ode(y = state, times=times, func=OpenSEIRS_qrtn, parms = test_param2)
plot(test1)
plot(test2)

# Start quarantine at day = 0/21 with coverage 0.5/0.9 at different R0s
t_start_q_list <- c(0,0,21) # 0 is for abscence of intervention
q_covg_list <- c(0,0.5, 0.9) # 0 is for abscence of intervention
beta_list <- c(beta_low,beta_high)
twsa_dt <- data.frame() # data to save the ode outcomes

for (this_beta in beta_list){
  temp_param <- parameters
  # Update beta
  temp_param['beta'] = this_beta
  for (i in 1:length(t_start_q_list)){
    # update quarantine timing and coverage
    temp_param['t_start_q'] = t_start_q_list[i]
    temp_param['q'] = q_covg_list[i]
    # Run ode solver
    this_output <- as.data.frame(ode(y = state, times=times, func=OpenSEIRS_qrtn, parms = temp_param))
    # Save beta
    this_output$beta = as.character(this_beta)
    # Save t_start_q
    this_output$t_start_q = as.character(t_start_q_list[i])
    # Save q
    this_output$q_covg = as.character(q_covg_list[i])
    
    # Stack the ode result
    twsa_dt <- rbind(twsa_dt, this_output)
  }
}

# Change the label for beta, t_start_q, and q_covg
beta_lb <- sapply(beta_list,function(x) paste0("beta=",round(x, 2)))
twsa_dt$beta <- factor(twsa_dt$beta, labels = beta_lb)
tstartq_lb <- sapply(t_start_q_list,function(x) paste0("t_startq=",x))
twsa_dt$t_start_q <- factor(twsa_dt$t_start_q, levels=as.character(t_start_q_list), labels = tstartq_lb)
qcovg_lb <- sapply(q_covg_list,function(x) paste0("q_covg=",x))
twsa_dt$q_covg <- factor(twsa_dt$q_covg, levels=as.character(q_covg_list), labels = qcovg_lb)
# Check if label is created correctly
head(twsa_dt)

# Plot the two-way sensitivity analysis result
# I compartment
ggplot(twsa_dt)+
  geom_line(aes(x=time, y=I))+
  facet_grid(beta~t_start_q+q_covg,scales="free")+
  ylab("Number of infections")+
  xlab("Time")+
  theme_bw()

ggplot(twsa_dt %>% filter(q_covg!="q_covg=0"))+
  geom_line(aes(x=time, y=C))+
  facet_grid(beta~t_start_q+q_covg)+
  ylab("Cumulative Infections")+
  xlab("Time")+
  theme_bw()

twsa_dt <- twsa_dt %>% 
  mutate(q_lab=if_else(q_covg=="q_covg=0", "None",
                     if_else(q_covg=="q_covg=0.5", "Early/Low",
                             "Late/High")))
ggplot(twsa_dt %>% filter(q_covg!="q_covg=0" & beta=="beta=0.71"))+
  geom_line(aes(x=time, y=C, color=q_lab))+
  ylab("Cumulative Infections")+
  xlab("Time")+
  labs(color="Quarantine Policy") + 
  theme_bw()

