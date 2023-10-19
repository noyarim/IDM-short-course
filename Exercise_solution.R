###############################################################################################

#Exercise solution

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)

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
    
    # return the rates of change as a list
    list(c(dS, dE, dI, dR))
  })
}

#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.2, #effective contact rate (aka transmission rate)
                gamma = 1/14, #recovery rate (1/duration infection)
                birth = 12/1000/365,#0.0000027, #birth rate (per capita)
                death = 12/1000/365, #all-cause mortality rate
                omega = 1/(30.5*3),#0.01, # waning immunity
                t_lat = 5 # latent period from E
                
)

state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           R = 0
)


T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 

#Run the base-case with beta = 0.2
parameters['beta'] <- 0.2
out_lowbeta <- data.frame(ode(y = state, times = times, func = OpenSEIRS, parms = parameters))
out_lowbeta$totI <- out_lowbeta$E + out_lowbeta$I
#Run the base-case with beta = 0.8
parameters['beta'] <- 0.8
out_highbeta <- data.frame(ode(y = state, times = times, func = OpenSEIRS, parms = parameters))
out_highbeta$totI <- out_highbeta$E + out_highbeta$I

# Plot S,E,I,R over time
out_lowbeta_t <- melt(out_lowbeta, id.vars='time')
ggplot(out_lowbeta_t)+
  geom_line(aes(time,value,color=variable))+
  ggtitle("Beta = 0.2")+
  theme_bw()

out_highbeta_t <- melt(out_highbeta, id.vars='time')
ggplot(out_highbeta_t)+
  geom_line(aes(time,value,color=variable))+
  ggtitle("Beta = 0.8")+
  theme_bw()

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
  
  
# Model quarantine
#1. Define model function
# SEIRS model with qaurantine
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
    
    # return the rates of change as a list
    list(c(dS, dE, dI, dQ, dR))
  })
}

#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.2, #effective contact rate (aka transmission rate)
                gamma = 1/14, #recovery rate (1/duration infection)
                birth = 12/1000/365,#0.0000027, birth rate (per capita)
                death = 12/1000/365, #all-cause mortality rate
                omega = 1/(30.5*3),#0.01, # waning immunity
                t_lat = 5, # latent period from E to I
                q = 0.1, # quarantine rate after implementing quarantine
                t_start_q = 50 # time to implement quarantine
                
)

state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           Q = 0,
           R = 0
)

T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 

# Check if quarantine is working correctly:
test_param1 <- parameters
test_param1['t_start_q'] <- 30
test1 <- ode(y = state, times=times, func=OpenSEIRS_qrtn, parms = test_param1)
test_param2 <- parameters
test_param2['t_start_q'] <- 150
test2 <- ode(y = state, times=times, func=OpenSEIRS_qrtn, parms = test_param2)
plot(test1)
plot(test2)

# Start quarantine at day = 10/100 and beta = 0.2/0.8
t_start_q_list <- c(30,150)
beta_list <- c(0.2,0.8)
twsa_dt <- data.frame() # data to save the ode outcomes

for (this_beta in beta_list){
  temp_param <- parameters
  # Update beta
  temp_param['beta'] = this_beta
  for (this_tstart_q in t_start_q_list){
    # update omega
    temp_param['t_start_q'] = this_tstart_q
    # Run ode solver
    this_output <- as.data.frame(ode(y = state, times=times, func=OpenSEIRS_qrtn, parms = temp_param))
    # Save beta
    this_output$beta = as.character(this_beta)
    # Save t_start_q
    this_output$t_start_q = as.character(this_tstart_q)
    # Stack the ode result
    twsa_dt <- rbind(twsa_dt, this_output)
  }
}

# Change the label for beta and t_start_q
beta_lb <- sapply(beta_list,function(x) paste0("beta=",x))
twsa_dt$beta <- factor(twsa_dt$beta, labels = beta_lb)
tstartq_lb <- sapply(t_start_q_list,function(x) paste0("t_startq=",x))
twsa_dt$t_start_q <- factor(twsa_dt$t_start_q, levels=c("30","150"), labels = tstartq_lb)
# Check if label is created correctly
head(twsa_dt)

# Plot the two-way sensitivity analysis result
# I compartment
ggplot(twsa_dt)+
  geom_line(aes(x=time, y=I))+
  facet_grid(beta~t_start_q)+
  ylab("Infected")+
  xlab("Time")+
  theme_bw()

# S, I, R, Q compartment
twsa_dt_t <- melt(twsa_dt, id.vars=c("time","beta","t_start_q"))

ggplot(twsa_dt_t)+
  geom_line(aes(time,value,color=variable))+
  facet_grid(beta~t_start_q)+
  theme_bw()

# low & early quarantine vs. high & late quarantine
# 1. early implementation of quarantine at low rate
param.q1 <- parameters
param.q1['t_start_q']<-10
param.q1['q'] <- 0.1
out.q1 <- as.data.frame(ode(y = state, times=times, func=OpenSEIRS_qrtn, parms = param.q1))

param.q1 <- parameters
param.q1['t_start_q']<-10
param.q1['q'] <- 0.1
out.q1 <- as.data.frame(ode(y = state, times=times, func=OpenSEIRS_qrtn, parms = param.q1))





