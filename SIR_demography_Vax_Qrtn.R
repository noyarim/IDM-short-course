###############################################################################################

#Adds to the basic SIR model (SIR with interventions)

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

## SIR model with vaccination ##

# 1. Define model function
# 1-1. SIR model with vaccination 
# where vaccination provides perfect protection from infection: S -> R
OpenSIR_Vax_pp<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + I + R
    
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + birth*N - death*S - mu*S
    dI <- beta*S*I/N - death*I - gamma*I
    dR <- gamma*I + mu*S - death*R 
    # cumulative incidence
    dC <- beta*S*I/N
    # return the rates of change as a list
    list(c(dS, dI, dR, dC))
  })
}

# 1-2. SIR model with vaccination
# where vaccination provides imperfect protection from infection: change in risk of infection and no change in risk of infecting others
OpenSIR_Vax_ip1<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V
    
    #compartments without vaccination
    dS_NotV <- -beta*S_NotV*(I_NotV+I_V)/N + birth*N - death*S_NotV - mu*S_NotV
    dI_NotV <- beta*S_NotV*(I_NotV+I_V)/N - death*I_NotV - gamma*I_NotV
    dR_NotV <- gamma*I_NotV - death*R_NotV 
    #compartments with vaccination
    dS_V <- -beta*S_V*(I_NotV+I_V)/N*(1-alpha) + birth*N - death*S_V + mu*S_NotV
    dI_V <- beta*S_V*(I_NotV+I_V)/N*(1-alpha) - death*I_V - gamma*I_V
    dR_V <- gamma*I_V - death*R_V 
    
    #cumulative number of cases
    dC <- beta*S_NotV*(I_NotV+I_V)/N + beta*S_V*(I_NotV+I_V)/N*(1-alpha)
    
    # return the rates of change as a list
    list(c(dS_NotV,dS_V,dI_NotV,dI_V, dR_NotV,dR_V,dC))
  })
}

# 1-3. SIR model with vaccination
# where vaccination provides imperfect protection from infection: change in risk of infection and  in risk of infecting others
OpenSIR_Vax_ip2<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V
    beta.2 = mbeta*beta.1 # transmission rate given vaccination
    
    #compartments without vaccination
    dS_NotV <- -beta.1*S_NotV*I_NotV/N - beta.2*S_NotV*I_V/N + birth*N - death*S_NotV - mu*S_NotV
    dI_NotV <- beta.1*S_NotV*I_NotV/N + beta.2*S_NotV*I_V/N - death*I_NotV - gamma*I_NotV
    dR_NotV <- gamma*I_NotV - death*R_NotV 
    #compartments with vaccination
    dS_V <- -beta.1*S_V*I_NotV/N*(1-alpha) - beta.2*S_V*I_V/N*(1-alpha) + birth*N - death*S_V + mu*S_V
    dI_V <- beta.1*S_V*I_NotV/N*(1-alpha) + beta.2*S_V*I_V/N*(1-alpha) - death*I_V - gamma*I_V
    dR_V <- gamma*I_V - death*R_V 
    
    #cumulative number of cases
    dC <- beta.1*S_NotV*I_NotV/N + beta.1*S_V*I_NotV/N*(1-alpha) + beta.2*S_NotV*I_V/N + beta.2*S_V*I_V/N*(1-alpha)
    
    # return the rates of change as a list
    list(c(dS_NotV,dS_V,dI_NotV,dI_V, dR_NotV,dR_V, dC))
  })
}


# 1-4. SIR model with quarantine
# where a proportion of the infected is under quarantine: I -> Q
OpenSIR_Qrtn<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + I + Q + R 
    
    #SIR with quarantine
    dS<- -beta*S*I/N + birth*N - death*S
    dI <- beta*S*I/N - death*I - q*I - gamma*I
    dQ <- q*I - death*Q - gamma*Q
    dR <- gamma*I + gamma*Q - death*R
    dC <- beta*S*I/N
    
    # return the rates of change as a list
    list(c(dS, dI, dQ, dR, dC))
  })
}


#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                beta.1 = 0.5, # transmission rate given no-vaccination
                mbeta = 0.5, # this should be lower than 1
                gamma = 0.3, #recovery rate (1/duration infection)
                birth = 0.03, #birth rate (per capita)
                death = 0.03, #all-cause mortality rate
                omega = 0, # waning immunity
                mu = 0.01, # vaccination rate
                alpha = 0.3, # vaccine effectiveness
                q = 0.02 # quarantine rate
            
)

# Initial state distribution
state.pp <- c(S = 99999, #population of 100,000, 1 person starts of infected
           I = 1, 
           R = 0,
           C = 0)

state.ip <- c(S_NotV = 99999,
             S_V = 0,
             I_NotV = 1,
             I_V = 0,
             R_NotV = 0,
             R_V = 0,
             C = 0)

state.q <- c(S = 99999, #population of 100,000, 1 person starts of infected
             I = 1, 
             Q = 0,
             R = 0,
             C = 0)

T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step 

#Run the base-case
output.vax.pp <- ode(y = state.pp, times = times, func = OpenSIR_Vax_pp, parms = parameters)
output.vax.ip1 <- ode(y = state.ip, times = times, func = OpenSIR_Vax_ip1, parms = parameters)
output.vax.ip2 <- ode(y = state.ip, times = times, func = OpenSIR_Vax_ip2, parms = parameters)
output.qrtn <- ode(y = state.q, times = times, func = OpenSIR_Qrtn, parms = parameters)

## Sensitivity Analysis
# 1. Vaccination rate (on OpenSIR_Vax_ip1)
mu_list <- seq(0,0.1,length.out=3)
mu_out_all <- data.frame()

for(this_mu in mu_list){
  temp_param <- parameters
  temp_param[['mu']] <- this_mu
  this_out <- data.frame(ode(y = state.ip, times = times, func = OpenSIR_Vax_ip1, parms = temp_param))
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V,
           mu = this_mu,
           S = S_NotV + S_V,#total susceptible
           I = I_NotV + I_V, # total infected
           R = R_NotV + R_V) # total recovered
  
  if(this_mu == mu_list[1]){
    mu_out_all <- this_out
  }else{
    mu_out_all <- rbind(mu_out_all, this_out)
  }
}
mu_out_all_t <- melt(mu_out_all %>% select(time,S,I,R,N,mu), id.vars=c("time","mu","N"))
# Plot S, I, R
ggplot(mu_out_all_t) +
  geom_line(aes(x=time,y=value/N*100,color=variable))+
  facet_wrap(.~mu)+
  theme_bw()
# Plot S, I, R stratified by vaccination status
ggplot(mu_out_all)+
  geom_point(aes(x=time,y=I_NotV),color='green')+
  geom_point(aes(x=time,y=I_V),color='tomato')+
  facet_wrap(.~mu)
# Plot prevalence of I (I/N)
ggplot(mu_out_all) + 
  geom_point(aes(x=time,y=I/N))+
  facet_wrap(.~mu)+
  theme_bw()
# Cumulative number of infections
ggplot(mu_out_all) + 
  geom_point(aes(x=time,y=C))+
  facet_wrap(.~mu)+
  theme_bw()


# 2. Vaccine effectiveness
alpha_list <- seq(0,0.3,length.out=3)
alpha_out_all <- data.frame()

for(this_alpha in alpha_list){
  temp_param <- parameters
  temp_param[['alpha']] <- this_alpha
  this_out <- data.frame(ode(y = state.ip, times = times, func = OpenSIR_Vax_ip1, parms = temp_param))
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V,
           alpha = this_alpha,
           S = S_NotV + S_V,
           I = I_NotV + I_V,
           R = R_NotV + R_V)
  
  if(this_alpha == alpha_list[1]){
    alpha_out_all <- this_out
  }else{
    alpha_out_all <- rbind(alpha_out_all, this_out)
  }
}
alpha_out_all_t <- melt(alpha_out_all%>% select(time,S,I,R,N,alpha), id.vars=c("time","alpha","N"))
# Plot S, I, R
ggplot(alpha_out_all_t) +
  geom_line(aes(x=time,y=value/N*100,color=variable))+
  facet_wrap(.~alpha)+
  theme_bw()
# Plot S, I, R stratified by vaccination status
ggplot(alpha_out_all)+
  geom_point(aes(x=time,y=I_NotV),color='green')+
  geom_point(aes(x=time,y=I_V),color='tomato')+
  facet_wrap(.~alpha)
# Plot prevalence of I
ggplot(alpha_out_all) + 
  geom_point(aes(x=time,y=I/N))+
  facet_wrap(.~alpha)+
  theme_bw()
# Plot cumulative number of infections
ggplot(alpha_out_all) + 
  geom_point(aes(x=time,y=C))+
  facet_wrap(.~alpha)+
  theme_bw()
# Breakthrough infection
ggplot(alpha_out_all)+
  geom_point(aes(x=time,y=I_V))+
  facet_wrap(.~alpha)+
  theme_bw()

# 3. Ratio of effective contact rate among infected individuals with and without vaccination 
mbeta_list <- seq(1,0.7,length.out=3) 
mbeta_out_all <- data.frame()

for(this_mbeta in mbeta_list){
  temp_param <- parameters
  temp_param[['mbeta']] <- this_mbeta
  this_out <- data.frame(ode(y = state.ip, times = times, func = OpenSIR_Vax_ip2, parms = temp_param))
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V,
           mbeta = this_mbeta,
           S = S_NotV + S_V,
           I = I_NotV + I_V,
           R = R_NotV + R_V)
  
  if(this_mbeta == mbeta_list[1]){
    mbeta_out_all <- this_out
  }else{
    mbeta_out_all <- rbind(mbeta_out_all, this_out)
  }
}
mbeta_out_all_t <- melt(mbeta_out_all%>% select(time,S,I,R,N,mbeta), id.vars=c("time","mbeta","N"))
# Plot S,I,R
ggplot(mbeta_out_all_t) +
  geom_line(aes(x=time,y=value/N*100,color=variable))+
  facet_wrap(.~mbeta)+
  theme_bw()
# plot cumulative number of infections
ggplot(mbeta_out_all) + 
  geom_point(aes(x=time,y=C))+
  facet_wrap(.~mbeta)+
  theme_bw()

# 4. Quarantine rate
q_list <- c(0,0.02,0.05)
q_out_all <- data.frame()

for(this_q in q_list){
  temp_param <- parameters
  temp_param[['q']] <- this_q
  this_out <- data.frame(ode(y = state.q, times = times, func = OpenSIR_Qrtn, parms = temp_param))
  this_out <- this_out %>%
    mutate(N = S+I+R+Q,
           q = this_q)  
  if(this_q == q_list[1]){
    q_out_all <- this_out
  }else{
    q_out_all <- rbind(q_out_all, this_out)
  }
}
q_out_all_t <- melt(q_out_all%>% select(time,S,I,R,N,q), id.vars=c("time","q","N"))
# Plot S, I, R
ggplot(q_out_all_t) +
  geom_line(aes(x=time,y=value/N*100,color=variable))+
  facet_wrap(.~q)+
  theme_bw()
# Plot cumulative number of infections
ggplot(q_out_all) + 
  geom_point(aes(x=time,y=C))+
  facet_wrap(.~q)+
  theme_bw()

