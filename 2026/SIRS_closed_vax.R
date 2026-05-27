###############################################################################################

#SIRS model with vaccination (closed population, no demography, with waning immunity)

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

## A. Vaccination provides perfect protection from infection: S -> R ##
# Direct flow from S to R indicates those who are vaccinated and protected from infection , with the same rate of immunological waning as natural infection.
# Backward flow R -> S at rate omega reflects waning of immunity (both natural and vaccine-induced).

#1. Define model function
SIRS_Vax_pp<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S + I + R

    #SIRS equations from lecture
    dS <- -beta*S*I/N - mu*S + omega*R
    dI <- beta*S*I/N - gamma*I
    dR <- gamma*I + mu*S - omega*R
    # cumulative incidence
    dC <- beta*S*I/N
    # return the rates of change as a list
    list(c(dS, dI, dR, dC))
  })
}

#2. Define parameters and starting compartment sizes
# Mu is vaccination rate. Omega is rate of waning immunity (1/duration of immunity).
params <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                mu = 0.01, # vaccination rate
                omega = 0.01 # waning immunity rate
)

# Initial state distribution
state.pp <- c(S = 99999, #population of 100,000, 1 person starts of infected
              I = 1,
              R = 0,
              C = 0)

#3. Run the SIRS model with perfect protection from vaccination
T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step

# Run the base-case
output.vax.pp <- ode(y = state.pp, times = times, func = SIRS_Vax_pp, parms = params)

# Plot the base-case result
plot(output.vax.pp)

# Compare the size of I compartment when mu = 0 (no vax) and mu = 0.01 (vax)
# Copy parameters
temp_param <- params
# Set vaccination rate = 0
temp_param['mu'] = 0
# Run ode solver and save the outcome in temporary output variable
temp_output <-  ode(y = state.pp, times = times, func = SIRS_Vax_pp, parms = temp_param)
# Create a dataset to compare no_vax and vax output
I_dt <- data.frame(time = temp_output[,'time'], no_vax = temp_output[,'I'], vax = output.vax.pp[,'I'])
# Plot the dataset
ggplot(I_dt)+
  geom_line(aes(time,vax,color='vax'))+
  geom_line(aes(time,no_vax,color='No vax'))+
  theme_bw()+
  scale_color_manual(values=c("blue","tomato"))

## B. Vaccination provides imperfect protection from infection ##
# Vaccination usually doesn't provide perfect protection. Here, we stratify susceptibles by vaccine status.
# Total I who can infect others = sum of infections across both groups.
# (1-alpha) indicates reduced susceptibility for vaccinated susceptibles; higher alpha -> lower susceptibility.
# This structure also allows to vary immunological waning by vaccination vs. natural infection.
# A single R compartment is shared by natural- and breakthrough-recoveries (vaccine status is lost on recovery).
# Waning immunity from natural infection: R -> S_NotV at rate omega.
# Waning vaccine-induced immunity: S_V -> S_NotV at rate omega_vx.
# We set omega_vx = omega in the params here, but they can be varied independently.

# B-1. No change in infectiousness
# In ip1, vaccinated and unvaccinated infected individuals are equally infectious, so a single
# I compartment suffices - the dynamics are identical with or without splitting I by vaccine
# status. We don't report breakthrough infections separately here; see ip2 below for that.
#1. Define model function
SIRS_Vax_ip1<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S_NotV + S_V + I + R

    #susceptible compartments
    dS_NotV <- -beta*S_NotV*I/N - mu*S_NotV + omega_vx*S_V + omega*R
    dS_V <- -beta*S_V*I/N*(1-alpha) + mu*S_NotV - omega_vx*S_V
    #single infected compartment (vaccinated and unvaccinated infected are equally infectious)
    dI <- beta*S_NotV*I/N + beta*S_V*I/N*(1-alpha) - gamma*I
    #single recovered compartment (both natural and breakthrough recoveries flow here)
    dR <- gamma*I - omega*R

    #cumulative number of cases (both naive and breakthrough infections)
    dC <- beta*S_NotV*I/N + beta*S_V*I/N*(1-alpha)

    # return the rates of change as a list
    list(c(dS_NotV, dS_V, dI, dR, dC))
  })
}

# B-2. Change in risk of infecting others with vaccination
# Apply separate beta values (beta.1 and beta.2) to I_NotV and I_V.
# Beta.2 is calculated as a fraction (mbeta) of beta.1.
# Because infectiousness now depends on vaccine status, we keep I_NotV and I_V as separate
# compartments. As a bonus, this structure lets us report breakthrough infections (I_V over time).
SIRS_Vax_ip2<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S_NotV + S_V + I_NotV + I_V + R
    beta.2 = mbeta*beta.1 # transmission rate given vaccination
    #compartments without vaccination
    dS_NotV <- -beta.1*S_NotV*I_NotV/N - beta.2*S_NotV*I_V/N - mu*S_NotV + omega_vx*S_V + omega*R
    dI_NotV <- beta.1*S_NotV*I_NotV/N + beta.2*S_NotV*I_V/N - gamma*I_NotV
    #compartments with vaccination
    dS_V <- -beta.1*S_V*I_NotV/N*(1-alpha) - beta.2*S_V*I_V/N*(1-alpha) + mu*S_NotV - omega_vx*S_V
    dI_V <- beta.1*S_V*I_NotV/N*(1-alpha) + beta.2*S_V*I_V/N*(1-alpha) - gamma*I_V
    #single recovered compartment (both natural and breakthrough recoveries flow here)
    dR <- gamma*I_NotV + gamma*I_V - omega*R

    #cumulative number of cases
    dC <- beta.1*S_NotV*I_NotV/N + beta.1*S_V*I_NotV/N*(1-alpha) + beta.2*S_NotV*I_V/N + beta.2*S_V*I_V/N*(1-alpha)

    # return the rates of change as a list
    list(c(dS_NotV,dS_V,dI_NotV,dI_V, dR, dC))
  })
}

#2. Define parameters and starting compartment sizes
# Mu is vaccination rate. Omega is rate of waning of natural-infection-induced immunity.
# Omega_vx is rate of waning of vaccine-induced immunity (set equal to omega by default).
# mbeta is 1 if no change in infectiousness with vaccination, or less than 1 if reduced infectiousness.
# alpha is vaccine effectiveness.
params <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                beta.1 = 0.5, # transmission rate given no-vaccination
                mbeta = 0.5, # 1 if no change in infectiousness with vaccination,
                           # less than 1 if reduced infectiousness with vaccination
                gamma = 0.3, #recovery rate (1/duration infection)
                mu = 0.01, # vaccination rate
                alpha = 0.3, # vaccine effectiveness
                omega = 0.01, # waning rate of natural-infection-induced immunity
                omega_vx = 0.01 # waning rate of vaccine-induced immunity
)

# Initial state vectors
# ip1 uses a single I compartment; ip2 splits I by vaccine status
state.ip1 <- c(S_NotV = 99999,
               S_V = 0,
               I = 1,
               R = 0,
               C = 0)
state.ip2 <- c(S_NotV = 99999,
               S_V = 0,
               I_NotV = 1,
               I_V = 0,
               R = 0,
               C = 0)

#3. Run the SIRS model with imperfect protection from vaccination
# Each model uses its corresponding state vector. We compute total S in both;
# for ip2 we also sum I_V + I_NotV to get total I.
T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step

# Run ODE solver
output.vax.ip1 <- ode(y = state.ip1, times = times, func = SIRS_Vax_ip1, parms = params)
output.vax.ip2 <- ode(y = state.ip2, times = times, func = SIRS_Vax_ip2, parms = params)
# Calculate total S (I is already a single column for ip1; for ip2, also calculate total I)
output.vax.ip1 <- as.data.frame(output.vax.ip1) %>%
  mutate(S = S_V + S_NotV)
output.vax.ip2 <- as.data.frame(output.vax.ip2) %>%
  mutate(S = S_V + S_NotV,
         I = I_V + I_NotV)

# Plot base-case outcomes and compare between the two models
g1<- ggplot(output.vax.ip1)+
        geom_line(aes(time,S,color="S"))+
        geom_line(aes(time,I,color="I"))+
        geom_line(aes(time,R,color="R"))+
        scale_color_manual(values=c("tomato","blue","green"))+
        ggtitle("No Change in Infectiousness")+
        ylab("Number of individuals")+
        theme_bw()
g2<- ggplot(output.vax.ip2)+
        geom_line(aes(time,S,color="S"))+
        geom_line(aes(time,I,color="I"))+
        geom_line(aes(time,R,color="R"))+
        scale_color_manual(values=c("tomato","blue","green"))+
        ggtitle("Reduced infectiousness with vaccination")+
        ylab("Number of individuals")+
       theme_bw()
# side-by-side plot
ggarrange(g1,g2,common.legend = TRUE)


## Sensitivity Analysis ##
# Perform sensitivity analysis on (1) vaccination rate, (2) vaccine effectiveness,
# (3) the ratio in effective contact rate between vaccinated and non-vaccinated individuals,
# and (4) the rate of waning immunity (varying omega and omega_vx together).
# (1), (2), and (4) use SIRS_Vax_ip1; (3) uses SIRS_Vax_ip2.

# 1. Vaccination rate (on SIRS_Vax_ip1)
# Discuss how changing vaccination rate changes the projected epidemic.
mu_list <- seq(0,0.1,length.out=3) # a vector of vaccination rate
mu_out_all <- data.frame() # an empty dataset to save outputs

for(this_mu in mu_list){
  temp_param <- params # copy parameters
  temp_param[['mu']] <- this_mu # replace mu
  this_out <- data.frame(ode(y = state.ip1, times = times, func = SIRS_Vax_ip1, parms = temp_param)) # run ode solver
  # compute total S and save mu in this cycle
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I + R,
           mu = this_mu,
           S = S_NotV + S_V) # total susceptible
  # stack the result
  if(this_mu == mu_list[1]){
    mu_out_all <- this_out
  }else{
    mu_out_all <- rbind(mu_out_all, this_out)
  }
}
# Change the output data to long-form data
mu_out_all_t <- melt(mu_out_all %>% select(time,S,I,R,N,mu), id.vars=c("time","mu","N"))
# Plot S, I, R
ggplot(mu_out_all_t) +
  geom_line(aes(x=time,y=value/N*100,color=variable))+
  ylab("Percentage of population")+
  facet_wrap(.~mu)+
  theme_bw()
# Plot prevalence of I (I/N)
ggplot(mu_out_all) +
  geom_point(aes(x=time,y=I/N))+
  facet_wrap(.~mu)+
  ylab("Proportion of population who are infected")+
  theme_bw()
# Cumulative number of infections
ggplot(mu_out_all) +
  geom_point(aes(x=time,y=C))+
  facet_wrap(.~mu)+
  ylab("Cumulative infections")+
  theme_bw()


# 2. Vaccine effectiveness
# How does increasing vaccine effectiveness change the projected epidemic?
alpha_list <- seq(0,0.3,length.out=3) # a vector of vaccine effectiveness
alpha_out_all <- data.frame() # an empty dataset to save outputs

for(this_alpha in alpha_list){
  temp_param <- params # copy parameters
  temp_param[['alpha']] <- this_alpha # replace alpha
  this_out <- data.frame(ode(y = state.ip1, times = times, func = SIRS_Vax_ip1, parms = temp_param)) # run ode solver
  # compute total S and save alpha in this cycle
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I + R,
           alpha = this_alpha,
           S = S_NotV + S_V)
  # Stack the result
  if(this_alpha == alpha_list[1]){
    alpha_out_all <- this_out
  }else{
    alpha_out_all <- rbind(alpha_out_all, this_out)
  }
}
# Change the outcome data to long-form data
alpha_out_all_t <- melt(alpha_out_all%>% select(time,S,I,R,N,alpha), id.vars=c("time","alpha","N"))
# Plot S, I, R
ggplot(alpha_out_all_t) +
  geom_line(aes(x=time,y=value/N*100,color=variable))+
  ylab("Percentage of population")+
  facet_wrap(.~alpha)+
  theme_bw()
# Plot prevalence of I
ggplot(alpha_out_all) +
  geom_point(aes(x=time,y=I/N))+
  facet_wrap(.~alpha)+
  ylab("Prevalence of infection")+
  theme_bw()
# Plot cumulative number of infections
ggplot(alpha_out_all) +
  geom_point(aes(x=time,y=C))+
  facet_wrap(.~alpha)+
  ylab("Cumulative infections")+
  theme_bw()


# 3. Ratio of effective contact rate among infected individuals with and without vaccination
# How does the projected epidemic change as we decrease the ratio of effective contact rate?
# This sensitivity uses SIRS_Vax_ip2, whose stratified I compartments let us also report
# breakthrough infections (I_V over time).
mbeta_list <- seq(1,0.7,length.out=3) # a vector of mbeta
mbeta_out_all <- data.frame() # an empty dataset to save outputs

for(this_mbeta in mbeta_list){
  temp_param <- params # copy parameters
  temp_param[['mbeta']] <- this_mbeta # replace mbeta
  this_out <- data.frame(ode(y = state.ip2, times = times, func = SIRS_Vax_ip2, parms = temp_param)) # run ode solver
  # compute total S, I and save mbeta in this cycle
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I_NotV + I_V + R,
           mbeta = this_mbeta,
           S = S_NotV + S_V,
           I = I_NotV + I_V)
  # Stack the result
  if(this_mbeta == mbeta_list[1]){
    mbeta_out_all <- this_out
  }else{
    mbeta_out_all <- rbind(mbeta_out_all, this_out)
  }
}
# change the data to long-form data
mbeta_out_all_t <- melt(mbeta_out_all%>% select(time,S,I,R,N,mbeta), id.vars=c("time","mbeta","N"))
# Plot S,I,R
ggplot(mbeta_out_all_t) +
  geom_line(aes(x=time,y=value/N*100,color=variable))+
  facet_wrap(.~mbeta)+
  ylab("Proportion of population")+
  theme_bw()
# plot cumulative number of infections
ggplot(mbeta_out_all) +
  geom_point(aes(x=time,y=C))+
  facet_wrap(.~mbeta)+
  ylab("Cumulative infections")+
  theme_bw()
# Breakthrough infections (only available with the ip2 structure)
ggplot(mbeta_out_all)+
  geom_point(aes(x=time,y=I_V))+
  facet_wrap(.~mbeta)+
  ylab("Number of breakthrough infections")+
  theme_bw()


# 4. Waning immunity rate (on SIRS_Vax_ip1)
# How does changing the rate of waning immunity affect the projected epidemic?
# Here we vary omega (natural-infection waning) and omega_vx (vaccine-induced waning) together in lockstep.
# Larger values -> faster loss of immunity -> larger and more sustained epidemic.
omega_list <- seq(0,0.1,by=0.025) # a vector of waning immunity rates
omega_out_all <- data.frame() # an empty dataset to save outputs

for(this_omega in omega_list){
  temp_param <- params # copy parameters
  temp_param[['omega']] <- this_omega # replace omega (natural waning)
  temp_param[['omega_vx']] <- this_omega # replace omega_vx in lockstep (vaccine-induced waning)
  this_out <- data.frame(ode(y = state.ip1, times = times, func = SIRS_Vax_ip1, parms = temp_param)) # run ode solver
  # compute total S and save omega in this cycle
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I + R,
           omega = this_omega,
           S = S_NotV + S_V)
  # Stack the result
  if(this_omega == omega_list[1]){
    omega_out_all <- this_out
  }else{
    omega_out_all <- rbind(omega_out_all, this_out)
  }
}
# Change the outcome data to long-form data
omega_out_all_t <- melt(omega_out_all%>% select(time,S,I,R,N,omega), id.vars=c("time","omega","N"))
# Plot S, I, R
ggplot(omega_out_all_t) +
  geom_line(aes(x=time,y=value/N*100,color=variable))+
  ylab("Percentage of population")+
  facet_wrap(.~omega)+
  theme_bw()
# Plot prevalence of I (I/N)
ggplot(omega_out_all) +
  geom_line(aes(x=time,y=I/N, color=as.character(omega), group=omega))+
  ylab("Proportion of population who are infected")+
  labs(color="omega")+
  theme_bw()
# Plot cumulative number of infections
ggplot(omega_out_all) +
  geom_line(aes(x=time,y=C, color=as.character(omega), group=omega))+
  ylab("Cumulative infections")+
  labs(color="omega")+
  theme_bw()
