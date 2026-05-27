###############################################################################################

#SIR model with vaccination (closed population, no demography)

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

## A. Vaccination provides perfect protection from infection: S -> R ##
# Direct flow from S to R indicates those who are vaccinated and protected from infection.

#1. Define model function
SIR_Vax_pp<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S + I + R

    #SIR equations from lecture
    dS <- -beta*S*I/N - mu*S
    dI <- beta*S*I/N - gamma*I
    dR <- gamma*I + mu*S
    # cumulative incidence
    dC <- beta*S*I/N
    # return the rates of change as a list
    list(c(dS, dI, dR, dC))
  })
}

#2. Define parameters and starting compartment sizes
# Mu is vaccination rate.
params <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                mu = 0.01 # vaccination rate
)

# Initial state distribution
state.pp <- c(S = 99999, #population of 100,000, 1 person starts of infected
              I = 1,
              R = 0,
              C = 0)

#3. Run the SIR model with perfect protection from vaccination
T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step

# Run the base-case
output.vax.pp <- ode(y = state.pp, times = times, func = SIR_Vax_pp, parms = params)

# Plot the base-case result
plot(output.vax.pp)

# Compare the size of I compartment when mu = 0 (no vax) and mu = 0.01 (vax)
# Copy parameters
temp_param <- params
# Set vaccination rate = 0
temp_param['mu'] = 0
# Run ode solver and save the outcome in temporary output variable
temp_output <-  ode(y = state.pp, times = times, func = SIR_Vax_pp, parms = temp_param)
# Create a dataset to compare no_vax and vax output
I_dt <- data.frame(time = temp_output[,'time'], no_vax = temp_output[,'I'], vax = output.vax.pp[,'I'])
# Plot the dataset
ggplot(I_dt)+
  geom_line(aes(time,vax,color='vax'))+
  geom_line(aes(time,no_vax,color='No vax'))+
  theme_bw()+
  scale_color_manual(values=c("blue","tomato"))


## B. Vaccination provides imperfect protection from infection ##
# Vaccination usually doesn't provide perfect protection. Here, we stratify by vaccine status.
# Non-vaccinated compartments move to vaccinated compartments at rate mu.
# Total I who can infect others = sum of I_NotV and I_V.
# (1-alpha) indicates reduced susceptibility for vaccinated susceptibles; higher alpha -> lower susceptibility.

# B-1. No change in infectiousness
#1. Define model function
SIR_Vax_ip1<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V

    #compartments without vaccination
    dS_NotV <- -beta*S_NotV*(I_NotV+I_V)/N - mu*S_NotV
    dI_NotV <- beta*S_NotV*(I_NotV+I_V)/N - gamma*I_NotV
    dR_NotV <- gamma*I_NotV
    #compartments with vaccination
    dS_V <- -beta*S_V*(I_NotV+I_V)/N*(1-alpha) + mu*S_NotV
    dI_V <- beta*S_V*(I_NotV+I_V)/N*(1-alpha) - gamma*I_V
    dR_V <- gamma*I_V

    #cumulative number of cases
    dC <- beta*S_NotV*(I_NotV+I_V)/N + beta*S_V*(I_NotV+I_V)/N*(1-alpha)

    # return the rates of change as a list
    list(c(dS_NotV,dS_V,dI_NotV,dI_V, dR_NotV,dR_V,dC))
  })
}

# B-2. Change in risk of infecting others with vaccination
# Apply separate beta values (beta.1 and beta.2) to I_NotV and I_V.
# Beta.2 is calculated as a fraction (mbeta) of beta.1.
SIR_Vax_ip2<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V
    beta.2 = mbeta*beta.1 # transmission rate given vaccination
    #compartments without vaccination
    dS_NotV <- -beta.1*S_NotV*I_NotV/N - beta.2*S_NotV*I_V/N - mu*S_NotV
    dI_NotV <- beta.1*S_NotV*I_NotV/N + beta.2*S_NotV*I_V/N - gamma*I_NotV
    dR_NotV <- gamma*I_NotV
    #compartments with vaccination
    dS_V <- -beta.1*S_V*I_NotV/N*(1-alpha) - beta.2*S_V*I_V/N*(1-alpha) + mu*S_NotV
    dI_V <- beta.1*S_V*I_NotV/N*(1-alpha) + beta.2*S_V*I_V/N*(1-alpha) - gamma*I_V
    dR_V <- gamma*I_V

    #cumulative number of cases
    dC <- beta.1*S_NotV*I_NotV/N + beta.1*S_V*I_NotV/N*(1-alpha) + beta.2*S_NotV*I_V/N + beta.2*S_V*I_V/N*(1-alpha)

    # return the rates of change as a list
    list(c(dS_NotV,dS_V,dI_NotV,dI_V, dR_NotV,dR_V, dC))
  })
}

#2. Define parameters and starting compartment sizes
# Mu is vaccination rate.
# mbeta is 1 if no change in infectiousness with vaccination, or less than 1 if reduced infectiousness.
# alpha is vaccine effectiveness.
params <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                beta.1 = 0.5, # transmission rate given no-vaccination
                mbeta = 0.5, # 1 if no change in infectiousness with vaccination,
                           # less than 1 if reduced infectiousness with vaccination
                gamma = 0.3, #recovery rate (1/duration infection)
                mu = 0.01, # vaccination rate
                alpha = 0.3 # vaccine effectiveness
)

# Initial state
state.ip <- c(S_NotV = 99999,
              S_V = 0,
              I_NotV = 1,
              I_V = 0,
              R_NotV = 0,
              R_V = 0,
              C = 0)

#3. Run the SIR model with imperfect protection from vaccination
# Because each compartment is stratified, we compute the total S, I, R compartments
T_end <- 500 #run model for 500 time steps (e.g. months)
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps (e.g. months), and computes output at each time step

# Run ODE solver
output.vax.ip1 <- ode(y = state.ip, times = times, func = SIR_Vax_ip1, parms = params)
output.vax.ip2 <- ode(y = state.ip, times = times, func = SIR_Vax_ip2, parms = params)
# Calculate S, I, and R
output.vax.ip1 <- as.data.frame(output.vax.ip1) %>%
  mutate(S = S_V + S_NotV,
         I = I_V + I_NotV,
         R = R_V + R_NotV)
output.vax.ip2 <- as.data.frame(output.vax.ip2) %>%
  mutate(S = S_V + S_NotV,
         I = I_V + I_NotV,
         R = R_V + R_NotV)

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
# Perform sensitivity analysis on (1) vaccination rate, (2) vaccine effectiveness, and
# (3) the ratio in effective contact rate between vaccinated and non-vaccinated individuals.
# (1) and (2) use SIR_Vax_ip1; (3) uses SIR_Vax_ip2.

# 1. Vaccination rate (on SIR_Vax_ip1)
# Discuss how changing vaccination rate changes the projected epidemic.
mu_list <- seq(0,0.1,length.out=3) # a vector of vaccination rate
mu_out_all <- data.frame() # an empty dataset to save outputs

for(this_mu in mu_list){
  temp_param <- params # copy parameters
  temp_param[['mu']] <- this_mu # replace mu
  this_out <- data.frame(ode(y = state.ip, times = times, func = SIR_Vax_ip1, parms = temp_param)) # run ode solver
  # compute S, I, R and save mu in this cycle
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V,
           mu = this_mu,
           S = S_NotV + S_V,#total susceptible
           I = I_NotV + I_V, # total infected
           R = R_NotV + R_V) # total recovered
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
# Plot I stratified by vaccination status
ggplot(mu_out_all)+
  geom_point(aes(x=time,y=I_NotV,color='Not vaccinated'))+
  geom_point(aes(x=time,y=I_V,color='Vaccinated'))+
  ylab("Infections")+
  facet_wrap(.~mu)+
  scale_color_manual(values=c("green","tomato"))+
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
  this_out <- data.frame(ode(y = state.ip, times = times, func = SIR_Vax_ip1, parms = temp_param)) # run ode solver
  # compute S, I, R and save mu in this cycle
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V,
           alpha = this_alpha,
           S = S_NotV + S_V,
           I = I_NotV + I_V,
           R = R_NotV + R_V)
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
# Plot S, I, R stratified by vaccination status
ggplot(alpha_out_all)+
  geom_point(aes(x=time,y=I_NotV,color='Not vaccinated'))+
  geom_point(aes(x=time,y=I_V,color='Vaccinated'))+
  ylab("Infections")+
  facet_wrap(.~alpha)+
  scale_color_manual(values=c("green","tomato"))+
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
# Breakthrough infection
ggplot(alpha_out_all)+
  geom_point(aes(x=time,y=I_V))+
  facet_wrap(.~alpha)+
  ylab("Number of breakthrough infections")+
  theme_bw()


# 3. Ratio of effective contact rate among infected individuals with and without vaccination
# How does the projected epidemic change as we decrease the ratio of effective contact rate?
mbeta_list <- seq(1,0.7,length.out=3) # a vector of mbeta
mbeta_out_all <- data.frame() # an empty dataset to save outputs

for(this_mbeta in mbeta_list){
  temp_param <- params # copy parameters
  temp_param[['mbeta']] <- this_mbeta # replace mbeta
  this_out <- data.frame(ode(y = state.ip, times = times, func = SIR_Vax_ip2, parms = temp_param)) # run ode solver
  # compute S,I,R and save mbeta in this cycle
  this_out <- this_out %>%
    mutate(N = S_NotV + S_V + I_NotV + I_V + R_NotV + R_V,
           mbeta = this_mbeta,
           S = S_NotV + S_V,
           I = I_NotV + I_V,
           R = R_NotV + R_V)
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
