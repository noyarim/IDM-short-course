#Set up and load packages
library(deSolve)
library(ggplot2)
library(dplyr)

##################
# STARTING OUT   #
##################

#We need a model that includes vaccination + waning immunity for both natural and
#vaccine-induced immunity. Starting from SIRS_Vax_ip1 in SIRS_closed_vax.R, which
#models imperfect protection with no change in infectiousness for breakthrough infections
SIRS_Vax_ip1 <- function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S_NotV + S_V + I + R

    #susceptible compartments - now stratified by vx status
    dS_NotV <- -beta*S_NotV*I/N - mu*S_NotV + omega_vx*S_V + omega*R
    dS_V <- -beta*S_V*I/N*(1-alpha) + mu*S_NotV - omega_vx*S_V
    #single infected compartment (vaccinated and unvaccinated infected are equally infectious)
    dI <- beta*S_NotV*I/N + beta*S_V*I/N*(1-alpha) - gamma*I
    #single recovered compartment (both natural and breakthrough recoveries flow here)
    dR <- gamma*I - omega*R

    #cumulative number of cases (both naive and breakthrough infections)
    dC <- beta*S_NotV*I/N + beta*S_V*I/N*(1-alpha)
    
    #cumulative number of vaccinations (use for costing later in the exercise)
    dV <- mu*S_NotV
    
    list(c(dS_NotV, dS_V, dI, dR, dC, dV))
  })
}

#Define parameters consistent with Parts 1 and 2
gamma <- 1/14                    #recovery rate (1/14 per day)
R0_bounds <- c(3, 10)
beta_bounds <- gamma*R0_bounds   #best- and worst-case beta (R0 = 3 and 10)
omega <- 1/(3*30.5)              #natural-infection immunity wanes after ~3 months
omega_vx <- 1/365                #vaccine-induced immunity wanes after 1 year

#Vaccine-specific parameters
N_pop <- 10000
muA <- 150/N_pop   #Vaccine A: ~150 people/day
muB <- 50/N_pop   #Vaccine B: ~50 people/day
alphaA <- 0.65    #Vaccine A: 65% effective
alphaB <- 0.90    #Vaccine B: 90% effective

#Common params (beta, mu, alpha get set per scenario below)
params <- c(gamma=gamma,
            omega=omega,
            omega_vx=omega_vx)

#Initial state - 1 infected, everyone else unvaccinated susceptible
state <- c(S_NotV = N_pop-1,
           S_V = 0,
           I = 1,
           R = 0,
           C = 0,
           V = 0)

T_end <- 300                     #300-day horizon
times <- seq(0, T_end, by = 1)

#No-vaccine baseline (mu=0, alpha=0 reduces SIRS_Vax_ip1 to a closed SIRS)
params['mu'] <- 0
params['alpha'] <- 0
params['beta'] <- beta_bounds[1]
out_novx_best <- data.frame(ode(y = state, times = times, func = SIRS_Vax_ip1, parms = params))
params['beta'] <- beta_bounds[2]
out_novx_worst <- data.frame(ode(y = state, times = times, func = SIRS_Vax_ip1, parms = params))

#Vaccine A: 65% effective, ~150/day
params['mu'] <- muA
params['alpha'] <- alphaA
params['beta'] <- beta_bounds[1]
out_vxA_best <- data.frame(ode(y = state, times = times, func = SIRS_Vax_ip1, parms = params))
params['beta'] <- beta_bounds[2]
out_vxA_worst <- data.frame(ode(y = state, times = times, func = SIRS_Vax_ip1, parms = params))

#Vaccine B: 90% effective, ~50/day
params['mu'] <- muB
params['alpha'] <- alphaB
params['beta'] <- beta_bounds[1]
out_vxB_best <- data.frame(ode(y = state, times = times, func = SIRS_Vax_ip1, parms = params))
params['beta'] <- beta_bounds[2]
out_vxB_worst <- data.frame(ode(y = state, times = times, func = SIRS_Vax_ip1, parms = params))

#################
#Question 1     #
#################
results_q1 <- data.frame(
  scenario = c("No vaccine", "Vaccine A", "Vaccine B"),
  best_case_R0_3 = c(round(out_novx_best[T_end+1, 'C']),
                     round(out_vxA_best[T_end+1, 'C']),
                     round(out_vxB_best[T_end+1, 'C'])),
  worst_case_R0_10 = c(round(out_novx_worst[T_end+1, 'C']),
                       round(out_vxA_worst[T_end+1, 'C']),
                       round(out_vxB_worst[T_end+1, 'C']))
)
print(results_q1)

#Plot cumulative infections over time to visualize the comparison
out_novx <- bind_rows(out_novx_best %>% mutate(scenario="No vaccine", case="Best case (R0=3)"),
                     out_novx_worst %>% mutate(scenario="No vaccine", case="Worst case (R0=10)"))
out_vxA <- bind_rows(out_vxA_best %>% mutate(scenario="Vaccine A", case="Best case (R0=3)"),
                     out_vxA_worst %>% mutate(scenario="Vaccine A", case="Worst case (R0=10)"))
out_vxB <- bind_rows(out_vxB_best %>% mutate(scenario="Vaccine B", case="Best case (R0=3)"),
                     out_vxB_worst %>% mutate(scenario="Vaccine B", case="Worst case (R0=10)"))
all_out <- bind_rows(out_novx, out_vxA, out_vxB)

ggplot(all_out, aes(x=time, y=C, color=scenario)) +
  geom_line(linewidth=1) +
  facet_wrap(~case) +
  labs(x="Day", y="Cumulative infections", color="") +
  ggtitle("Cumulative infections over 300 days by vaccine option")

#################
#Question 2     #
#################

#specify additional parameters: costs
#note: not needed in "params" vectors since these aren't going into the ODE solver
c_vxA <- 100    #$100/dose
c_vxB <- 30    #$30/dose
c_vx_initial <- 10000  #start-up costs (e.g., comms, training, etc.)
c_illness <- 50   #$50 per day spent sick

#specify additional parameters: QALYs
q_vxA <- 1
q_vxB <- 0.95
q_illness <- 0.75

#No vaccination: calculate costs and QALYs
out_novx_best <- out_novx_best %>% 
  mutate(N = S_NotV + S_V + I + R,
         inc_vx = 0,
         cum_costs = cumsum(c_illness*I),
         cum_qalys = cumsum(1*N - (1-q_illness)*I))

#Vaccine A: calculate costs and QALYs
out_vxA_best <- out_vxA_best %>% 
  mutate(N = S_NotV + S_V + I + R,
         inc_V = V - lag(V, 1, default=0),
         cum_costs = cumsum(c_illness*I + c_vxA*inc_V) + c_vx_initial,
         cum_qalys = cumsum(1*N - ((1-q_illness)*I + (1-q_vxA)*inc_V*3)))

#Vaccine B: calculate costs and QALYs
out_vxB_best <- out_vxB_best %>% 
  mutate(N = S_NotV + S_V + I + R,
         inc_V = V - lag(V, 1, default=0),
         cum_costs = cumsum(c_illness*I + c_vxB*inc_V) + c_vx_initial,
         cum_qalys = cumsum(1*N - ((1-q_illness)*I + (1-q_vxB)*inc_V*3)))

#Report total cost and QALYs under each policy
results_q2 <- data.frame(
  scenario = c("No vaccine", "Vaccine A", "Vaccine B"),
  costs_pc = round(c(out_novx_best[T_end+1, 'cum_costs'],
                     out_vxA_best[T_end+1, 'cum_costs'],
                     out_vxB_best[T_end+1, 'cum_costs'])/N_pop),
  qalys_pc = round(c(out_novx_best[T_end+1, 'cum_qalys'],
                     out_vxA_best[T_end+1, 'cum_qalys'],
                     out_vxB_best[T_end+1, 'cum_qalys'])/(365*N_pop), 3)
)
#calculated per capita for easier interpretation
#note we divide QALYs by 365 since QALYs are in years and our timestep is days
print(results_q2)

#################
#Question 3     #
#################

#arrange results by cost and calculate inc costs and QALYs
results_q3 <- results_q2 %>% arrange(costs_pc)
lowest_cost_strategy <- results_q3[1, "scenario"]
results_q3 <- results_q3 %>% 
  mutate(inc_costs_pc=costs_pc-lag(costs_pc, 1, default=as.numeric(NA)),
         inc_qalys_pc=qalys_pc-lag(qalys_pc, 1 ,default=as.numeric(NA)))
print(results_q3)
#we clearly don't want to do "no vaccine" - it costs the most for the least QALYs
results_q3 <- results_q3 %>% filter(scenario!="No vaccine")

#calculate ICERs
results_q3 <- results_q3 %>%
  mutate(ICER=if_else(scenario==lowest_cost_strategy, as.numeric(NA),
                      inc_costs_pc/inc_qalys_pc))
print(results_q3)

# The ICER of Vaccine A vs. Vaccine B is $7000
# We prefer Vaccine B at a WTP of $5000
# But we prefer Vaccine A at a WTP of $10,000

#################
#Question 4     #
#################

#we need to update our model function to allow for lower infectiousness of breakthrough infections
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
    
    #cumulative numbers vaccinated
    dV <- mu*S_NotV
    
    # return the rates of change as a list
    list(c(dS_NotV,dS_V,dI_NotV,dI_V, dR, dC, dV))
  })
}

#update population starting distribution
state <- c(S_NotV = N_pop-1,
           S_V = 0,
           I_NotV = 1,
           I_V = 0,
           R = 0,
           C = 0,
           V = 0)

#add mbeta to params
params <- c(params, 
            "mbeta" = 0.5)

#rerun vaccine A best case scenario
params['mu'] <- muA
params['alpha'] <- alphaA
params['beta.1'] <- beta_bounds[1]
out_vxA_best_mbeta <- data.frame(ode(y = state, 
                                     times = times, 
                                     func = SIRS_Vax_ip2, 
                                     parms = params))

#plot the difference
out_vxA_best_all <- bind_rows(out_vxA_best %>% mutate(label="No Reduced Infectiousness"),
                              out_vxA_best_mbeta %>% mutate(label="50% Reduced Infectiousness"))

ggplot(out_vxA_best_all, aes(x=time, y=C, color=label)) +
  geom_line(linewidth=1) +
  labs(x="Day", y="Cumulative infections", color="") +
  ggtitle("Cumulative infections over 300 days by scenario")

#redo the cost-effectiveness
out_vxA_best_mbeta <- out_vxA_best_mbeta %>% 
  mutate(N = S_NotV + S_V + I_V + I_NotV + R,
         inc_V = V - lag(V, 1, default=0),
         cum_costs = cumsum(c_illness*(I_V + I_NotV) + c_vxA*inc_V) + c_vx_initial,
         cum_qalys = cumsum(1*N - ((1-q_illness)*(I_V + I_NotV) + (1-q_vxA)*inc_V*3)))
results_q4 <- data.frame(
  scenario = c("No vaccine", "Vaccine A", "Vaccine B"),
  costs_pc = round(c(out_novx_best[T_end+1, 'cum_costs'],
                     out_vxA_best_mbeta[T_end+1, 'cum_costs'],
                     out_vxB_best[T_end+1, 'cum_costs'])/N_pop),
  qalys_pc = round(c(out_novx_best[T_end+1, 'cum_qalys'],
                     out_vxA_best_mbeta[T_end+1, 'cum_qalys'],
                     out_vxB_best[T_end+1, 'cum_qalys'])/(365*N_pop), 3)
)
print(results_q4)

#now vaccine A is cost-saving and preferred across all WTP thresholds
#no need to calculate ICERs