#Set up and load packages
library(deSolve)
library(ggplot2)
library(dplyr)


# SEIRS model with quarantine
SEIRS_qrtn<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S + E + I + R + Q
    # Differential equations
    dS <- -beta*S*I/N + omega*R
    dE <- beta*S*I/N - sigma*E 
    dI <- sigma*E  - gamma*I - q*I
    dQ <- q*I  - gamma*Q
    dR <- gamma*I + gamma*Q - omega*R
    # Cumulative number of infection
    dC <- beta*S*I/N
    # return the rates of change as a list
    list(c(dS, dE, dI, dQ, dR, dC))
  })
}


#################
#Question 1     #
#################
# Define model parameters (gamma, beta, omega, sigma same as in part 2)
gamma <- 1/14
R0_bounds <- c(3, 10)
beta_bounds <- gamma*R0_bounds
omega <- 1/(3*30.5) #1 over time to immunological waning - in days, since gamma is in days
sigma <- 1/5 
# quarantine rates of two policies
q_bounds <- c(0.2, 0.5) 

# Define parameters and starting compartment sizes
params <- c(beta = beta_bounds[1], #effective contact rate (aka transmission rate)
            gamma = gamma, #recovery rate (1/duration infection)
            omega = omega,#0.01, # waning immunity
            sigma = sigma, # latent period from E to I
            q = q_bounds[1] # quarantine rate after implementing quarantine
)

#Also need to update "state" vector - to add Q
state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           Q = 0,
           R = 0,
           C = 0
)

T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 

# Best case with low quarantine coverage
params['beta'] <- beta_bounds[1] 
params['q'] <- q_bounds[1]
out_best4 <- data.frame(ode(y = state, times=times, func=SEIRS_qrtn, parms = params))
# Best case with high quarantine coverage
params['q'] <- q_bounds[2]
out_best5 <- data.frame(ode(y = state, times=times, func=SEIRS_qrtn, parms = params))
# Worst case with low quarantine coverage
params['beta'] <- beta_bounds[2] 
params['q'] <- q_bounds[1]
out_worst4 <- data.frame(ode(y = state, times=times, func=SEIRS_qrtn, parms = params))
# Worst case with high quarantine coverage
params['q'] <- q_bounds[2]
out_worst5 <- data.frame(ode(y = state, times=times, func=SEIRS_qrtn, parms = params))

# Cumulative infections at 300 days
print(paste("Best case with low quarantine coverage:", round(out_best4[T_end+1, 'C']), "cumulative cases"))
print(paste("Best case with high quarantine coverage:", round(out_best5[T_end+1, 'C']), "cumulative cases"))
print(paste("Worst case with low quarantine coverage:", round(out_worst4[T_end+1, 'C']), "cumulative cases"))
print(paste("Worst case with high quarantine coverage:", round(out_worst5[T_end+1, 'C']), "cumulative cases"))


# We can also plot trajectories over time, to compare 
ggplot() + 
  geom_line(data=out_best4, aes(x=time, y=I, color="Best Case, low quarantine coverage")) +
  geom_line(data=out_worst4, aes(x=time, y=I, color="Worst Case, low quarantine coverage")) +
  geom_line(data=out_best5, aes(x=time, y=I, color="Best Case, high quarantine coverage")) +
  geom_line(data=out_worst5, aes(x=time, y=I, color="Worst Case, high quarantine coverage")) +
  labs(x="Day", y="Number of people infected", color="") +
  ggtitle("Cases over time")

#################
#Question 2     #
#################
# An isolation policy with high coverage can be less preferred due to several factors. 
# Firstly, the costs associated with implementing high coverage isolation are significantly higher than those for low coverage. 
# Achieving high isolation coverage can strain the current healthcare system capacity and may necessitate renovations or expansions, 
# which are expensive and time-consuming. Secondly, compliance with quarantine policies can vary, as not everyone may agree to follow them. 
# This non-compliance can delay the implementation of a high coverage isolation policy, thus hampering the early control of an outbreak.
# The code below shows the outcome with high coverage quarantine when there is a delay in implementation

# Incorporate different timing of isolation policies (i.e. time-dependent q)
SEIRS_qrtn2<-function(t, state, params) {
  with(as.list(c(state, params)),{
    N = S + E + I + R + Q
    q_t = ifelse(t<t_start_q, 0, q) # time-varying q: before implementing quarantine, q(t) = 0
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + omega*R
    dE <- beta*S*I/N - sigma*E 
    dI <- sigma*E - gamma*I - q_t*I
    dQ <- q_t*I - gamma*Q
    dR <- gamma*I + gamma*Q - omega*R
    # Cumulative number of infection
    dC <- beta*S*I/N
    # return the rates of change as a list
    list(c(dS, dE, dI, dQ, dR, dC))
  })
}

# Parameters
params <- c(beta = beta_bounds[1], #effective contact rate (aka transmission rate)
            gamma = gamma, #recovery rate (1/duration infection)
            omega = omega, # waning immunity
            sigma = sigma, # latent period from E to I
            q = q_bounds[2], # rate for high coverage
            t_start_q = 50 # time to implement the high coverage policy
)

# Initial state
state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           Q = 0,
           R = 0,
           C = 0
)

# Worst case with high quarantine coverage but with delay
params['beta'] <- beta_bounds[2]
params['q'] <- q_bounds[2]
out_worst5_delay <- data.frame(ode(y = state, times=times, func=SEIRS_qrtn2, parms = params))

# Plot I to compare low coverage, high coverage without delay, and high coverage with delay in the best case scenario
ggplot() + 
  geom_line(data=out_worst4, aes(x=time, y=I, color="Worst Case, low quarantine coverage")) +
  geom_line(data=out_worst5, aes(x=time, y=I, color="Worst Case, high quarantine coverage, no delay")) +
  geom_line(data=out_worst5_delay, aes(x=time, y=I, color="Worst Case, high quarantine coverage, with delay")) +
  labs(x="Day", y="Number of people infected", color="") +
  ggtitle("Cases over time")


#######################
#### PART 4 - EXT ####
######################

#################
#Question 1     #
#################
# Define parameters for CEA
CEA_params <- c(
  # total cost
  c_lowcov <- 5000,
  c_highcov <- 500000,
  # disutility
  duS <- 0,
  duE <- 0,
  duI <- 0.3,
  duR <- 0,
  duQ <- 0.3,
  # daily cost
  dcS <- 0,
  dcE <- 0,
  dcI <- 50,
  dcR <- 0,
  dcQ <- 100
)

# A function to calculate the total QALY and Cost given the output table and strategy
cal_CEA <- function(out, strategy){
  with(as.list(CEA_params), {
    out_cea <- as.data.frame(out) %>%
      mutate(
        uS = S * (1-duS), #e.g. S: the number of susceptibles, (1-duS): QALY per cycle per individual
        uE = E * (1-duE),
        uI = I * (1-duI),
        uR = R * (1-duR),
        uQ = Q * (1-duQ),
        uTotal_percycle = uS + uE + uI + uR + uQ, # total utility per day (column sum)
        cS = S * dcS,
        cE = E * dcE,
        cI = I * dcI,
        cR = R * dcR,
        cQ = Q * dcQ,
        cTotal_percycle = cS + cE + cI + cR + cQ # total cost per day (column sum)
      )
    # Total QALYs over 300 days (row sum)
    tot_QALY = sum(out_cea$uTotal_percycle)
    # Total cost over 300 days (row sum and add one-time cost)
    if (strategy == "low"){
      tot_Cost = sum(out_cea$cTotal_percycle) + c_lowcov
    }else{
      tot_Cost = sum(out_cea$cTotal_percycle) + c_highcov
    }
    
    return(c(tot_QALY = tot_QALY, tot_Cost = tot_Cost))
  })
}

# Best case, low coverage
cea_best4 <- cal_CEA(out_best4, 'low')
# Best case, high coverage
cea_best5 <- cal_CEA(out_best5, 'high')
# Worst case, low coverage
cea_worst4 <- cal_CEA(out_worst4,'low')
# Worst case, high coverage
cea_worst5 <- cal_CEA(out_worst5, 'high')

data.frame( Scenario = c("Best", "Best", "Worst", "Worst"),
            Policy = c("low","high","low","high"),
            QALY = c(cea_best4[1], cea_best5[1], cea_worst4[1], cea_worst5[1]),
            Cost = c(cea_best4[2], cea_best5[2], cea_worst4[2], cea_worst5[2])
)


#################
#Question 2     #
#################
# Calculate ICER
icer_best <- (cea_best5['tot_Cost'] - cea_best4['tot_Cost'])/(cea_best5['tot_QALY']-cea_best4['tot_QALY'])
icer_worst <- (cea_worst5['tot_Cost'] - cea_worst4['tot_Cost'])/(cea_worst5['tot_QALY']-cea_worst4['tot_QALY'])

print(paste0("Best scenario: ICER of high coverage policy compared to low coverage policy is:", round(icer_best,2)))
print(paste0("Worst scenario: ICER of high coverage policy compared to low coverage policy is:", round(icer_worst,2)))
