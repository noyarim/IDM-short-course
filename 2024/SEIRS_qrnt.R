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

