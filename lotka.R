model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    
    dR <- r*R*(1 - R/K) - a*R*N
    dN <- c*a*R*N - delta*N
    
    return(list(c(dR, dN)))  
  }) 
}  

p <- c(r=1,K=1,a=1,c=1,delta=0.5) # p is a named vector of parameters
s <- c(R=1,N=0.01)                # s is the state
