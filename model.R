model <- function(t, state, parms){
  with(as.list(c(state, parms)),{
    rho <- r.2*c*x.1 + r.4*c*x.2
    
    dC = x.1*x.2*r.1 -r.3*c - rho*c
    dX1 = x.1*c*r.2 - x.1*x.2*r.1 + r.3*c - rho*x.1
    dX2 = x.2*c*r.4 - x.1*x.2*r.1 + r.3*c - rho*x.2
    
    
    return(list(c(dX1, dX2, dC)))
  })
}

p <- c(r.1 = 0.5, r.2 = 0.50002, r.3 = 0.5, r.4 = 0.5)
s <- c(x.1 = 0.2, x.2 = 0.2, c=.4)

run(tmax = 1000)