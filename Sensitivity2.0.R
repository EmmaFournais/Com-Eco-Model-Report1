
#Author: Emma Fournais Næsheim
#Computational Marine Ecological Modelling - report 1
#04-03-2024

library(deSolve)


par(mfrow=c(2,2))
u <- 0.5 #

d <- 100
n <- 100 # 
DeltaZ <- d/n # 
kw <- 0.2
kp <- 0.05 #

gmax <- 1 #/day
I0 <- 350 #



kn <- 0.3
m <- 0.24 #/day
alpha <- 0.1

D <- 4.32 # m^2/day
D1 <- 20
D2 <- 60
D3 <- 80
z_values <- seq( 0.5*DeltaZ, d, by = DeltaZ)
grazing <- 0.1 #
remin <- 1.5 # 
remin1 <- 2
remin2 <- 2.5
remin3 <- 3

# Initial conditions 
Nb <- 5
P <- rep(5,n)
N <- rep(0.2,n)
De <- rep(0,n)


#--------------------------------------------------------
Light <- function(kw,kp,DeltaZ,P,I0){
  
  #Integral calculations - cumulative sums of the light in depth
  integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
  
  #calculation the light values 
  I = I0*exp(-integral)
  
  return(I)
}

#------------------------------------------------------------


G.n <- N/(kn+N)
G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
g <- gmax * pmin(G.n, G.l)- m


#------------------------------------------------
#Doing the diffusion and advection 
#------------------------------------------------

derivative <- function(times, y, params) {
  P <- y[1:n]
  N <- y[(n+1):(2*n)]
  De <- y[(2*n+1):(3*n)]
  
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  J_n <- numeric(n + 1)
  J_D <- numeric(n + 1)
  
  for(i in 2:n){
    J_a[i] <-  u*P[i-1]
    J_d[i] <- -D*(P[i]-P[i-1])/DeltaZ
    J_n[i] <- -D*(N[i]-N[i-1])/DeltaZ 
    J_D[i] <- -D*(De[i]-De[i-1])/DeltaZ 
  }
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_n[1] <- 0 
  J_D[1] <- 0 
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  J_D[n + 1] <- u*De[n]
  J_n[n + 1] <- -D * ((Nb-N[n])/DeltaZ)
  
  J = J_a + J_d
  
  jP= (-(J[2:(n+1)] - J[1:n])) / DeltaZ
  jN= -(J_n[2:(n+1)] - J_n[1:n]) / DeltaZ
  jD= -(J_D[2:(n+1)] - J_D[1:n]) / DeltaZ
  
  
  Light <- function(kw,kp,DeltaZ,P,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  
  dP_dt <- numeric(n)
  dN_dt <- numeric(n)
  dD_dt <- numeric(n)
  
  for (i in 1:n){
    dP_dt[i] <- (g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + jP[i]
    dN_dt[i] <- (-g[i]*P[i])+ (remin*De[i]) + jN[i]
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin*De[i]) - (u*De[i]) + jD[i]
  }
  
  return(list(c(dP_dt, dN_dt,dD_dt)))
}


params <-
  list(
    n = n,
    DeltaZ = DeltaZ,
    u = u,
    D = D,
    m = m,
    g = g, 
    gmax = gmax, 
    I = I, 
    Nb = Nb, 
    yield = yield, 
    Hi = Hi,
    Hn = Hn, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 100, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out <- result[, 2:(n+1)]
N_out <- result[, (n+2):(2*n+1)]
D_out <- result[, (2*n+2):(3*n+1)]
Time <- result [,1]
z_values <- z_values



#-----------D1--------------------------------------------------

derivative <- function(times, y, params) {
  P <- y[1:n]
  N <- y[(n+1):(2*n)]
  De <- y[(2*n+1):(3*n)]
  
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  J_n <- numeric(n + 1)
  J_D <- numeric(n + 1)
  
  for(i in 2:n){
    J_a[i] <-  u*P[i-1]
    J_d[i] <- -D1*(P[i]-P[i-1])/DeltaZ
    J_n[i] <- -D1*(N[i]-N[i-1])/DeltaZ 
    J_D[i] <- -D1*(De[i]-De[i-1])/DeltaZ 
  }
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_n[1] <- 0 
  J_D[1] <- 0 
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  J_D[n + 1] <- u*De[n]
  J_n[n + 1] <- -D1 * ((Nb-N[n])/DeltaZ)
  
  J = J_a + J_d
  
  jP= (-(J[2:(n+1)] - J[1:n])) / DeltaZ
  jN= -(J_n[2:(n+1)] - J_n[1:n]) / DeltaZ
  jD= -(J_D[2:(n+1)] - J_D[1:n]) / DeltaZ
  
  
  Light <- function(kw,kp,DeltaZ,P,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  
  dP_dt <- numeric(n)
  dN_dt <- numeric(n)
  dD_dt <- numeric(n)
  
  for (i in 1:n){
    dP_dt[i] <- (g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + jP[i]
    dN_dt[i] <- (-g[i]*P[i])+ (remin*De[i]) + jN[i]
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin*De[i]) - (u*De[i]) + jD[i]
  }
  
  return(list(c(dP_dt, dN_dt,dD_dt)))
}


params <-
  list(
    n = n,
    DeltaZ = DeltaZ,
    u = u,
    D = D,
    m = m,
    g = g, 
    gmax = gmax, 
    I = I, 
    Nb = Nb, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 100, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out1 <- result[, 2:(n+1)]
N_out1 <- result[, (n+2):(2*n+1)]
D_out1 <- result[, (2*n+2):(3*n+1)]
Time <- result [,1]
z_values <- z_values


#---------D2------------------------------------------------------
derivative <- function(times, y, params) {
  P <- y[1:n]
  N <- y[(n+1):(2*n)]
  De <- y[(2*n+1):(3*n)]
  
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  J_n <- numeric(n + 1)
  J_D <- numeric(n + 1)
  
  for(i in 2:n){
    J_a[i] <-  u*P[i-1]
    J_d[i] <- -D2*(P[i]-P[i-1])/DeltaZ
    J_n[i] <- -D2*(N[i]-N[i-1])/DeltaZ 
    J_D[i] <- -D2*(De[i]-De[i-1])/DeltaZ 
  }
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_n[1] <- 0 
  J_D[1] <- 0 
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  J_D[n + 1] <- u*De[n]
  J_n[n + 1] <- -D2 * ((Nb-N[n])/DeltaZ)
  
  J = J_a + J_d
  
  jP= (-(J[2:(n+1)] - J[1:n])) / DeltaZ
  jN= -(J_n[2:(n+1)] - J_n[1:n]) / DeltaZ
  jD= -(J_D[2:(n+1)] - J_D[1:n]) / DeltaZ
  
  
  Light <- function(kw,kp,DeltaZ,P,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  
  dP_dt <- numeric(n)
  dN_dt <- numeric(n)
  dD_dt <- numeric(n)
  
  for (i in 1:n){
    dP_dt[i] <- (g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + jP[i]
    dN_dt[i] <- (-g[i]*P[i])+ (remin*De[i]) + jN[i]
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin*De[i]) - (u*De[i]) + jD[i]
  }
  
  return(list(c(dP_dt, dN_dt,dD_dt)))
}


params <-
  list(
    n = n,
    DeltaZ = DeltaZ,
    u = u,
    D = D,
    m = m,
    g = g, 
    gmax = gmax, 
    I = I, 
    Nb = Nb, 
    yield = yield, 
    Hi = Hi,
    Hn = Hn, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 100, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out2 <- result[, 2:(n+1)]
N_out2 <- result[, (n+2):(2*n+1)]
D_out2 <- result[, (2*n+2):(3*n+1)]
Time <- result [,1]
z_values <- z_values


#----------D3-------------------------------------------------------

derivative <- function(times, y, params) {
  P <- y[1:n]
  N <- y[(n+1):(2*n)]
  De <- y[(2*n+1):(3*n)]
  
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  J_n <- numeric(n + 1)
  J_D <- numeric(n + 1)
  
  for(i in 2:n){
    J_a[i] <-  u*P[i-1]
    J_d[i] <- -D3*(P[i]-P[i-1])/DeltaZ
    J_n[i] <- -D3*(N[i]-N[i-1])/DeltaZ 
    J_D[i] <- -D3*(De[i]-De[i-1])/DeltaZ 
  }
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_n[1] <- 0 
  J_D[1] <- 0 
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  J_D[n + 1] <- u*De[n]
  J_n[n + 1] <- -D3 * ((Nb-N[n])/DeltaZ)
  
  J = J_a + J_d
  
  jP= (-(J[2:(n+1)] - J[1:n])) / DeltaZ
  jN= -(J_n[2:(n+1)] - J_n[1:n]) / DeltaZ
  jD= -(J_D[2:(n+1)] - J_D[1:n]) / DeltaZ
  
  
  Light <- function(kw,kp,DeltaZ,P,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  
  dP_dt <- numeric(n)
  dN_dt <- numeric(n)
  dD_dt <- numeric(n)
  
  for (i in 1:n){
    dP_dt[i] <- (g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + jP[i]
    dN_dt[i] <- (-g[i]*P[i])+ (remin*De[i]) + jN[i]
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin*De[i]) - (u*De[i]) + jD[i]
  }
  
  return(list(c(dP_dt, dN_dt,dD_dt)))
}


params <-
  list(
    n = n,
    DeltaZ = DeltaZ,
    u = u,
    D = D,
    m = m,
    g = g, 
    gmax = gmax, 
    I = I, 
    Nb = Nb, 
    yield = yield, 
    Hi = Hi,
    Hn = Hn, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 100, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out3 <- result[, 2:(n+1)]
N_out3 <- result[, (n+2):(2*n+1)]
D_out3 <- result[, (2*n+2):(3*n+1)]
Time <- result [,1]
z_values <- z_values



plot(P_out[ dim(P_out)[2],],-z_values, type = "l", col = "darkgreen", xlim = c(0,1.5),lwd=3, main = "Sensitivity test for D in Phytoplankton", ylab = "Depth[m]", xlab = "[mmol N/m^3]", cex.lab=1.25)
lines(P_out1[ dim(P_out1)[2],], -z_values, type = "l", col = "seagreen", lwd=3)
lines( P_out2[ dim(P_out2)[2],],-z_values, type = "l",col = "green", lwd=3)
lines(P_out3[ dim(P_out3)[2],],-z_values, type = "l",col = "lightgreen",  lwd=3)
legend("topright",legend = c("D = 4.32 ", "D = 20", "D = 60", "D = 80"),lty = c(1, 1, 1, 1),col = c("darkgreen", "seagreen", "green", "lightgreen"),lwd = 2)

plot( N_out[ dim(N_out)[2],],-z_values, type = "l", col = "brown", lwd=3, main = "Sensitivity test for D in Nutrients", ylab = "Depth[m]", xlab = "[mmol N/m^3]", cex.lab=1.25)  
lines(N_out1[ dim(N_out1)[2],], -z_values,  type = "l", col = "darkorange3", lwd=3)
lines(N_out2[ dim(N_out2)[2],], -z_values,  type = "l", col = "orange", lwd=3)
lines(N_out3[ dim(N_out3)[2],], -z_values,type = "l", col = "yellow4", lwd=3)
legend("topright",legend = c("D = 4.32 ", "D = 20", "D = 60", "D = 80"),lty = c(1, 1, 1, 1),col = c("brown", "darkorange", "orange", "yellow4"),lwd = 2)


plot(D_out[ dim(D_out)[2],],-z_values,  type = "l", col = "black", lwd=3, xlim=c(0,0.2), main = "Sensitivity test for D in Detritus", ylab = "Depth[m]", xlab = "[mmol N/m^3]", cex.lab=1.25)
lines( D_out1[ dim(D_out1)[2],], -z_values, type = "l", col = "dimgray", lwd=3)
lines(D_out2[ dim(D_out2)[2],], -z_values,  type = "l", col = "darkgrey", lwd=3)
lines( D_out3[ dim(D_out3)[2],], -z_values, type = "l", col = "azure3", lwd=3)
legend("bottomright",legend = c("D = 4.32 ", "D = 20", "D = 60", "D = 80"),lty = c(1, 1, 1, 1),col = c("black", "dimgray", "darkgrey", "azure3"),lwd = 2)



#-----------Changing--remin---------------------------

derivative <- function(times, y, params) {
  P <- y[1:n]
  N <- y[(n+1):(2*n)]
  De <- y[(2*n+1):(3*n)]
  
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  J_n <- numeric(n + 1)
  J_D <- numeric(n + 1)
  
  for(i in 2:n){
    J_a[i] <-  u*P[i-1]
    J_d[i] <- -D*(P[i]-P[i-1])/DeltaZ
    J_n[i] <- -D*(N[i]-N[i-1])/DeltaZ 
    J_D[i] <- -D*(De[i]-De[i-1])/DeltaZ 
  }
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_n[1] <- 0 
  J_D[1] <- 0 
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  J_D[n + 1] <- u*De[n]
  J_n[n + 1] <- -D * ((Nb-N[n])/DeltaZ)
  
  J = J_a + J_d
  
  jP= (-(J[2:(n+1)] - J[1:n])) / DeltaZ
  jN= -(J_n[2:(n+1)] - J_n[1:n]) / DeltaZ
  jD= -(J_D[2:(n+1)] - J_D[1:n]) / DeltaZ
  
  
  Light <- function(kw,kp,DeltaZ,P,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  
  dP_dt <- numeric(n)
  dN_dt <- numeric(n)
  dD_dt <- numeric(n)
  
  for (i in 1:n){
    dP_dt[i] <- (g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + jP[i]
    dN_dt[i] <- (-g[i]*P[i])+ (remin*De[i]) + jN[i]
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin*De[i]) - (u*De[i]) + jD[i]
  }
  
  return(list(c(dP_dt, dN_dt,dD_dt)))
}


params <-
  list(
    n = n,
    DeltaZ = DeltaZ,
    u = u,
    D = D,
    m = m,
    g = g, 
    gmax = gmax, 
    I = I, 
    Nb = Nb, 
    yield = yield, 
    Hi = Hi,
    Hn = Hn, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 100, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out <- result[, 2:(n+1)]
N_out <- result[, (n+2):(2*n+1)]
D_out <- result[, (2*n+2):(3*n+1)]
Time <- result [,1]
z_values <- z_values



#-----------remin1--------------------------------------------------

derivative <- function(times, y, params) {
  P <- y[1:n]
  N <- y[(n+1):(2*n)]
  De <- y[(2*n+1):(3*n)]
  
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  J_n <- numeric(n + 1)
  J_D <- numeric(n + 1)
  
  for(i in 2:n){
    J_a[i] <-  u*P[i-1]
    J_d[i] <- -D*(P[i]-P[i-1])/DeltaZ
    J_n[i] <- -D*(N[i]-N[i-1])/DeltaZ 
    J_D[i] <- -D*(De[i]-De[i-1])/DeltaZ 
  }
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_n[1] <- 0 
  J_D[1] <- 0 
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  J_D[n + 1] <- u*De[n]
  J_n[n + 1] <- -D * ((Nb-N[n])/DeltaZ)
  
  J = J_a + J_d
  
  jP= (-(J[2:(n+1)] - J[1:n])) / DeltaZ
  jN= -(J_n[2:(n+1)] - J_n[1:n]) / DeltaZ
  jD= -(J_D[2:(n+1)] - J_D[1:n]) / DeltaZ
  
  
  Light <- function(kw,kp,DeltaZ,P,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  
  dP_dt <- numeric(n)
  dN_dt <- numeric(n)
  dD_dt <- numeric(n)
  
  for (i in 1:n){
    dP_dt[i] <- (g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + jP[i]
    dN_dt[i] <- (-g[i]*P[i])+ (remin1*De[i]) + jN[i]
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin1*De[i]) - (u*De[i]) + jD[i]
  }
  
  return(list(c(dP_dt, dN_dt,dD_dt)))
}


params <-
  list(
    n = n,
    DeltaZ = DeltaZ,
    u = u,
    D = D,
    m = m,
    g = g, 
    gmax = gmax, 
    I = I, 
    Nb = Nb, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 100, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out1 <- result[, 2:(n+1)]
N_out1 <- result[, (n+2):(2*n+1)]
D_out1 <- result[, (2*n+2):(3*n+1)]
Time <- result [,1]
z_values <- z_values


#---------remin2------------------------------------------------------
derivative <- function(times, y, params) {
  P <- y[1:n]
  N <- y[(n+1):(2*n)]
  De <- y[(2*n+1):(3*n)]
  
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  J_n <- numeric(n + 1)
  J_D <- numeric(n + 1)
  
  for(i in 2:n){
    J_a[i] <-  u*P[i-1]
    J_d[i] <- -D*(P[i]-P[i-1])/DeltaZ
    J_n[i] <- -D*(N[i]-N[i-1])/DeltaZ 
    J_D[i] <- -D*(De[i]-De[i-1])/DeltaZ 
  }
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_n[1] <- 0 
  J_D[1] <- 0 
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  J_D[n + 1] <- u*De[n]
  J_n[n + 1] <- -D * ((Nb-N[n])/DeltaZ)
  
  J = J_a + J_d
  
  jP= (-(J[2:(n+1)] - J[1:n])) / DeltaZ
  jN= -(J_n[2:(n+1)] - J_n[1:n]) / DeltaZ
  jD= -(J_D[2:(n+1)] - J_D[1:n]) / DeltaZ
  
  
  Light <- function(kw,kp,DeltaZ,P,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  
  dP_dt <- numeric(n)
  dN_dt <- numeric(n)
  dD_dt <- numeric(n)
  
  for (i in 1:n){
    dP_dt[i] <- (g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + jP[i]
    dN_dt[i] <- (-g[i]*P[i])+ (remin2*De[i]) + jN[i]
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin2*De[i]) - (u*De[i]) + jD[i]
  }
  
  return(list(c(dP_dt, dN_dt,dD_dt)))
}


params <-
  list(
    n = n,
    DeltaZ = DeltaZ,
    u = u,
    D = D,
    m = m,
    g = g, 
    gmax = gmax, 
    I = I, 
    Nb = Nb, 
    yield = yield, 
    Hi = Hi,
    Hn = Hn, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 100, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out2 <- result[, 2:(n+1)]
N_out2 <- result[, (n+2):(2*n+1)]
D_out2 <- result[, (2*n+2):(3*n+1)]
Time <- result [,1]
z_values <- z_values


#----------remin3-------------------------------------------------------

derivative <- function(times, y, params) {
  P <- y[1:n]
  N <- y[(n+1):(2*n)]
  De <- y[(2*n+1):(3*n)]
  
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  J_n <- numeric(n + 1)
  J_D <- numeric(n + 1)
  
  for(i in 2:n){
    J_a[i] <-  u*P[i-1]
    J_d[i] <- -D*(P[i]-P[i-1])/DeltaZ
    J_n[i] <- -D*(N[i]-N[i-1])/DeltaZ 
    J_D[i] <- -D*(De[i]-De[i-1])/DeltaZ 
  }
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_n[1] <- 0 
  J_D[1] <- 0 
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  J_D[n + 1] <- u*De[n]
  J_n[n + 1] <- -D * ((Nb-N[n])/DeltaZ)
  
  J = J_a + J_d
  
  jP= (-(J[2:(n+1)] - J[1:n])) / DeltaZ
  jN= -(J_n[2:(n+1)] - J_n[1:n]) / DeltaZ
  jD= -(J_D[2:(n+1)] - J_D[1:n]) / DeltaZ
  
  
  Light <- function(kw,kp,DeltaZ,P,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  
  dP_dt <- numeric(n)
  dN_dt <- numeric(n)
  dD_dt <- numeric(n)
  
  for (i in 1:n){
    dP_dt[i] <- (g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + jP[i]
    dN_dt[i] <- (-g[i]*P[i])+ (remin3*De[i]) + jN[i]
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin3*De[i]) - (u*De[i]) + jD[i]
  }
  
  return(list(c(dP_dt, dN_dt,dD_dt)))
}


params <-
  list(
    n = n,
    DeltaZ = DeltaZ,
    u = u,
    D = D,
    m = m,
    g = g, 
    gmax = gmax, 
    I = I, 
    Nb = Nb, 
    yield = yield, 
    Hi = Hi,
    Hn = Hn, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 100, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out3 <- result[, 2:(n+1)]
N_out3 <- result[, (n+2):(2*n+1)]
D_out3 <- result[, (2*n+2):(3*n+1)]
Time <- result [,1]
z_values <- z_values



plot(P_out[ dim(P_out)[2],],-z_values, type = "l", col = "darkgreen", xlim = c(0,1),lwd=3, main = "Sensitivity test for τ in Phytoplankton", ylab = "Depth[m]", xlab = "[mmol N/m^3]", cex.lab=1.25)
lines(P_out1[ dim(P_out1)[2],], -z_values, type = "l", col = "seagreen", lwd=3)
lines( P_out2[ dim(P_out2)[2],],-z_values, type = "l",col = "green", lwd=3)
lines(P_out3[ dim(P_out3)[2],],-z_values, type = "l",col = "lightgreen",  lwd=3)
legend("bottomright",legend = c("τ = 1.5 ", "τ = 2", "τ = 2.5", "τ = 3"),lty = c(1, 1, 1, 1),col = c("darkgreen", "seagreen", "green", "lightgreen"),lwd = 2)

plot( N_out[ dim(N_out)[2],],-z_values, type = "l", col = "brown", lwd=3, main = "Sensitivity test for τ  in Nutrients", ylab = "Depth[m]", xlab = "[mmol N/m^3]", cex.lab=1.25)  
lines(N_out1[ dim(N_out1)[2],], -z_values,  type = "l", col = "darkorange3", lwd=3)
lines(N_out2[ dim(N_out2)[2],], -z_values,  type = "l", col = "orange", lwd=3)
lines(N_out3[ dim(N_out3)[2],], -z_values,type = "l", col = "yellow4", lwd=3)
legend("topright",legend = c("τ = 1.5 ", "τ = 2", "τ = 2.5", "τ = 3"),lty = c(1, 1, 1, 1),col = c("brown", "darkorange", "orange", "yellow4"),lwd = 2)


plot(D_out[ dim(D_out)[2],],-z_values,  type = "l", col = "black", lwd=3, xlim=c(0,0.08), main = "Sensitivity test for τ in Detritus", ylab = "Depth[m]", xlab = "[mmol N/m^3]", cex.lab=1.25)
lines( D_out1[ dim(D_out1)[2],], -z_values, type = "l", col = "dimgray", lwd=3)
lines(D_out2[ dim(D_out2)[2],], -z_values,  type = "l", col = "darkgrey", lwd=3)
lines( D_out3[ dim(D_out3)[2],], -z_values, type = "l", col = "azure3", lwd=3)
legend("bottomright",legend = c("τ = 1.5 ", "τ = 2", "τ = 2.5", "τ = 3"),lty = c(1, 1, 1, 1),col = c("black", "dimgray", "darkgrey", "azure3"),lwd = 2)



#-----------Changing--remin---------------------------
