# Hybrid-stoachstic-deterministic algorithm


library(deSolve)

##################################################

Hybrid <- function(t0,Tend,X0,param,propensities,partitioning,ODEs,prop_threshold=10){
  # t0				starttime
  # Tend				endtime
  # param				parameters
  # propensities		function for calculation of the propensities
  # partitioning		function for calculation of the partitioning
  # ODEs				ODE-Funktion, which is handedover the the ODE solver
  # prop_threshold	propensitiesThreshold, small -> stochastic; big -> deterministic
  
  param <- c(param,prop_threshold=prop_threshold)
  
  # 1. Set initial time t=t0 and number of molecules X(t) = X0
  t <- t0
  X <- X0
  t_vector <- t
  X_matrix <- matrix(X,ncol=1)
  counter <- 0
  
  while (t<Tend){
    
    # Determine partitioning
    partition <- partitioning(X,param)
    
    # 2. Draw xi~Exp(1)...
    xi <- -log(runif(1,min=0,max=1))
    
    # ... and solve ODEs
    initial <- c(X,g=0)								# g(s0)=0
    param2 <- list(param=c(param,xi=xi),partition=partition)
    timepoints <- seq(t,Tend+0.1,by=0.1)			# Simulation till Tend, in case no reaction occured till then
    out <- ode(y=initial,time=timepoints,func=ODEs,parms=param2,rootfun=rootfun)
    X <- out[nrow(out),2:(ncol(out)-1)]				# update number of molekueles
    X[X<0] <- 0
    
    # 3. Determine next reaction type
    if (out[nrow(out),'time']<Tend){				# in case reaction occured before reaching Tend
      prop <- propensities(X,param2)				# Propensities at the time the reaction occured
      y <- runif(1,min=0,max=1)
      mu <- 0
      fractions <- cumsum(prop$a)/sum(prop$a)
      for (i in 1:length(fractions)){
        if (y<=fractions[i]){
          mu <- i
          break
        }
      }
      X <- X+prop$v[,mu]							# update number molekueles
      X[X<0] <- 0
    }
    # 4. Update time ...
    t <- out[nrow(out),'time']
    t_vector <- c(t_vector,out[,'time'],t)
    # ... and number of molecules
    X_matrix <- cbind(X_matrix,t(out[,2:(ncol(out)-1)]),X)
  }	# 5. Go to 2. as long as t < Tend
  return(list(t=t_vector,X=X_matrix))
}

##################################################

# Rootfun fuer ODE-Solver: stops if rootfun==0
rootfun <- function(time,state,parameters){
  return(state[length(state)]-parameters$param['xi'])
}

# function for calculation of the partitioning
partitioning <- function(state,param){
  with(as.list(c(state,param)),{
    a    <- rep(0,12)
    a[1] <- lambdaT
    a[2] <- deltaT*TU
    a[3] <- CLT*TU*VI
    a[4] <- betaT*TU*VI
    a[5] <- CL*VI
    a[6] <- deltaT1*T1
    a[7] <- deltaPICT*T1
    a[8] <- kT*T1
    a[9] <- deltaT2*T2
    a[10]<- NT*T2
    a[11]<- (NThut-NT)*T2
    a[12]<- CL*VNI
    partition <- rep(0,length(a))
    partition[a<param['prop_threshold']] <- 1 		# 1: stochastic, 0: deterministic
    return(partition)
  })
}

# function for calculation of the propensities
propensities <- function(state,parameters){
  param		<- parameters$param
  partition 	<- parameters$partition
  with(as.list(c(state,param,partition)),{
    # propensities
    TU <- state[1]
    T1 <- state[2]
    T2 <- state[3]
    VI <- state[4]
    VNI<- state[5]
    a1 <- lambdaT		*partition[1]
    a2 <- deltaT*TU		*partition[2]
    a3 <- CLT*TU*VI		*partition[3]
    a4 <- betaT*TU*VI	*partition[4]
    a5 <- CL*VI			*partition[5]
    a6 <- deltaT1*T1	*partition[6]
    a7 <- deltaPICT*T1	*partition[7]
    a8 <- kT*T1			*partition[8]
    a9 <- deltaT2*T2	*partition[9]
    a10<- NT*T2			*partition[10]
    a11<- (NThut-NT)*T2	*partition[11]
    a12<- CL*VNI		*partition[12]
    # stochiometric vectors
    v  <- matrix(c(
      1, 0, 0, 0, 0,		# mu=1
      -1, 0, 0, 0, 0,		# mu=2
      0, 0, 0,-1, 0, 	# mu=3
      -1, 1, 0,-1, 0,		# mu=4
      0, 0, 0,-1, 0, 	# mu=5
      0,-1, 0, 0, 0,		# mu=6
      1,-1, 0, 0, 0,		# mu=7
      0,-1, 1, 0, 0,		# mu=8
      0, 0,-1, 0, 0,		# mu=9
      0, 0, 0, 1, 0,		# mu=10
      0, 0, 0, 0, 1,		# mu=11
      0, 0, 0, 0,-1		# mu=12
    ),byrow=FALSE,ncol=12)
    return(list(a=unname(c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)),v=v))
  })
}

# ODE-Funktion, which is handedover the the ODE solver
ODEs <- function(t,state,parameters){
  param		<- parameters$param
  partition 	<- parameters$partition
  with(as.list(c(state,param,partition)),{
    p 		<- 1-partition
    TU  	<- state[1]
    T1  	<- state[2]
    T2  	<- state[3]
    VI  	<- state[4]
    VNI 	<- state[5]
    dTU 	<- lambdaT*p[1] -deltaT*TU*p[2] -betaT*VI*TU*p[4] +deltaPICT*T1*p[7]
    dT1 	<- betaT*VI*TU*p[4] -deltaT1*T1*p[6] -deltaPICT*T1*p[7] -kT*T1*p[8]
    dT2		<- kT*T1*p[8] -deltaT2*T2*p[9]
    dVI 	<- NT*T2*p[10] -CL*VI*p[5] -CLT*VI*TU*p[3] -betaT*VI*TU*p[4]
    dVNI 	<- (NThut-NT)*T2*p[11] -CL*VNI*p[12]
    dg 		<- sum(propensities(state,parameters)$a)
    return(list(c(dTU,dT1,dT2,dVI,dVNI,dg)))
  })
}

# parameter
param <- c( lambdaT = 2e9,
            betaT	= 8e-12,
            NThut	= 1000,
            deltaT	= 0.02,
            deltaT1	= 0.02,
            deltaT2	= 1,
            CL		= 23,
            kT		= 0.35,
            deltaPICT	= 0.35,
            NT		= 0.67*1000,
            CLT		= (1/0.33-1)*8e-12)

##################################################

nr <- 100						# number of realisiations
infection_threshold <- 1		# VI > threshold at t_end defines succesful infection
VI0 <- c(1,10,100,1000)			# Viruslast VI(0)
prob <- rep(NA,length(VI0))
for (j in 1:length(VI0)){
  X0 <- c(TU=2e9/0.02,T1=0,T2=0,VI=VI0[j],VNI=0)
  inf_nr <- 0					# number of realisiations with succesful infection
  for (i in 1:nr){
    G <- Hybrid(t0=0,Tend=60,X0=X0,param=param,propensities=propensities,partitioning=partitioning,ODEs=ODEs)
    if (G$X['VI',ncol(G$X)] > infection_threshold) inf_nr <- inf_nr+1
  }
  prob[j] <- inf_nr/nr
  cat('Probability of infection with viral load VI =',VI0[j],'is',prob[j],'.\n')
}
