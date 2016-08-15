# Gillespie-Algorithm
# catalyitcl degragation

########################################################################

Gillespie <- function(t0,Tend,X0,param,propensities){
  # 1. Set initial time t=t0 and number of molecules X(t) = X0
  t 			<- t0
  X 			<- X0
  t_vector 	<- t
  X_matrix 	<- matrix(X,ncol=1)
  while (t<Tend){
    prop 	<- propensities(X,param)
    a0 		<- sum(prop$a)
    # 2. Determine time to next reaction event ...
    tau 	<- -1/a0*log(runif(1,min=0,max=1))
    # ... and type of next reaction event
    y 		<- runif(1,min=0,max=1)
    mu 		<- 0
    a_sum 	<- 0
    while (a_sum<y){
      mu 		<- mu+1
      a_sum 	<- a_sum+prop$a[mu]/sum(prop$a)
    }
    # 4. Update time ...
    t 			<- t+tau
    t_vector 	<- c(t_vector,t)
    # ... and number of molecules
    X 			<- X+prop$nu[,mu]
    X_matrix 	<- cbind(X_matrix,X)
  }	# 5. Go to 2. as long as t < Tend
  return(list(t=t_vector,X=X_matrix))
}

########################################################################

start 		<- c(10,1)				# (A,B)
params 	<- c(c1=0.5,c2=2)

propensitiesUsed <- function(state,param){
  with(as.list(c(state,param)),{
    # propensities
    A 	<- state[1]
    B 	<- state[2]
    a1 	<- c1*A*B
    a2 	<- c2
    # stochiometric matrix
    nu 	<- matrix(c(-1,0,1,0),ncol=2)
    return(list(a=c(a1,a2),nu=nu))
  })
}

G <- Gillespie(t0=0,Tend=200,X0=start,param=params,propensities=propensitiesUsed)
plot(G$t,G$X[1,],'s',xlab='time [sec]',ylab='number of molecules A')

########################################################################
