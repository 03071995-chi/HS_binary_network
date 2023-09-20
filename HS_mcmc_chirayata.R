
########################
# Arguments required
# n -- Number of traders
# T -- Number of timepoints
# R -- Number of classes
# theta_true -- true signal (for computing MSE)
# y -- data
#################################
#load packages
library("nimble")
n = 2
T = 2
R = 4
set.seed(1234)

######################################
###### theta_true simulation #########
######################################
library("smoothmest")
simulate_true_theta<- function(n,T,R){
  theta_true = array(0,c(n,n,(T+1),R))
  for (i in 1:n){
    for(j in 1:n){
      if(i!=j){
        for (t in 2:T){
          for (r in 1:(R-1)){
            theta_true[i,j,t,r] = rdoublex(1,mu=theta_true[i,j,t-1,r],lambda=1/12)
          }
        }
      }
    }
  }
  return(theta_true)
}

#theta_true = simulate_theta_true(n,T,R)

###########################################
######## data simulation ##################
######## getting the adjacency matrix######
###########################################
simulate_data <- function(n,T,R){
  theta_true = simulate_true_theta(n,T,R)
  y = array(0,c(n,n,(T+1),R))
  for (i in 1:n){
    for(j in 1:n){
      if(i!=j){
        for(r in 1:(R-1)){
          for(t in 1:T){
            c  = (1+exp(theta_true[i,j,t,1])+exp(theta_true[i,j,t,2])+(exp(theta_true[i,j,t,1])+(exp(theta_true[i,j,t,2]))+(exp(theta_true[i,j,t,3]))))
            p1 = 1/c
            p2 = exp(theta_true[i,j,t,1])/c
            p3 = exp(theta_true[i,j,t,2])/c
            p4 = 1-(p1+p2+p3)
            
            #generate a uniform(0,1)
            
            u = runif(1,0,1)
            if(u<=p1){
              y[i,j,t,r] = 0
              y[j,i,t,r] = 0
            }else if ((p1< u ) & (u<=p1+p2)){
              y[i,j,t,r] = 1
              y[j,i,t,r] = 0
            }else if ((p1+p2<u) & (u<=p1+p2+p3)){
              y[i,j,t,r] = 0
              y[j,i,t,r] = 1
            }else{
              y[i,j,t,r] = 1
              y[j,i,t,r] = 1
            }
          }  
        }
      }
    }
  }
  return(y)
}

#### here is my data ####
##################
data = simulate_data(n,T,R)

#########################################
#### MCMC ###############################
#### Posterior sampling of theta ########
#########################################
library("pgdraw")
update_theta <- function (n, T, R, data){
  y = data
  ystar = array(0, c(n, n, T + 1, R))
  theta = array(0, c(n, n, T + 1, R))
  zeta = array(1, c(n, n, T + 1, R))
  mu = array(0.5, c(n, n,T + 2,R))
  omega = array(0.5, c(n, n, T + 2,R))
  lambda_sq  = array(0.5, c(n, n, T + 2, R))
  nu =  array(0.5, c(n, n, T + 2, R))
  tau2 = 0.01
  xi = array(0.5, c(n, n, T + 2, R))
  
  for (i in 1:n){
    for(j in 1:n){
      if(i!=j){
        for (r in 1:(R-1)){
          for(t in 1:T){
            
            ##################### for t=1 ###################
            if (t == 1) {
              C_1 = 0
              if(r == 1){
                C_1 = log(1+exp(theta[i,j,1,r+1]+theta[i,j,1,r+2])) - log(1+exp(theta[i,j,1,r+1]))
                ystar[i,j,1,r] = ((y[i,j,1,r] -0.5)/omega[i,j,1,r]) -C_1
              }else if(r == 2){
                C_1 = log(1+exp(theta[i,j,1,r-1]+theta[i,j,1,r+1])) - log(1+exp(theta[i,j,1,r-1]))
                ystar[i,j,1,r] = ((y[j,i,1,r] -0.5)/omega[i,j,1,r]) - C_1
              }else {
                C_1 = log(exp(theta[i,j,1,r-2]+theta[i,j,1,r-1] + theta[i,j,t,r])) - log(1+exp(theta[i,j,1,r-2]) + exp(theta[i,j,1,r-1]))
                ystar[i,j,1,r] = ((y[i,j,1,r]*y[j,i,1,r] -0.5)/omega[i,j,1,r]) - C_1
              }
              
              zeta[i,j,1,r] = omega[i,j,1,r] + (omega[i,j,1,r]/(tau2*lambda_sq[i,j,1,r])) + (omega[i,j,2,r]/(tau2*lambda_sq[i,j,2,r]))
              mu[i,j,1,r] = (1/zeta[i,j,1,r])* ((ystar[i,j,1,r]*omega[i,j,1,r])+(theta[i,j,2,r]*omega[i,j,2,r])/(tau2*lambda_sq[i,j,2,r]))
              
              ##sampling theta1 and omega1
              theta[i,j,1,r] =  rnorm(1,mu[i,j,1,r],(1/zeta[i,j,1,r]))
              omega[i, j, 1, r] <- pgdraw(1, (theta[i,j,1,r]+ C_1))
            }
            ####################### for t>1 #################
            if (t > 1) {
              C_t = 0
              if(r == 1){
                C_t = log(1+exp(theta[i,j,t,r+1]+theta[i,j,t,r+2])) - log(1+exp(theta[i,j,t,r+1]))
                ystar[i,j,t,r] = ((y[i,j,t,r] -0.5)/omega[i,j,t,r]) -C_t
              }else if(r==2){
                C_t = log(1+exp(theta[i,j,t,r-1]+theta[i,j,t,r+1])) - log(1+exp(theta[i,j,t,r-1]))
                ystar[i,j,t,r] = ((y[j,i,t,r] -0.5)/omega[i,j,t,r]) - C_t
              }else {
                C_t = log(exp(theta[i,j,t,r-2]+theta[i,j,t,r-1]+theta[i,j,t,r])) - log(1+exp(theta[i,j,t,r-2])+exp(theta[i,j,t,r-1]))
                ystar[i,j,t,r] = ((y[i,j,t,r]*y[j,i,t,r] -0.5)/omega[i,j,t,r]) - C_t 
              }
              #sampling omega
              omega[i,j,t,r] = pgdraw(1,(theta[i,j,t,r]+C_t))
              
              ## sampling theta's
              
              ## step-2: sampling rest of the theta's except theta1
              zeta[i,j,t,r] = omega[i,j,t,r]+(omega[i,j,t,r]/(tau2*lambda_sq[i,j,t,r]))+(omega[i,j,t+1,r]/(tau2*lambda_sq[i,j,t+1,r]))
              mu[i,j,t,r] = (1/zeta[i,j,t,r])* ((ystar[i,j,t,r]*omega[i,j,t,r])+(theta[i,j,t-1,r]*omega[i,j,t,r])/tau2*lambda_sq[i,j,t,r]+(theta[i,j,t+1,r]*omega[i,j,t+1,r])/tau2*lambda_sq[i,j,t+1,r])
              theta[i,j,t,r] =  rnorm(1,mu[i,j,t,r],(1/zeta[i,j,t,r]))
              
              ## sampling lambda_sq and nu
              lambda_sq[i,j,t,r] = 1/rgamma(1,1,(1/nu[i,j,t,r]) +(theta[i,j,t+1,r]-theta[i,j,t,r])^2 * omega[i,j,t+1,r]/(2*tau2))
              nu[i,j,t,r] = 1/rgamma(1,1,1+1/lambda_sq[i,j,t,r])
              
              ## sampling xi
              xi[i,j,t,r] = 1/rgamma(1,1,1+1/tau2)
            }
          }
        }
      }
    }
  }
  
  return(c(theta))
}


## Here is my true theta ##
theta_true = simulate_true_theta(n,T,R)

#########################################
############## HS-MCMC ##################
#########################################
hs_mcmc <- function(n,T,R,theta_true,data){
  start = Sys.time()
  ## posterior sampling 
  n.iter = 550 # number of samples
  burnin = 100 # burn in samples 
 
  keep_theta = matrix(0,n.iter, n*n*(T+1)*R)##(i,j,t,r) combinations can take all such possible combinations and store these n*(n-1)*T*R values in each row of this matrix
  
  ## running iterations and storing them row-wise of the keep_theta matrix 
  for (iter in 1:n.iter){
    keep_theta[iter,] = update_theta(n,T,R,data)
  }
  bayes_theta = as.vector(apply(keep_theta[(burnin+1):n.iter,],2,mean))
  bayes_theta.q = apply(keep_theta[(burnin+1):n.iter,],2,quantile, probs = c(0.025,0.975))
  #print(bayes_theta)
  theta_true =  as.vector(theta_true)
  MSE = sum((bayes_theta - theta_true)^2)/ (n*(n-1)*T*(R-1))
  adj.MSE <- sum((bayes_theta - theta_true)^2)/sum(theta_true^2)
  
  out <- NULL
  out$bayes_theta <- bayes_theta
  out$bayes_theta.q <- bayes_theta.q
  out$MSE <- MSE
  out$adj.MSE <- adj.MSE
  
  end= Sys.time()
  
  execution_time = end-start
  
  return(list(out,execution_time))
  
}
##################################################### 
# 
# fun1<- function(x){
#   (1-exp(-x))/(2*sqrt(pi)*x^(3/2))
# }
# 
# fun2<-function(x){
#   2/(pi*(1+x^2))
# }
# 
# x = seq(0,5,0.01)
# y1 = sapply(x,fun1)
# y2 = sapply(x,fun2)
# 
# plot(x,y1,ylim = c(0,1),type="l",col="red")
# lines(x,y2,ylim=c(0,1),type="l",col = "blue")
# legend(x = "topright", lty = c(4,6), text.font = 4, 
#        col= c("blue","red"), 
#        legend=c("half cauchy", "hs-like"))
# 
# #integrate(fun,lower=0.00001,upper= 10000)

u = matrix(1,1e4,1)
x = pgdraw(1,0.5*u)
mean(x)
var(x)
pgdraw.moments(1,0.5)

x = pgdraw(2,2*u)
mean(x)
var(x)
pgdraw.moments(2,2)
