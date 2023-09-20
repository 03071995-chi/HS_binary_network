
########################
# Arguments required
# n -- Number of traders
# T -- Number of timepoints
# R -- Number of classes
# theta_true -- true signal
# y -- data
#################################
#load packages
library("nimble")
n = 4
T = 4
R = 4
set.seed(1234)

######################################
###### theta_true simulation #########
######################################
library("smoothmest")
simulate_true_theta<- function(n,T,R){
  theta_true = array(0,c(n,n,T,(R-1)))
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


###########################################
######## data simulation ##################
######## getting the adjacency matrix######
###########################################
simulate_data <- function(n,T,R){
  theta_true = simulate_true_theta(n,T,R)
  y = array(0,c(n,n,T,(R-1)))
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
data = simulate_data(n,T,R)

#########################################
#### EM Update for latent params ###############################
#### ECM Algorithm to get updated theta ########
#########################################
em_update_theta <- function (n, T, R, data){
  a = 0.0005 ## param for hs-like density
  y = data
  ystar = array(0, c(n, n, (T + 1), (R-1)))
  theta = array(0.001, c(n, n, (T + 1), (R-1))) ### a little greater than 0 otherwise initial estimate of nu1 can be infinity
  omega = array(0.5, c(n, n, (T + 2),(R-1)))
  nu =  array(1, c(n, n, (T + 2), (R-1)))
  
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
              
              
              ##estimates of nu1 , omega1 , theta1
              ### estimate of omega1 
              #omega[i,j,1,r] = (1/(2*(theta[i,j,1,r]+C_1)))*tanh((theta[i,j,1,r]+C_1)/2)
              omega[i,j,1,r] = pgdraw.moments(1,(theta[i,j,1,r]+C_1))$mu
              
              ### estimate of nu1 
              nu[i,j,1,r] = (sqrt(a)/(2*pi))*((1/(theta[i,j,1,r])^2)-(1/(a+(theta[i,j,1,r])^2)))
              
              ### estimate of theta1
              num = (ystar[i,j,1,r]*omega[i,j,1,r]/2) + (theta[i,j,2,r]*nu[i,j,2,r]/a)
              denom = (omega[i,j,1,r]/2) + (nu[i,j,1,r]/a) + (nu[i,j,2,r]/a)
              theta[i,j,1,r] =  num/denom
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
              
              ##estimates of nu , omega , theta
              ### estimate of omega 
              #omega[i,j,t,r] = (1/(2*(theta[i,j,t,r]+C_t)))*tanh((theta[i,j,t,r]+C_t)/2)
              omega[i,j,t,r] = pgdraw.moments(1,(theta[i,j,t,r]+C_t))$mu
              
              ### estimate of nu
              nu[i,j,t,r] = (sqrt(a)/(2*pi))*((1/(theta[i,j,t,r]-theta[i,j,(t-1),r])^2)-(1/(a+(theta[i,j,t,r]-theta[i,j,(t-1),r])^2)))
              
              ### estimate of theta
              num_t = (ystar[i,j,t,r]*omega[i,j,1,r]/2) + (theta[i,j,(t-1),r]*nu[i,j,t,r]/a)+(theta[i,j,(t+1),r]*nu[i,j,(t+1),r]/a)
              denom_t = (omega[i,j,t,r]/2) + (nu[i,j,t,r]/a) + (nu[i,j,(t+1),r]/a)
              theta[i,j,t,r] =  num_t/denom_t
            }
          }
        }
      }
    }
  }
  ### assigning 50000 to i=j entries in theta###
  for (i in 1:n){
    for(j in 1:n){
      if(i==j){
        for (r in 1:(R-1)){
          for(t in 1:T){
            theta[i,j,t,r] = 50000
          }
        }
      }
    }
  }
  
  theta_final = theta[,,-(T+1),] ### removing(T+1)entries.
  theta_final = c(theta_final)
  theta_final = theta_final[! theta_final %in% c(50000)] ### removing i=j entries
  
  return(c(theta_final))
}


## Here is my true theta ##

theta_true = simulate_true_theta(n,T,R)
### removing i=j entries###
### this should be done outside the function because simulate_data function needs theta_true of similar dimension to y(data) ### 
for (i in 1:n){
  for(j in 1:n){
    if(i==j){
      for(r in 1:(R-1)){
        for(t in 1:T){
          theta_true[i,j,t,r] = 50000
        }
      }
    }
  }
}
theta_true = c(theta_true)
theta_true = theta_true[! theta_true %in% c(50000)] ### final theta_true, after removing unnecessary values.

#########################################
# #########################################
# em_hs_like <- function(n,T,R,data,theta_true){
#   start = Sys.time()
#   # a = 0.01 ## param for hs-like density
#   # y = data
#   # ystar = array(0, c(n, n, T + 1, R))
#   # theta = array(0.001, c(n, n, T + 1, R)) ### a little greater than 0 otherwise inititial estimate of nu1 will be infinity
#   # omega = array(0.5, c(n, n, T + 2,R))
#   # nu =  array(1, c(n, n, T + 2, R))
#   # ####################################
#   # keep_theta = matrix(0,k, n*n*(T+1)*(R))##(i,j,t,r) combinations can take all such possible combinations and store these n*(n-1)*T*R values in each row of this matrix
#   # keep_theta[2,] = 0.1
#   # diff_matrix = matrix(0,(k-1),n*n*(T+1)*(R))
#   s = array(NA,1000)
#   s[1] = 0.001
#   s_init = em_update_theta(n,T,R,data)
#   s[2] = 
#     k=1
#   while(s[k] > 0.001){
#     keep_theta[k,] = em_update_theta(n,T,R,data)
#     diff_matrix[k,] = abs(keep_theta[(k+1),]-keep_theta[k,])
#     s[k] = max(diff_matrix[k,])
#     k = k + 1
#   }
#   
#   MSE = sum((c(keep_theta[k,]) -c(theta_true))^2)/ (n*(n-1)*T*(R-1))
#   
#   end= Sys.time()
#   execution_time = end-start
#   
#   return(list(keep_theta,MSE))
#   
# }
##################################################### 

run_em_algo<-function(n,T,R,data){
  theta_init = rep(0.001,n*n*(T+1)*(R))
  while(TRUE){
    theta_update = em_update_theta(n,T,R,data)
    s = max(abs(theta_update-theta_init))
    if (s< 0.0000001){
      break
    }
    else{
      theta_init = em_update_theta(n,T,R,data)
      #print("Theta estimates:")
    }
  }
  return(c(theta_update))
}

##### how many zeros in eta #####
l = c(diff(run_em_algo(n,T,R,data)))
p = ifelse(l==0,1,0)
sum(p)
