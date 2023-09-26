
########################
# Arguments required
# n -- Number of traders
# T -- Number of timepoints (train + test)
# R -- Number of classes
# theta_true -- true value of theta
# y -- data
#########################
#load packages
library("nimble")
library("pROC")
library("caret")
library("smoothmest")
library("pgdraw")
#########################

n = 6
T_train = 20
T_test = 5
T = T_train + T_test
R = 4

#set.seed(1234)

######################################
###### theta_true simulation #########
######################################
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
  theta_true = theta_true[,,-(T+1),-R]
  return(theta_true)
}


## Here is my true theta ##
theta_true = simulate_true_theta(n,T,R)

###########################################
######## data simulation ##################
######## getting the adjacency matrix######
###########################################
simulate_data <- function(n,T,R){
  theta_true = simulate_true_theta(n,T,R)
  y = array(0,c(n,n,(T+1),R-1))
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
#########################
data = simulate_data(n,T,R)

### split the data into data_train and data_test
data_train = data[,,1:T_train,]
data_test = data[,,(T_train+1):T,]


#########################################
#### MCMC ###############################
#### Posterior sampling of theta ########
#########################################
update_theta <- function (n, T, R, theta, data){
  y = data
  ystar = array(0, c(n, n, T + 1, R-1))
  theta = array(0, c(n, n, T + 1, R-1))
  zeta = array(1, c(n, n, T + 1, R-1))
  mu = array(0.5, c(n, n,T + 2,R-1))
  omega = array(0.5, c(n, n, T + 2,R-1))
  lambda_sq  = array(0.5, c(n, n, T + 2, R-1))
  nu =  array(0.5, c(n, n, T + 2, R-1))
  tau2 = 0.1
  xi = array(0.5, c(n, n, T + 2, R-1))
  
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
              
              ##step-1: sampling theta1 and omega1
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
              nu[i,j,t,r] = 1/rgamma(1,1,1+(1/lambda_sq[i,j,t,r]))
              
              ## sampling xi
              xi[i,j,t,r] = 1/rgamma(1,1,1+(1/tau2))
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
  
  theta_final = theta[,,-(T+1),] ### removing(T+1) entries.
  #theta_final = c(theta_final)
  #theta_final = theta_final[! theta_final %in% c(50000)] ### removing i=j entries
  
  return(theta_final)
}


#########################################
############## HS-MCMC ##################
#########################################
hs_mcmc <- function(n,T,R,theta,data){
  start = Sys.time()
  ## posterior sampling 
  n.iter = 10000 # number of samples
  burnin = 1000 # 00burn in samples 
  
  keep_theta = matrix(0,n.iter, n*n*T*(R-1))##(i,j,t,r) combinations can take all such possible combinations and store these n*(n-1)*T*R values in each row of this matrix
  
  ## running iterations and storing them row-wise in the keep_theta matrix 
  for (iter in 1:n.iter){
    keep_theta[iter,] = c(update_theta(n,T,R,theta,data))
  }
  bayes_theta = as.vector(apply(keep_theta[(burnin+1):n.iter,],2,mean))
  bayes_theta.q = apply(keep_theta[(burnin+1):n.iter,],2,quantile, probs = c(0.025,0.975))
  #theta_true =  as.vector(theta_true)
  # MSE = sum((bayes_theta - theta_true)^2)/ (n*(n-1)*T*(R-1))
  # adj.MSE <- sum((bayes_theta - theta_true)^2)/sum(theta_true^2)
  out <- NULL
  out$bayes_theta <- bayes_theta
  out$bayes_theta.q <- bayes_theta.q
  # out$MSE <- MSE
  # out$adj.MSE <- adj.MSE
  
  end= Sys.time()
  
  execution_time = end-start
  print(execution_time)
  return(out)
  
}

##################################################
####prediction after estimating the thetas########
##################################################
## based on training data and theta_init(=0) we predict theta's based on hs_mcmc algo
## we arrange those theta values into a 4 dimensional array. Note that the algo gives me theta of the same timepoint length as T_train.
## we subset 1:T_test from the estimated theta and get the estimated theta.
## we predict Y's from the estimated/predicted theta. If P(y[i,j,t,r]=1)>0.5, we assign 1 to y.
## thus we get predicted Y matrices for each time points.

prediction<- function(n,T_train,T_test,R,theta_init,data_train)
{
  theta_init = 0
  theta_for_prediction_of_theta_test = hs_mcmc(n,T_train,R,theta_init,data_train)$bayes_theta
  theta_1 = theta_for_prediction_of_theta_test
  theta_1 = array(theta_1,c(n,n,T_train,R-1))
  theta_estimated = theta_1[,,1:T_test,]

  y_pred = array(0.001,c(n,n,T_test,(R-1)))
  for (i in 1:n){
    for(j in 1:n){
      if(i!=j){
        for (t in 1:T_test){
          for(r in 1:(R-1)){
            y_pred[i,j,t,r] = ifelse(((exp(theta_estimated[i,j,t,1])+exp(theta_estimated[i,j,t,1]+theta_estimated[i,j,t,2]+theta_estimated[i,j,t,3]))/(1+exp(theta_estimated[i,j,t,1])+exp(theta_estimated[i,j,t,2])+exp(theta_estimated[i,j,t,1]+theta_estimated[i,j,t,2]+theta_estimated[i,j,t,3])))>0.5,1,0)
            y_pred[j,i,t,r] = ifelse(((exp(theta_estimated[i,j,t,2])+exp(theta_estimated[i,j,t,1]+theta_estimated[i,j,t,2]+theta_estimated[i,j,t,3]))/(1+exp(theta_estimated[i,j,t,1])+exp(theta_estimated[i,j,t,2])+exp(theta_estimated[i,j,t,1]+theta_estimated[i,j,t,2]+theta_estimated[i,j,t,3])))>0.5,1,0)
            
          }
        }
      }
    }
  }
  
  return(y_pred)
}

#########################################
## Now we have predicted Y's. We remove i=j entries by assigning 0.01 values and dropping it afterwards.
## we also have actual Y's(from test_data). 
## we will plot ROC curves for each time point t_(train+1),...t_test. (5 points for now)
## we also plot auc values for each curve.


pred =  prediction(n,T_train,T_test,R,theta_init,data_train)

for (i in 1:n){
  for(j in 1:n){
    if(i==j){
      for( r in 1:(R-1)){
        for ( t in 1:T_test){
          data_test[i,j,t,r] = 0.001
        }
      }
    }
  }
}


data_test_1 = c(data_test[,,1,])
data_test_1 = data_test_1[! data_test_1 %in% c(0.001)] ## removing i=j entries
pred_1 = c(pred[,,1,])
pred_1 = pred_1[! pred_1 %in% c(0.001)]
roc_score1 = roc(data_test_1 , pred_1) 
#conf_matrix_1 = confusionMatrix(factor(pred_1), factor(data_test_1))
##################################
data_test_2 = c(data_test[,,2,])
data_test_2 = data_test_2[! data_test_2 %in% c(0.001)]
pred_2 = c(pred[,,2,])
pred_2 = pred_2[! pred_2 %in% c(0.001)]
roc_score2 = roc(data_test_2 , pred_2) 
#conf_matrix_2 = confusionMatrix(factor(pred_2), factor(data_test_2))
#####################################
data_test_3 = c(data_test[,,3,])
data_test_3 = data_test_3[! data_test_3 %in% c(0.001)]
pred_3 = c(pred[,,3,])
pred_3 = pred_3[! pred_3 %in% c(0.001)]
roc_score3 = roc(data_test_3 , pred_3) 
#conf_matrix_3 = confusionMatrix(factor(pred_3), factor(data_test_3))
######################################
data_test_4 = c(data_test[,,4,])
data_test_4 = data_test_4[! data_test_4 %in% c(0.001)]
pred_4 = c(pred[,,4,])
pred_4 = pred_4[! pred_4 %in% c(0.001)]
roc_score4 = roc(data_test_4 , pred_4) 
#conf_matrix_4 = confusionMatrix(factor(pred_4), factor(data_test_4))
#######################################
data_test_5 = c(data_test[,,5,])
data_test_5 = data_test_5[! data_test_5 %in% c(0.001)]
pred_5 = c(pred[,,5,])
pred_5 = pred_5[! pred_5 %in% c(0.001)]
roc_score5 = roc(data_test_5 , pred_5) 
#conf_matrix_5 = confusionMatrix(factor(pred_5), factor(data_test_5))


##################################
#####plotting ROC Curves##########
##################################
plot(roc_score1, col = 'red', lty = 2, main = "ROC")
plot(roc_score2, col = 'green', lty = 4, add = TRUE)
plot(roc_score3, col = 'purple', lty = 6, add = TRUE)
plot(roc_score4, col = 'blue', lty = 8, add = TRUE)
plot(roc_score5, col = 'black', lty = 10, add = TRUE)

legend(x = "topleft",  lty = c(2,4,6,8,10), box.lwd = 2, text.font = 4, 
       col= c("red","green","purple","blue","black"), cex=0.7, text.col = "black", 
       legend=c("roc_curve1","roc_curve2","roc_curve3","roc_curve4","roc_curve5"))


##############################
##### plotting AUC values ####
##############################
auc_values = c(roc_score1$auc,roc_score2$auc,roc_score3$auc,roc_score4$auc,roc_score5$auc)
plot(auc_values,xlab = "test_data_timepoints",ylab = "auc_values",type="o",lty = 2, ylim=c(0,1),pch=25)
legend(x = "bottomleft",  lty = 4,  box.lwd = 2, text.font = 4,  cex=0.6, text.col = "black", title = "AUC values from left to right",
       legend=c(round(roc_score1$auc,3),round(roc_score2$auc,3),round(roc_score3$auc,3),round(roc_score4$auc,3),round(roc_score5$auc,3)))





