library(nlme)

#Function repVar() fits a weighted regression model to the input data, 
#and then returns the variance structure estimated from the model.
#Currently only varExp() function is used as the variance structure
varRep <- function(tf1_x, tf1_y, plot = F, xlab = "tf1_x", ylab = "tf1_y"){
  data_fit <- data.frame(x = tf1_x, y = tf1_y)
  n <- nrow(data_fit)
  
  #Use weighted least square as the fitting criteria
  #wls1 is a linear model
  wls1 <- gls(y ~ x, data = data_fit, weights = varExp(), 
              method = "ML") #Use "ML" so we can do anova test
  wls1_fit <- predict(wls1)
  t_wls1 <- as.numeric(wls1$modelStruct)
  sigma_wls1 <- wls1$sigma
  var_wls1 <- exp(2*t_wls1*data_fit$x)*(sigma_wls1)^2 #variance structure
  #This band is not a confidence interval. It just shows the variance structure
  wls1_l <- wls1_fit - sqrt(var_wls1)
  wls1_u <- wls1_fit + sqrt(var_wls1)
  
  #wls2 is a quadratic model
  wls2 <- gls(y ~ I(x^2) + x, data = data_fit, weights = varExp(form = ~x), 
              method = "ML")
  wls2_fit <- predict(wls2)
  t_wls2 <- as.numeric(wls2$modelStruct)
  sigma_wls2 <- wls2$sigma
  var_wls2 <- exp(2*t_wls2*data_fit$x)*(sigma_wls2)^2 #variance structure
  wls2_l <- wls2_fit - sqrt(var_wls2)
  wls2_u <- wls2_fit + sqrt(var_wls2)
  
  #Use anova test to determine which model is better
  #Anova may be invalid since the variance structure is different
  #Maybe there is a better way to compare the two models
  res_anova <- anova(wls1, wls2)
  if(res_anova$'p-value'[2] < 0.05){ #Currently we use a fixed value here
    t_str <- t_wls2
    sigma_str <- sigma_wls2
  } else{
    t_str <- t_wls1
    sigma_str <- sigma_wls1
  }
  res_str <- list(t = t_str, sigma = sigma_str, anova.test = res_anova)
  
  #Plot the two weighted regression models with prediction bands
  if(plot == T){
    plot(data_fit, col = rgb(0,0,0,0.3), cex=.3, xlab = xlab, ylab = ylab)
    lines(sort(data_fit$x), wls1_fit[order(data_fit$x)], lwd = 2, col = 2)
    lines(sort(data_fit$x), wls1_l[order(data_fit$x)], lty = 2, col = 2)
    lines(sort(data_fit$x), wls1_u[order(data_fit$x)], lty = 2, col = 2)
    lines(sort(data_fit$x), wls2_fit[order(data_fit$x)], lwd = 2, col = 3)
    lines(sort(data_fit$x), wls2_l[order(data_fit$x)], lty = 2, col = 3)
    lines(sort(data_fit$x), wls2_u[order(data_fit$x)], lty = 2, col = 3)
    if(res_anova$`p-value`[2] < 0.05){
      lines(sort(data_fit$x), wls2_fit[order(data_fit$x)], lwd = 2, col = 4)
      lines(sort(data_fit$x), wls2_l[order(data_fit$x)], lty = 2, col = 4)
      lines(sort(data_fit$x), wls2_u[order(data_fit$x)], lty = 2, col = 4)
    } else{
      lines(sort(data_fit$x), wls1_fit[order(data_fit$x)], lwd = 2, col = 4)
      lines(sort(data_fit$x), wls1_l[order(data_fit$x)], lty = 2, col = 4)
      lines(sort(data_fit$x), wls1_u[order(data_fit$x)], lty = 2, col = 4)
    }
  }
  
  return(res_str)
}

#Function classTF() detects the binding preference between two different
#transcription factors for each DNA site. The function takes in the binding
#intensity data of the two transcription factors obtained from PBM 
#experiments, fits a weighted regression model to the data, evaluates the
#deviation of each point from a hypothetical pair of TF replicas by adapting
#a given variance structure of TF replicas, and at last returns the binding
#specificity by a label. 
# >0: preferred by TF1, 0: no significant preference, <0: preferred by TF2

classTF <- function(tf1_x, tf2_y, t, sigma, plot = F, 
                    xlab = "tf1_x", ylab = "tf1_y", level = 0.99){
  data_fit <- data.frame(x = tf1_x, y = tf2_y)
  n <- nrow(data_fit)
  var_adapted <- exp(2*t*data_fit$x)*sigma^2 #adapted variance structure
  
  #wls1 is a linear model
  wls1 <- gls(y ~ x, data = data_fit, weights = varExp(), 
              method = "ML")
  wls1_fit <- predict(wls1)
  t_wls1 <- as.numeric(wls1$modelStruct)
  sigma_wls1 <- wls1$sigma
  var_wls1 <- exp(2*t_wls1*data_fit$x)*(sigma_wls1)^2 #model variance structure
  #var_str is the variance of the prediction, a sum of the uncertainty in 
  #the regression line and in the new hypothetical replica prediction
  var_str <- var_adapted + var_wls1*
    (1/n + (data_fit$x-mean(data_fit$x))^2/((n-1)*var(data_fit$x)))
  wls1_l <- wls1_fit - qnorm(0.5 + level/2)*sqrt(var_str)
  wls1_u <- wls1_fit + qnorm(0.5 + level/2)*sqrt(var_str)
  #the function form of wls1 and its prediction band
  wls1_fun <- function(x)
    wls1$coefficients[2]*x + wls1$coefficients[1]
  
  wls1_fun_l <- function(x)
    wls1$coefficients[2]*x + wls1$coefficients[1] - 
    qnorm(0.5 + level/2)*
    sqrt(exp(2*t*x)*sigma^2) # + exp(2*t_wls1*x)*(sigma_wls1)^2*
  #         (1/n + (x-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
  #The part under square root represents the variance of y given any x
  #We could use this information to calculate how much standard deviation
  #there is between y and the regression line
  
  wls1_fun_u <- function(x)
    wls1$coefficients[2]*x + wls1$coefficients[1] + 
    qnorm(0.5 + level/2)*
    sqrt(exp(2*t*x)*sigma^2) # + exp(2*t_wls1*x)*(sigma_wls1)^2*
  #         (1/n + (x-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
  
  #wls2 is a quadratic model
  wls2 <- gls(y ~ I(x^2) + x, data = data_fit, weights = varExp(form = ~x), 
              method = "ML")
  wls2_fit <- predict(wls2)
  t_wls2 <- as.numeric(wls2$modelStruct)
  sigma_wls2 <- wls2$sigma
  var_wls2 <- exp(2*t_wls2*data_fit$x)*(sigma_wls2)^2 #model variance structure
  #var_str is the variance of the prediction, a sum of the uncertainty in 
  #the regression line and in the new hypothetical replica prediction
  var_str <- var_adapted + var_wls2*
    (1/n + (data_fit$x-mean(data_fit$x))^2/((n-1)*var(data_fit$x)))
  wls2_l <- wls2_fit - qnorm(0.5 + level/2)*sqrt(var_str)
  wls2_u <- wls2_fit + qnorm(0.5 + level/2)*sqrt(var_str)
  #the function form of wls2 and its prediction band
  wls2_fun <- function(x)
    wls2$coefficients[2]*x^2 + wls2$coefficients[3]*x + 
    wls2$coefficients[1]
  
  wls2_fun_l <- function(x)
    wls2$coefficients[2]*x^2 + wls2$coefficients[3]*x + 
    wls2$coefficients[1] - 
    qnorm(0.5 + level/2)*
    sqrt(exp(2*t*x)*sigma^2) # + exp(2*t_wls2*x)*(sigma_wls2)^2*
    #       (1/n + (x-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
  
  wls2_fun_u <- function(x)
    wls2$coefficients[2]*x^2 + wls2$coefficients[3]*x + 
    wls2$coefficients[1] + 
    qnorm(0.5 + level/2)*
    sqrt(exp(2*t*x)*sigma^2) #+ exp(2*t_wls2*x)*(sigma_wls2)^2*
    #       (1/n + (x-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
  
  #Use anova test to determine which model is better
  res_anova <- anova(wls1, wls2)
  pref <- numeric(n) #a vector for recording the preference
  pref_score <- numeric(n) 
  #preference score is defined as the number of standard deviations from the regression line
  sd <- numeric(n) #for storing standard deviations
  if(res_anova$'p-value'[2] < 0.05){ #use polynomial model
    for(i in 1:n){
      sd[i] <- sqrt(exp(2*t*data_fit$x[i])*sigma^2) # + exp(2*t_wls2*data_fit$x[i])*(sigma_wls2)^2*
      #                (1/n + (data_fit$x[i]-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
      if(data_fit$y[i] < wls2_fun_l(data_fit$x[i])){
        pref_score[i] <- (wls2_fun(data_fit$x[i]) - data_fit$y[i])/sd[i]   
      } else{
        if(data_fit$y[i] > wls2_fun_u(data_fit$x[i])){
          pref_score[i] <- (wls2_fun(data_fit$x[i]) - data_fit$y[i])/sd[i]   
        }
      }  
    }
  }else{ #use linear model
    for(i in 1:n){
      sd[i] <- sqrt(exp(2*t*data_fit$x[i])*sigma^2) # + exp(2*t_wls1*data_fit$x[i])*(sigma_wls1)^2*
      #                (1/n + (data_fit$x[i]-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
      if(data_fit$y[i] < wls1_fun_l(data_fit$x[i])){
        pref_score[i] <- (wls1_fun(data_fit$x[i]) - data_fit$y[i])/sd[i]
      } else{
        if(data_fit$y[i] > wls1_fun_u(data_fit$x[i])){
          pref_score[i] <- (wls1_fun(data_fit$x[i]) - data_fit$y[i])/sd[i]
        }
      }
    }
  }
  
  #Plot the two weighted regression models with prediction bands
  if(plot == T){
    plot(data_fit, col = rgb(0,0,0,0.3), xlab = xlab, ylab = ylab,cex=.3)
    lines(sort(data_fit$x), wls1_fit[order(data_fit$x)], lwd = 2, col = 2)
    lines(sort(data_fit$x), wls1_l[order(data_fit$x)], lty = 2, col = 2)
    lines(sort(data_fit$x), wls1_u[order(data_fit$x)], lty = 2, col = 2)
    lines(sort(data_fit$x), wls2_fit[order(data_fit$x)], lwd = 2, col = 3)
    lines(sort(data_fit$x), wls2_l[order(data_fit$x)], lty = 2, col = 3)
    lines(sort(data_fit$x), wls2_u[order(data_fit$x)], lty = 2, col = 3)
    if(res_anova$`p-value`[2] < 0.05){
      curve(wls2_fun, col = 4, add = T)
      curve(wls2_fun_l, col = 4, lty = 2, add = T)
      curve(wls2_fun_u, col = 4, lty = 2, add = T)
    } else{
      curve(wls1_fun, col = 4, add = T)
      curve(wls1_fun_l, col = 4, lty = 2, add = T)
      curve(wls1_fun_u, col = 4, lty = 2, add = T)
    }
  }
  
  return(pref_score) #returns pref as a vector
}

#Function class_predTF() detects the binding preference between two different
#transcription factors for each DNA site. The function takes in several sources 
#of input including: 1) the occupancy data of the two transcription 
#factors obtained from PBM experiments, 2) the predicted occupancy matrix for 
#each TF for all genomic regions, 3) the variance structure of TF replicates.
#The function fits a weighted regression model to the measured gcPBM data, 
#evaluates the deviation of each point from a hypothetical pair of TF replicas 
#by adapting a given variance structure of TF replicas, and at last returns the 
#occupancy matrix by preference score. 
# First three columns of Pred_matrix_tf are genomic coordinates and ChIP-seq signal 
# information, remove first before assigning preference
#<0 = preferred by TF1, 0 = no significant preference, >0 = preferred by TF2

class_predTF <- function(tf1_x, tf2_y, Pred_matrix_tf1, Pred_matrix_tf2, t, sigma, plot = F, 
                         xlab = "tf1_x", ylab = "tf2_y", level = 0.99){
  # Convert the Pred_matrix to make easy for assigning preference
  Pred_genomic_coord = Pred_matrix_tf1[,1:3] # Keep info about genomic coordinates predicted
  Pred_matrix_tf1 = as.matrix(Pred_matrix_tf1[,-(1:3)]) # Remove first 3 columns with information of genomic coordinates
  Pred_matrix_tf2 = as.matrix(Pred_matrix_tf2[,-(1:3)])
  n_pred_pos = length(Pred_matrix_tf1[1,])
  Pred_matrix_combined = data.frame(x = as.vector(t(Pred_matrix_tf1),mode='numeric'), y = as.vector(t(Pred_matrix_tf2),mode='numeric')) # Combine the predictions into one data frame
  n_line = nrow(Pred_matrix_combined)
  
  data_fit <- data.frame(x = tf1_x, y = tf2_y) #occupancy data for paralogous TFs measured from gcPBM to fit
  n <- nrow(data_fit)
  var_adapted <- exp(2*t*data_fit$x)*sigma^2 #adapted variance structure
  
  #wls1 is a linear model
  wls1 <- gls(y ~ x, data = data_fit, weights = varExp(), 
              method = "ML")
  wls1_fit <- predict(wls1)
  t_wls1 <- as.numeric(wls1$modelStruct)
  sigma_wls1 <- wls1$sigma
  var_wls1 <- exp(2*t_wls1*data_fit$x)*(sigma_wls1)^2 #model variance structure
  #var_str is the variance of the prediction, a sum of the uncertainty in 
  #the regression line and in the new hypothetical replica prediction
  var_str <- var_adapted + var_wls1*
    (1/n + (data_fit$x-mean(data_fit$x))^2/((n-1)*var(data_fit$x)))
  wls1_l <- wls1_fit - qnorm(0.5 + level/2)*sqrt(var_str)
  wls1_u <- wls1_fit + qnorm(0.5 + level/2)*sqrt(var_str)
  #the function form of wls1 and its prediction band
  wls1_fun <- function(x)
    wls1$coefficients[2]*x + wls1$coefficients[1]
  
  wls1_fun_l <- function(x)
    wls1$coefficients[2]*x + wls1$coefficients[1] - 
    qnorm(0.5 + level/2)*
    sqrt(exp(2*t*x)*sigma^2) # + exp(2*t_wls1*x)*(sigma_wls1)^2*
  #         (1/n + (x-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
  #The part under square root represents the variance of y given any x
  #We could use this information to calculate how much standard deviation
  #there is between y and the regression line
  
  wls1_fun_u <- function(x)
    wls1$coefficients[2]*x + wls1$coefficients[1] + 
    qnorm(0.5 + level/2)*
    sqrt(exp(2*t*x)*sigma^2) # + exp(2*t_wls1*x)*(sigma_wls1)^2*
  #         (1/n + (x-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
  
  #wls2 is a quadratic model
  wls2 <- gls(y ~ I(x^2) + x, data = data_fit, weights = varExp(form = ~x), 
              method = "ML")
  wls2_fit <- predict(wls2)
  t_wls2 <- as.numeric(wls2$modelStruct)
  sigma_wls2 <- wls2$sigma
  var_wls2 <- exp(2*t_wls2*data_fit$x)*(sigma_wls2)^2 #model variance structure
  #var_str is the variance of the prediction, a sum of the uncertainty in 
  #the regression line and in the new hypothetical replica prediction
  var_str <- var_adapted + var_wls2*
    (1/n + (data_fit$x-mean(data_fit$x))^2/((n-1)*var(data_fit$x)))
  wls2_l <- wls2_fit - qnorm(0.5 + level/2)*sqrt(var_str)
  wls2_u <- wls2_fit + qnorm(0.5 + level/2)*sqrt(var_str)
  #the function form of wls2 and its prediction band
  wls2_fun <- function(x)
    wls2$coefficients[2]*x^2 + wls2$coefficients[3]*x + 
    wls2$coefficients[1]
  
  wls2_fun_l <- function(x)
    wls2$coefficients[2]*x^2 + wls2$coefficients[3]*x + 
    wls2$coefficients[1] - 
    qnorm(0.5 + level/2)*
    sqrt(exp(2*t*x)*sigma^2) # + exp(2*t_wls2*x)*(sigma_wls2)^2*
  #         (1/n + (x-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
  
  wls2_fun_u <- function(x)
    wls2$coefficients[2]*x^2 + wls2$coefficients[3]*x + 
    wls2$coefficients[1] + 
    qnorm(0.5 + level/2)*
    sqrt(exp(2*t*x)*sigma^2) # + exp(2*t_wls2*x)*(sigma_wls2)^2*
  #         (1/n + (x-mean(data_fit$x))^2/((n-1)*var(data_fit$x))))
  
  #Use anova test to determine which model is better
  res_anova <- anova(wls1, wls2)
  pref_score <- numeric(n_line) #a vector for recording the preference
  #preference score is defined as the number of standard deviations from the regression line
  sd <- numeric(n_line) #for storing standard deviations
  if(res_anova$'p-value'[2] < 0.05){ #use polynomial model
    for(i in 1:n_line){
      if(Pred_matrix_combined[i,1] == 0){ # No preference if no prediction occurred in this position
        pref_score[i] <- 0
      }else{ # Assign preference to predicted values
        sd[i] <- sqrt(exp(2*t*Pred_matrix_combined$x[i])*sigma^2) # + exp(2*t_wls2*Pred_matrix_combined$x[i])*(sigma_wls2)^2*
        #                (1/n + (Pred_matrix_combined$x[i]-mean(Pred_matrix_combined$x))^2/((n-1)*var(Pred_matrix_combined$x))))
        if(Pred_matrix_combined$y[i] < wls2_fun_l(Pred_matrix_combined$x[i])){
          pref_score[i] <- (wls2_fun(Pred_matrix_combined$x[i]) - Pred_matrix_combined$y[i])/sd[i]   
        } else{
          if(Pred_matrix_combined$y[i] > wls2_fun_u(Pred_matrix_combined$x[i])){
            pref_score[i] <- (wls2_fun(Pred_matrix_combined$x[i]) - Pred_matrix_combined$y[i])/sd[i]   
          }
        }
      }
    }
  }else{ #use linear model
    for(i in 1:n_line){
      if(Pred_matrix_combined[i,1] == 0){ # No preference if no prediction occurred in this position
        pref_score[i] <- 0
      }else{ # Assign preference to predicted values
        sd[i] <- sqrt(exp(2*t*Pred_matrix_combined$x[i])*sigma^2) # + exp(2*t_wls1*Pred_matrix_combined$x[i])*(sigma_wls1)^2*
        #                (1/n + (Pred_matrix_combined$x[i]-mean(Pred_matrix_combined$x))^2/((n-1)*var(Pred_matrix_combined$x))))
        if(Pred_matrix_combined$y[i] < wls1_fun_l(Pred_matrix_combined$x[i])){
          pref_score[i] <- (wls1_fun(Pred_matrix_combined$x[i]) - Pred_matrix_combined$y[i])/sd[i]
        } else{
          if(Pred_matrix_combined$y[i] > wls1_fun_u(Pred_matrix_combined$x[i])){
            pref_score[i] <- (wls1_fun(Pred_matrix_combined$x[i]) - Pred_matrix_combined$y[i])/sd[i]
          }
        }
      }
    }
  }
  # Convert the prediction vector into prediction matrix matching each position of genomic regions
  Pred_matrix = matrix(pref_score,ncol=n_pred_pos,byrow = T)
  Pred_matrix = round(Pred_matrix,2) # Define decimal points 
  Pred_matrix = data.frame(Pred_genomic_coord,Pred_matrix) # add genomic coordinate info
  
  #Plot the two weighted regression models with prediction bands
  if(plot == T){
    plot(Pred_matrix_combined, col = rgb(0,0,0,0.1), xlab = xlab, ylab = ylab,cex=.3)
    lines(sort(data_fit$x), wls1_fit[order(data_fit$x)], lwd = 2, col = 2)
    lines(sort(data_fit$x), wls1_l[order(data_fit$x)], lty = 2, col = 2)
    lines(sort(data_fit$x), wls1_u[order(data_fit$x)], lty = 2, col = 2)
    lines(sort(data_fit$x), wls2_fit[order(data_fit$x)], lwd = 2, col = 3)
    lines(sort(data_fit$x), wls2_l[order(data_fit$x)], lty = 2, col = 3)
    lines(sort(data_fit$x), wls2_u[order(data_fit$x)], lty = 2, col = 3)
    if(res_anova$`p-value`[2] < 0.05){
      curve(wls2_fun, col = 4, add = T)
      curve(wls2_fun_l, col = 4, lty = 2, add = T)
      curve(wls2_fun_u, col = 4, lty = 2, add = T)
    } else{
      curve(wls1_fun, col = 4, add = T)
      curve(wls1_fun_l, col = 4, lty = 2, add = T)
      curve(wls1_fun_u, col = 4, lty = 2, add = T)
    }
  }
  
  return(Pred_matrix) #returns pref as a data frame
}
