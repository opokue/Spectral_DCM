rm(list = ls())

## Recall packages
library(fda)
library(fdapace)
library(graphics)
library(boot)
library(dplyr)

# Read data 
Gene <- read.csv("~/SNPs_top100.csv", header = T, fill = T, stringsAsFactors = F) 
DMN6_Long <- read.csv("~/DMN6_Long.csv", header = T, fill = T, stringsAsFactors = F) 
DMN6_base <- read.csv("~/DMN6_Baseline.csv", header = T, fill = T, stringsAsFactors = F)

## -------------------------------------------------
## Data Cleaning ##
# Filling out all missing Subject_ID
for (i in 1:nrow(DMN6_Long)){
  if (DMN6_Long$Subject_ID[i] == "")
  {DMN6_Long$Subject_ID[i]=DMN6_Long$Subject_ID[i-1]}
  else
  {DMN6_Long$Subject_ID[i]=DMN6_Long$Subject_ID[i]}
}

# Remove records for subject "002_S_4225" because it's not in the Gene data
DMN6_base <- DMN6_base[-(which(DMN6_base$Subject_ID == "002_S_4225")),] 
DMN6_Long <- DMN6_Long[-(which(DMN6_Long$Subject_ID == "002_S_4225")),]

# Count days for each subject with the 1st measurement date as the subject's initial day = 0
DMN6date <- DMN6_Long$DMN6date <- as.Date(DMN6_Long$EXAMDATE)
DMN6_Long$days[1] = 0
DMN6_Long$initialdate[1] = DMN6date[1]
for (i in 2:nrow(DMN6_Long)){
  if (DMN6_Long$Subject_ID[i] == DMN6_Long$Subject_ID[i-1])
  {DMN6_Long$initialdate[i] = DMN6_Long$initialdate[i-1]}
  else
  {DMN6_Long$initialdate[i] = DMN6date[i]}
  
  if (DMN6_Long$Subject_ID[i] == DMN6_Long$Subject_ID[i-1])
  {DMN6_Long$days[i] = DMN6date[i] - DMN6_Long$initialdate[i]}
  else
  {DMN6_Long$days[i] = 0}
}

# Remove the 2nd measurement for subject "130_S_4660" (deplicate measurements on same day for this subject)
DMN6_Long.new <- DMN6_Long[-which(DMN6_Long$Subject_ID == "130_S_4660" & DMN6_Long$EXAMDATE == "2013-08-23" & DMN6_Long$PCC_PCC == "-0.46063"),]

DMN6_days = as.numeric(DMN6_Long.new$days)

## -------------------------------------------------
## Missing Data Imputation ##
# Dealing with missing values in Gene data
anyNA(Gene)
SNPs.top100 <- Gene[,c(2:102)]
Mat <- apply(Gene[,c(2:102)], 2, as.numeric)
anyNA(Mat)
Temp <- c()
# Impute the missing SNPs 
for (i in which(apply(Mat, 2, function(x) {any(is.na(x))}))){
  Temp <- Mat[,i]
  Temp[is.na(Temp)] <- median(Temp, na.rm = T)
  Mat[,i] <- Temp
}
anyNA(Mat)
Gene[,c(2:102)] <- Mat
anyNA(Gene)

## -------------------------------------------------
## Merging Datasets ##
merge_DMN6_base <- inner_join(x = DMN6_base, y = Gene, by = "Subject_ID")
merge_DMN6_Long <- merge(x = DMN6_Long.new, y = Gene, by.x = "Subject_ID", by.y = "Subject_ID")

################################################################
#####------------------------------------------------------#####
##### PACE Analysis of all edge weights on the covariates  #####
#####------------------------------------------------------#####
################################################################

## PACE Analysis for total 64 edge weights within 500 days ##
DMN6_edgeweights.names <- colnames(DMN6_Long.new[,3:38])
DMN6_Long500 <- DMN6_Long.new[which(DMN6_Long.new$days <= 500),]
length(unique(DMN6_Long500$Subject_ID)) #111

# PACE Analysis for total 64 edge weights
DMN6_Files500 <- vector("list", length = 36)
DMN6_fpcaobjFiles500 <- vector("list", length = 36)
names(DMN6_fpcaobjFiles500) <- DMN6_edgeweights.names

for (j in 3:(3+36-1)){
  DMN6_Files500[[j-2]] <- MakeFPCAInputs(DMN6_Long500$Subject_ID, DMN6_Long500$days, DMN6_Long500[,j])
  DMN6_fpcaobjFiles500[[j-2]] <- FPCA(DMN6_Files500[[j-2]]$Ly, DMN6_Files500[[j-2]]$Lt,
                                      list(plot = F, maxK = 3, methodMuCovEst = 'smooth'))
}

# Getting fitted value
DMN6_est_weights500 <- vector("list", length = 36)
names(DMN6_est_weights500) <- DMN6_edgeweights.names

for (i in 1:36){
  DMN6_est_weights500[[i]] = fitted(DMN6_fpcaobjFiles500[[i]], K = 3)
  names(DMN6_est_weights500)[[i]] = DMN6_edgeweights.names[i]
}

# Visualising the fitted trajectories (36 plots)
for (i in 1:length(DMN6_fpcaobjFiles500)){
  CreatePathPlot(DMN6_fpcaobjFiles500[[i]], 
                 K = 3,
                 showMean = T,
                 xlab = "Days",
                 ylab = colnames(DMN6_Long.new[,3:38])[i],
                 main = c(colnames(DMN6_Long.new[,3:38])[i]," based on GCV bandwidth"), 
                 pch = 20)
}

###################################################
#####-----------------------------------------#####
##### FFunctional on Scalar Regression Model  #####
#####-----------------------------------------#####
###################################################

## Covariate Matrix 
# 111 = # of subjects
# 100 = # of SNPs
# 7 = intercept + SNP + Age + Gender + Handiness + Education Length + APOEe4
DMN6_Covariate <- array(NA, dim = c(111, 7, 100))
DMN6_Covariate[,1,] <- rep(1, 111)
for (j in 1:100){
  DMN6_Covariate[,2,j] <- merge_DMN6_base[,j+44]
}
DMN6_Covariate[,3,] <- DMN6_base$Age

for(i in 1:nrow(DMN6_base)){
  if (DMN6_base$PTGENDER[i] == "Female")
  {DMN6_base$Gender[i] = 2}
  else {DMN6_base$Gender[i] = 1}
}
DMN6_Covariate[,4,] <- DMN6_base$Gender

for(i in 1:nrow(DMN6_base)){
  if (DMN6_base$PTHAND[i] == "Right")
  {DMN6_base$Hand[i] = 2}
  else {DMN6_base$Hand[i] = 1}
}
DMN6_Covariate[,5,] <- DMN6_base$Hand
DMN6_Covariate[,6,] <- DMN6_base$PTEDUCAT
DMN6_Covariate[,7,] <- merge_DMN6_base$APOEe4

rownames(DMN6_Covariate) <- DMN6_base$Subject_ID
colnames(DMN6_Covariate) <- Cov.names <- c("Intercept", "SNPs", "Age", "Gender", "Handiness", "Education Length", "APOEe4")

## Compute beta
# 7 = intercept + 6 covariates
# 51 = time pts
# 36 = # of edge weights
# 100 = # of SNPs
# Create features contain all edge weights names and SNPs names
DMN6_SNPs.names <- colnames(merge_DMN6_base[45:144])

# Compute Beta(t) and get plots of estimated Beta(t)
DMN6_est_beta <- array(NA, dim = c(7, 51, 100, 36))

for (j in 1:36){
  for(i in 1:100){
    X <- as.matrix(DMN6_Covariate[,,i])
    DMN6_est_beta[,,i,j] <- solve(t(X)%*%X)%*%t(X)%*%DMN6_est_weights500[[j]]

    matplot(DMN6_fpcaobjFiles500[[j]]$workGrid, t(DMN6_est_beta[,,i,j]), 
            xlab="Days", ylab="est_beta",
            main = paste("Estimated Beta(t) for", DMN6_edgeweights.names[j], "on",
                         DMN6_SNPs.names[i]),
            type="l",cex.lab=1,cex.axis=1, xlim =  c(0,500), col = 1:7, lty = 1:7)
    
    par(xpd=TRUE)
    legend("bottomright",
           c("Intercept", DMN6_SNPs.names[i], 
             "Age", "Gender", "Handiness", "Education Length", "APOEe4"),
           col = 1:7, lty = 1:7,cex = 0.6, bty = "n")
  }
}

# Plots of Coefficients for 6 Covariates without Intercept
for (j in 1:36){
  for(i in 1:100){
    X <- as.matrix(DMN6_Covariate[,,i])
    DMN6_est_beta[,,i,j] <- solve(t(X)%*%X)%*%t(X)%*%DMN6_est_weights500[[j]]

    matplot(DMN6_fpcaobjFiles500[[j]]$workGrid, t(DMN6_est_beta[-1,,i,j]), 
            xlab="Days", ylab="est_beta",
            main = paste("Estimated Beta(t) for", DMN6_edgeweights.names[j], "on",
                         DMN6_SNPs.names[i]),
            type="l",cex.lab=1,cex.axis=1, xlim =  c(0,500), col = 1:7, lty = 1:7)
    
    par(xpd=TRUE)
    legend("bottomright",
           c(DMN6_SNPs.names[i], 
             "Age", "Gender", "Handiness", "Education Length", "APOEe4"),
           col = 1:7, lty = 1:7,cex = 0.6, bty = "n")
  }
}

##########################
#####----------------#####
##### F-test for FSR #####
#####----------------#####
##########################

## -------------------------------------------------
##### Full Model - M_1

## Compute SSE1 (for the Full Model)
# Get y_hat based on est_beta and covariates
DMN6_y_hat <- array(NA, dim = c(111, 51, 100, 36))
for (j in 1:36){
  for (i in 1:100){
    X <- as.matrix(DMN6_Covariate[,,i])
    DMN6_y_hat[,,i,j] <- X%*%DMN6_est_beta[,,i,j]
  }
}

# Compute SSE1
DMN6_SSE1 <- array(NA, dim = c(100, 36))
for (j in 1:36){
  for (i in 1:100){
    DMN6_SSE1[i,j] <- sum((DMN6_y_hat[,,i,j] - DMN6_est_weights500[[j]])^2)
  }
}

# Convariance function of the residuals from  the full model (100*36 Cov functions)
# Get residuals based on est_beta and covariates
DMN6_res_hat <- array(NA, dim = c(111, 51, 100, 36))
for (j in 1:36){
  for (i in 1:100){
    DMN6_res_hat[,,i,j] <- DMN6_y_hat[,,i,j] - DMN6_est_weights500[[j]]
  }
}

# Get 3600 Cov functions (each one is a 51*51 matrix)
DMN6_timepts <- seq(1, 51)
DMN6_rcov_hat <- array(NA, dim = c(51, 51, 100, 36))
for(j in 1:36){
  for (i in 1:100){
    for (t in DMN6_timepts){
      for (s in DMN6_timepts){
        DMN6_rcov_hat[t,s,i,j] <- var(DMN6_res_hat[,t,i,j], DMN6_res_hat[,s,i,j])
      }
    }
  }
}

# Compute the d.f.(f1 = trace(E)^2/trace(E%*%E)) for the functional F distribution which the F statistic follows
DMN6_f1 <- array(NA, dim = c(100, 36))
for (j in 1:36){
  for (i in 1:100){
    E.temp <- DMN6_rcov_hat[,,i,j]
    tr.E <- sum(diag(E.temp))
    tr.EE <- sum(diag(E.temp%*%E.temp))
    DMN6_f1[i,j] <- tr.E^2/tr.EE
  }
}

## -------------------------------------------------
##### Reduced Model - M_0 (without covariate SNP)

#Covariates Matrix for the reduced model M_0
DMN6_Covariate_rm <- DMN6_Covariate[,-2,1]

## Compute beta for the reduced model
# Compute Beta(t) and get plots of estimated Beta(t)
DMN6_est_beta_rm <- array(NA, dim = c(6, 51, 36))
DMN6_X_rm <- as.matrix(DMN6_Covariate_rm)

for (j in 1:36){
  DMN6_est_beta_rm[,,j] <- solve(t(DMN6_X_rm)%*%DMN6_X_rm)%*%t(DMN6_X_rm)%*%DMN6_est_weights500[[j]]
  
  matplot(DMN6_fpcaobjFiles500[[j]]$workGrid, t(DMN6_est_beta_rm[,,j]), xlab="Days", ylab="est_beta",
          main = paste("Estimated Beta(t) for", DMN6_edgeweights.names[j]),
          type="l",cex.lab=1,cex.axis=1, xlim = c(0,500))
  
  par(xpd = T)
  legend("bottomright",
         c("Intercept", "Age", "Gender", "Handiness", "Education Length", "APOEe4"),
         col = 1:6, lty = 1:6,cex = 0.6, bty = "n")
}


# Plots of Coefficients for 5 Covariates without Intercept
for (j in 1:36){
  DMN6_est_beta_rm[,,j] <- solve(t(DMN6_X_rm)%*%DMN6_X_rm)%*%t(DMN6_X_rm)%*%DMN6_est_weights500[[j]]
  
  matplot(DMN6_fpcaobjFiles500[[j]]$workGrid, t(DMN6_est_beta_rm[-1,,j]), xlab="Days", ylab="est_beta",
          main = paste("Estimated Beta(t) for", DMN6_edgeweights.names[j]),
          type="l",cex.lab=1,cex.axis=1, xlim = c(0,500))
  
  par(xpd = T)
  legend("bottomright",
         c("Age", "Gender", "Handiness", "Education Length", "APOEe4"),
         col = 1:5, lty = 1:5,cex = 0.6, bty = "n")
}

## Compute SSE0 (for the Reduced Model)
# Get y_hat based on est_beta and covariates
DMN6_y_hat_rm <- array(NA, dim = c(111, 51, 36))

for (j in 1:36){
  DMN6_y_hat_rm[,,j] <- DMN6_X_rm%*%DMN6_est_beta_rm[,,j]
}

# Compute SSE0
DMN6_SSE0 <- vector("logical", length = 36)

for (j in 1:36){
  DMN6_SSE0[j] <- sum((DMN6_y_hat_rm[,,j] - DMN6_est_weights500[[j]])^2)
}

# Convariance function of the residuals from the reduced model (36 Cov functions)
# Get residuals based on est_beta and covariates
DMN6_res_hat_rm <- array(NA, dim = c(111, 51, 36))
for (j in 1:36){
  DMN6_res_hat_rm[,,j] <- DMN6_y_hat_rm[,,j] - DMN6_est_weights500[[j]]
}

DMN6_rcov_hat_rm <- array(NA, dim = c(51, 51, 36))
DMN6_rmu_hat_rm <- vector("list", length = 36)
for(j in 1:36){
  for (t in DMN6_timepts){
    for (s in DMN6_timepts){
      DMN6_rcov_hat_rm[t,s,j] <- var(DMN6_res_hat_rm[,t,j], DMN6_res_hat_rm[,s,j])
    }
  }
  DMN6_rmu_hat_rm[[j]] <- colMeans(DMN6_res_hat_rm[,,j])
}

## -------------------------------------------------
## Compute p-value, q-value of F test

# Compute F statistics
DMN6_F.stat <- DMN6_pval.ftest <- array(NA, dim = c(100, 36))

DMN6_n = 111
DMN6_df0 <- DMN6_n-6 # df of reduced model
DMN6_df1 <- DMN6_n-7 # df of full model

Ftest.fct <- function(x, y, df_fm, df_rm){
  result.ftest <- ((y-x)/(df_rm - df_fm))/(x/df_fm)
  return(result.ftest)
}

for (j in 1:36){
  for (i in 1:100){
    DMN6_F.stat[i,j] <- Ftest.fct(DMN6_SSE1[i,j], DMN6_SSE0[j], DMN6_df1, DMN6_df0)
  }
}

# Compute p-values and put them into a matrix
for (j in 1:36){
  for (i in 1:100){
    DMN6_pval.ftest[i,j] <- pf(DMN6_F.stat[i,j], DMN6_f1[i,j], DMN6_df1*DMN6_f1[i,j], lower.tail = F)
  }
}

DMN6_pval.ftest.matrix <- as.matrix(DMN6_pval.ftest)
colnames(DMN6_pval.ftest.matrix) <- DMN6_edgeweights.names
rownames(DMN6_pval.ftest.matrix) <- DMN6_SNPs.names

p_DMN6 = as.vector(DMN6_pval.ftest.matrix)

# Distribution of p-values obtained from the null F-dist for FSR for DMN6
hist(p_DMN6, xlab = "P-value", ylab = "Frequency", main = "Distribution of P-value")

# Compute q-values and put them into a matrix
q_DMN6 = p.adjust(p_DMN6, method = "fdr")
q_matrix_DMN6 = matrix(q_DMN6, nrow = 100, ncol = 36)
colnames(q_matrix_DMN6) <- DMN6_edgeweights.names
rownames(q_matrix_DMN6) <- DMN6_SNPs.names

# Create a data.frame containing the p-value and q-vlaue for all Edge Weight-SNP pairs
DMN6_pval.ftest.data <- data.frame(Edgeweights = NA, SNPs = NA, Pvalue = NA, Qvalue = NA)
for (j in 1:36){
  for (i in 1:100){
    DMN6_pval.ftest.data <- rbind(DMN6_pval.ftest.data, 
                                  data.frame(Edgeweights = colnames(DMN6_pval.ftest.matrix)[j], 
                                             SNPs = rownames(DMN6_pval.ftest.matrix)[i], 
                                             Pvalue = DMN6_pval.ftest.matrix[i,j],
                                             Qvalue = q_matrix_DMN6[i,j]))
  }
}

# Get all the pairs with p-value < 0.05 in increasing order
DMN6_export <- DMN6_pval.ftest.data[which(DMN6_pval.ftest.data$Pvalue < 0.05),]
DMN6_export[order(DMN6_export$Pvalue, decreasing = F),]

################################
#####----------------------#####
##### Bootstrap for F test #####
#####----------------------#####
################################

library(MASS)

# Simulate 10000 replicates w/ original Covariates for the top 20 Edge Weight-SNP pairs
DMN6_timepts <- sim_timepts <- seq(1, 51)
sim_n = 111
sim_df0 <- sim_n-6 # df of reduced model
sim_df1 <- sim_n-7 # df of full model

# Create Names for time, simulations and subjects
time.names <- vector("logical", length = 51)
for (i in 1:51){
  time.names[i] <- noquote(paste("time", i, sep = ""))
}

sub.names <- vector("logical", length = 111)
for (j in 1:111){
  sub.names[j] <- paste("subject", j, sep = "")
}

# Set top 20 Edgeweights-SNP pairs list
DMN6_edge.list <- DMN6_SNP.list <- vector("logical", length = 20)

for (i in 1:20){
  DMN6_edge.list[i] <- which(DMN6_edgeweights.names == DMN6_export[order(DMN6_export$Pvalue, decreasing = F),]$Edgeweights[i])
  DMN6_SNP.list[i] <- which(DMN6_SNPs.names == DMN6_export[order(DMN6_export$Pvalue, decreasing = F),]$SNPs[i])
}

# Set seeds
set.seed(12312)
pair.seeds <- sample(1:500, size = 20, replace = F)

DMN6_beta_rm <- DMN6_sim.F.stat.20pairs <- DMN6_sim.pval.ftest.20pairs <- vector("list", length=20)

# Simulation Start
for (i in 1:20){
  # Set Coviariate Lists
  res.sim <- edge.sim <- sim.est.beta <- sim.est.beta.rm <- sim.y.hat <- sim.y.hat.rm <-  sim.res.hat <- sim.res.hat.rm <- sim.rcov.hat <- sim.rcov.hat.rm <- sim.covariate  <-  sim.covariate.rm <- vector("list", length = 10000)
  
  sim.SSE1 <- sim.SSE0 <- sim.f1 <- sim.F.stat <- sim.pval.ftest <- vector("logical", length = 10000)
  
  set.seed(pair.seeds[i])
  res.seeds <- sample(x = 1:50000, size = 10000, replace = F)
  
  for (j in 1:10000){
    # simulate residuals
    set.seed(res.seeds[j])
    res.sim[[j]] <- mvrnorm(n=111, mu = DMN6_rmu_hat_rm[[DMN6_edge.list[i]]], 
                            Sigma = DMN6_rcov_hat_rm[,,DMN6_edge.list[i]], 
                            empirical = T)

    # Original covariates for the top 20 pairs
    sim.covariate[[j]] <- DMN6_Covariate[,,DMN6_SNP.list[i]]
    sim.covariate.rm[[j]] <- DMN6_Covariate_rm
    
    # Obs data computed using the original covariates with redeuce model (Null Hypothesis)
    edge.sim[[j]] <- matrix(NA, nrow = 111, ncol = 51)
    
    DMN6_beta_rm[[i]] <- DMN6_est_beta_rm[,,DMN6_edge.list[i]]
    
    edge.sim[[j]] <- sim.covariate.rm[[j]]%*%DMN6_beta_rm[[i]] + res.sim[[j]]
    colnames(edge.sim[[j]]) <- time.names
    # rownames(edge.sim[[j]]) <- sub.names
    
    # Re-compute beta for both full & reduced model
    sim.est.beta[[j]] <- solve(t(sim.covariate[[j]])%*%sim.covariate[[j]])%*%
      t(sim.covariate[[j]])%*%edge.sim[[j]] 
    sim.est.beta.rm[[j]] <- solve(t(sim.covariate.rm[[j]])%*%sim.covariate.rm[[j]])%*%
      t(sim.covariate.rm[[j]])%*%edge.sim[[j]] 
    
    # Get y_hat based on est_beta and covariates
    sim.y.hat[[j]] <- sim.covariate[[j]]%*%sim.est.beta[[j]]
    sim.y.hat.rm[[j]] <- sim.covariate.rm[[j]]%*%sim.est.beta.rm[[j]]
    
    # Compute SSE1
    sim.SSE1[j] <- sum((sim.y.hat[[j]] - edge.sim[[j]])^2)
    sim.SSE0[j]<- sum((sim.y.hat.rm[[j]] - edge.sim[[j]])^2)
    
    # convariance function of the residuals from full model (100*36 cov functions)
    sim.res.hat[[j]] <- sim.y.hat[[j]] - edge.sim[[j]]
    
    sim.rcov.hat[[j]] <- sim.rcov.hat.rm[[j]] <- matrix(NA, nrow = 51, ncol = 51)
    for (t in sim_timepts){
      for (s in sim_timepts){
        sim.rcov.hat[[j]][t,s] <- var(sim.res.hat[[j]][,t], sim.res.hat[[j]][,s])
      }
    }
    
    E.temp <- sim.rcov.hat[[j]]
    tr.E <- sum(diag(E.temp))
    tr.EE <- sum(diag(E.temp%*%E.temp))
    sim.f1[j] <- tr.E^2/tr.EE
    
    # Cov function of the residuals from the reduced model for bootstrap
    sim.res.hat.rm[[j]] <- sim.y.hat.rm[[j]] - edge.sim[[j]]
    for (t in DMN6_timepts){
      for (s in DMN6_timepts){
        sim.rcov.hat.rm[[j]][t,s] <- var(sim.res.hat.rm[[j]][,t], sim.res.hat.rm[[j]][,s])
      }
    }
    
    # Compute F statistics for simualtions
    sim.F.stat[j] <- Ftest.fct(sim.SSE1[j], sim.SSE0[j], sim_df1, sim_df0)
    # Compute p-values for simulations
    sim.pval.ftest[j] <- pf(sim.F.stat[j], sim.f1[j], sim_df1*sim.f1[j], lower.tail = F)
  }
  
  DMN6_sim.F.stat.20pairs[[i]] <- sim.F.stat
  DMN6_sim.pval.ftest.20pairs[[i]] <- sim.pval.ftest

  print(paste("pair", i, "done"))
}

# Get the Bootstrap p-values for the top 20 pairs
DMN6_b.value <- vector("logical", length = 20)
for (i in 1:20){
  DMN6_b.value[i] = sum(DMN6_sim.pval.ftest.20pairs[[i]] < DMN6_pval.ftest[DMN6_SNP.list[i], DMN6_edge.list[i]])/10000
}
