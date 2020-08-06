rm(list = ls())

## Recall packages
library(fda)
library(fdapace)
library(graphics)
library(boot)
library(dplyr)

# Read data
Gene <- read.csv("~/SNPs_top100.csv", header = T, fill = T, stringsAsFactors = F) 
DMN4_base <- read.csv("~/DMN4_BASE_OLDSPM.csv", header = T, fill = T, stringsAsFactors = F) 
DMN4_Long <- read.csv("~/DMN4_LONG_OLDSPM.csv", header = T, fill = T, stringsAsFactors = F) 

## -------------------------------------------------
## Data Cleaning ##
# Filling out all missing Subject_ID
for (i in 1:nrow(DMN4_Long)){
  if (DMN4_Long$Subject_ID[i] == "")
  {DMN4_Long$Subject_ID[i]=DMN4_Long$Subject_ID[i-1]}
  else
  {DMN4_Long$Subject_ID[i]=DMN4_Long$Subject_ID[i]}
}

# Remove records for subject "002_S_4225" because it's not in the Gene data
DMN4_base <- DMN4_base[-(which(DMN4_base$Subject_ID == "002_S_4225")),]
DMN4_Long <- DMN4_Long[-(which(DMN4_Long$Subject_ID == "002_S_4225")),]

# Count days for each subject with the 1st measurement date as the subject's initial day = 0
DMN4date <- DMN4_Long$DMN4date <- as.Date(DMN4_Long$EXAMDATE)
DMN4_Long$days[1] = 0
DMN4_Long$initialdate[1] = DMN4date[1]
for (i in 2:nrow(DMN4_Long)){
  if (DMN4_Long$Subject_ID[i] == DMN4_Long$Subject_ID[i-1])
  {DMN4_Long$initialdate[i] = DMN4_Long$initialdate[i-1]}
  else
  {DMN4_Long$initialdate[i] = DMN4date[i]}
  
  if (DMN4_Long$Subject_ID[i] == DMN4_Long$Subject_ID[i-1])
  {DMN4_Long$days[i] = DMN4date[i] - DMN4_Long$initialdate[i]}
  else
  {DMN4_Long$days[i] = 0}
}

# Remove the 2nd measurement for subject "130_S_4660" 
DMN4_Long.new <- DMN4_Long[-which(DMN4_Long$Subject_ID == "130_S_4660" & DMN4_Long$EXAMDATE == "2013-08-23" & DMN4_Long$PCC..PCC == "-0.620794607"),]

DMN4_days = as.numeric(DMN4_Long.new$days)

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
merge_DMN4_base <- inner_join(x = DMN4_base, y = Gene, by = "Subject_ID")
merge_DMN4_Long <- merge(x = DMN4_Long.new, y = Gene, by.x = "Subject_ID", by.y = "Subject_ID")

################################################################
#####------------------------------------------------------#####
##### PACE Analysis of all edge weights on the covariates  #####
#####------------------------------------------------------#####
################################################################

## PACE Analysis for total 16 edge weights within 500 days
DMN4_edgeweights.names <- colnames(DMN4_Long.new[,4:19])
DMN4_Long500 <- DMN4_Long.new[which(DMN4_Long.new$days <= 500),]
length(unique(DMN4_Long500$Subject_ID)) #111

# PACE Analysis for total 16 edge weights
DMN4_Files500 <- vector("list", length = 16)
DMN4_fpcaobjFiles500 <- vector("list", length = 16)
names(DMN4_fpcaobjFiles500) <- DMN4_edgeweights.names

for (j in 4:19){
  DMN4_Files500[[j-3]] <- MakeFPCAInputs(DMN4_Long500$Subject_ID, DMN4_Long500$days,
                                         DMN4_Long500[,j])
  DMN4_fpcaobjFiles500[[j-3]] <- FPCA(DMN4_Files500[[j-3]]$Ly, DMN4_Files500[[j-3]]$Lt,
                                      list(plot = F, maxK = 3, methodMuCovEst = 'smooth'))
}


# Getting fitted value
DMN4_est_weights500 <- vector("list", length = 16)
names(DMN4_est_weights500) <- DMN4_edgeweights.names

for (i in 1:16){
  DMN4_est_weights500[[i]] = fitted(DMN4_fpcaobjFiles500[[i]], K = 3)
  names(DMN4_est_weights500)[[i]] = DMN4_edgeweights.names[i]
}

# Visualising the fitted trajectories
for (i in 1:length(DMN4_fpcaobjFiles500)){
  CreatePathPlot(DMN4_fpcaobjFiles500[[i]], 
                 K = 3,
                 showMean = T,
                 xlab = "Days",
                 ylab = colnames(DMN4_Long.new[,4:19])[i],
                 main = c(colnames(DMN4_Long.new[,4:19])[i]," based on GCV bandwidth"), 
                 pch = 16)
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
DMN4_Covariate <- array(NA, dim = c(111, 7, 100))
DMN4_Covariate[,1,] <- rep(1, 111)
for (j in 1:100){
  DMN4_Covariate[,2,j] <- merge_DMN4_base[,j+24]
}
DMN4_Covariate[,3,] <- DMN4_base$Age

for(i in 1:nrow(DMN4_base)){
  if (DMN4_base$PTGENDER[i] == "Female")
  {DMN4_base$Gender[i] = 2}
  else {DMN4_base$Gender[i] = 1}
}
DMN4_Covariate[,4,] <- DMN4_base$Gender

for(i in 1:nrow(DMN4_base)){
  if (DMN4_base$PTHAND[i] == "Right")
  {DMN4_base$Hand[i] = 2}
  else {DMN4_base$Hand[i] = 1}
}
DMN4_Covariate[,5,] <- DMN4_base$Hand
DMN4_Covariate[,6,] <- DMN4_base$PTEDUCAT
DMN4_Covariate[,7,] <- merge_DMN4_base$APOEe4


rownames(DMN4_Covariate) <- DMN4_base$Subject_ID
colnames(DMN4_Covariate) <- Cov.names <- c("Intercept", "SNPs", "Age", "Gender", "Handiness", "Education Length", "APOEe4")

## Compute beta
# 7 = intercept + 6 covariates
# 51 = time pts
# 16 = # of edge weights
# 101 = # of SNPs
# Create features contain all edge weights names and SNPs names

DMN4_SNPs.names <- colnames(merge_DMN4_base[25:124])

# Compute Beta(t) and get plots of estimated Beta(t)
DMN4_est_beta <- array(NA, dim = c(7, 51, 100, 16))

for (j in 1:16){
  for(i in 1:100){
    X <- as.matrix(DMN4_Covariate[,,i])
    DMN4_est_beta[,,i,j] <- solve(t(X)%*%X)%*%t(X)%*%DMN4_est_weights500[[j]]

    matplot(DMN4_fpcaobjFiles500[[j]]$workGrid, t(DMN4_est_beta[,,i,j]), 
            xlab="Days", ylab="est_beta",
            main = paste("Estimated Beta(t) for", DMN4_edgeweights.names[j], "on",
                         DMN4_SNPs.names[i]),
            type="l",cex.lab=1,cex.axis=1, xlim =  c(0,500), col = 1:7, lty = 1:7)
    
    par(xpd=TRUE)
    legend("bottomright",
           c("Intercept", DMN4_SNPs.names[i], 
             "Age", "Gender", "Handiness", "Education Length", "APOEe4"),
           col = 1:7, lty = 1:7,cex = 0.6, bty = "n")
  }
}

# Plots of Coefficients for 6 Covariates without Intercept
for (j in 1:16){
  for(i in 1:100){
    X <- as.matrix(DMN4_Covariate[,,i])
    DMN4_est_beta[,,i,j] <- solve(t(X)%*%X)%*%t(X)%*%DMN4_est_weights500[[j]]

    matplot(DMN4_fpcaobjFiles500[[j]]$workGrid, t(DMN4_est_beta[-1,,i,j]), 
            xlab="Days", ylab="est_beta",
            main = paste("Estimated Beta(t) for", DMN4_edgeweights.names[j], "on",
                         DMN4_SNPs.names[i]),
            type="l",cex.lab=1,cex.axis=1, xlim =  c(0,500), col = 1:6, lty = 1:6)
    
    par(xpd=TRUE)
    legend("bottomright",
           c(DMN4_SNPs.names[i], 
             "Age", "Gender", "Handiness", "Education Length", "APOEe4"),
           col = 1:6, lty = 1:6,cex = 0.6, bty = "n")
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
DMN4_y_hat <- array(NA, dim = c(111, 51, 100, 16))

for (j in 1:16){
  for (i in 1:100){
    X <- as.matrix(DMN4_Covariate[,,i])
    DMN4_y_hat[,,i,j] <- X%*%DMN4_est_beta[,,i,j]
  }
}

# Compute SSE1
DMN4_SSE1 <- array(NA, dim = c(100, 16))

for (j in 1:16){
  for (i in 1:100){
    DMN4_SSE1[i,j] <- sum((DMN4_y_hat[,,i,j] - DMN4_est_weights500[[j]])^2)
  }
}

# Convariance function of the residuals from  the full model (100*16 Cov functions)
# Get residuals based on est_beta and covariates
DMN4_res_hat <- array(NA, dim = c(111, 51, 100, 16))
for (j in 1:16){
  for (i in 1:100){
    DMN4_res_hat[,,i,j] <- DMN4_y_hat[,,i,j] - DMN4_est_weights500[[j]]
  }
}

# Get 1600 Cov functions (each one is a 51*51 matrix)
DMN4_timepts <- seq(1, 51)
DMN4_rcov_hat <- array(NA, dim = c(51, 51, 100, 16))
for(j in 1:16){
  for (i in 1:100){
    for (t in DMN4_timepts){
      for (s in DMN4_timepts){
        DMN4_rcov_hat[t,s,i,j] <- var(DMN4_res_hat[,t,i,j], DMN4_res_hat[,s,i,j])
      }
    }
  }
}

# Compute the d.f.(f1 = trace(E)^2/trace(E%*%E)) for the functional F distribution which the F statistic follows
DMN4_f1 <- array(NA, dim = c(100, 16))
for (j in 1:16){
  for (i in 1:100){
    E.temp <- DMN4_rcov_hat[,,i,j]
    tr.E <- sum(diag(E.temp))
    tr.EE <- sum(diag(E.temp%*%E.temp))
    DMN4_f1[i,j] <- tr.E^2/tr.EE
  }
}

## -------------------------------------------------
##### Reduced Model - M_0 (without covariate SNP)

#Covariates Matrix for the reduced model M_0
DMN4_Covariate_rm <- DMN4_Covariate[,-2,1]

## Compute beta for the reduced model
# Compute Beta(t) and get plots of estimated Beta(t)
DMN4_est_beta_rm <- array(NA, dim = c(6, 51, 16))
DMN4_X_rm <- as.matrix(DMN4_Covariate_rm)

for (j in 1:16 ){
  DMN4_est_beta_rm[,,j] <- solve(t(DMN4_X_rm)%*%DMN4_X_rm)%*%t(DMN4_X_rm)%*%DMN4_est_weights500[[j]]
  
  matplot(DMN4_fpcaobjFiles500[[j]]$workGrid, t(DMN4_est_beta_rm[,,j]), xlab="Days", ylab="est_beta",
          main = paste("Estimated Beta(t) for", DMN4_edgeweights.names[j]),
          type="l",cex.lab=1,cex.axis=1, xlim =  c(0,500), col = 1:6, lty = 1:6)
  
  legend("bottomright",
         c("Intercept", "Age", "Gender", "Handiness", "Education Length"),
         col = 1:6, lty = 1:6, cex = 0.6, bty = "n")
}

# Plots of Coefficients for 5 Covariates without Intercept
for (j in 1:16 ){
  DMN4_est_beta_rm[,,j] <- solve(t(DMN4_X_rm)%*%DMN4_X_rm)%*%t(DMN4_X_rm)%*%DMN4_est_weights500[[j]]
  
  matplot(DMN4_fpcaobjFiles500[[j]]$workGrid, t(DMN4_est_beta_rm[-1,,j]), xlab="Days", ylab="est_beta",
          main = paste("Estimated Beta(t) for", DMN4_edgeweights.names[j]),
          type="l",cex.lab=1,cex.axis=1, xlim =  c(0,500), col = 1:5, lty = 1:5)
  
  legend("bottomright",
         c("Age", "Gender", "Handiness", "Education Length"),
         col = 1:5, lty = 1:6, cex = 0.5, bty = "n")
}

## Compute SSE0 (for the Reduced Model)
# Get y_hat based on est_beta and covariates
DMN4_y_hat_rm <- array(NA, dim = c(111, 51, 16))

for (j in 1:16){
  DMN4_y_hat_rm[,,j] <- DMN4_X_rm%*%DMN4_est_beta_rm[,,j]
}

# Compute SSE0
DMN4_SSE0 <- vector("logical", length = 16)

for (j in 1:16){
  DMN4_SSE0[j] <- sum((DMN4_y_hat_rm[,,j] - DMN4_est_weights500[[j]])^2)
}

# Convariance function of the residuals from the reduced model (36 Cov functions)
# Get residuals based on est_beta and covariates
DMN4_res_hat_rm <- array(NA, dim = c(111, 51, 16))
for (j in 1:16){
  DMN4_res_hat_rm[,,j] <- DMN4_y_hat_rm[,,j] - DMN4_est_weights500[[j]]
}

DMN4_rcov_hat_rm <- array(NA, dim = c(51, 51, 16))
DMN4_rmu_hat_rm <- vector("list", length = 16)
for(j in 1:16){
  for (t in DMN4_timepts){
    for (s in DMN4_timepts){
      DMN4_rcov_hat_rm[t,s,j] <- var(DMN4_res_hat_rm[,t,j], DMN4_res_hat_rm[,s,j])
    }
  }
  DMN4_rmu_hat_rm[[j]] <- colMeans(DMN4_res_hat_rm[,,j])
}

## -------------------------------------------------
## Compute p-value, q-value of F test

# Compute F statistics
DMN4_F.stat <- DMN4_pval.ftest <- array(NA, dim = c(100, 16))

DMN4_n = 111
DMN4_df0 <- DMN4_n-6 # df of reduced model
DMN4_df1 <- DMN4_n-7 # df of full model

Ftest.fct <- function(x, y, df_fm, df_rm){
  result.ftest <- ((y-x)/(df_rm-df_fm))/(x/df_fm)
  return(result.ftest)
}

for (j in 1:16){
  for (i in 1:100){
    DMN4_F.stat[i,j] <- Ftest.fct(DMN4_SSE1[i,j], DMN4_SSE0[j], DMN4_df1, DMN4_df0)
  }
}

# Compute p-values and put them into a matrix
for (j in 1:16){
  for (i in 1:100){
    DMN4_pval.ftest[i,j] <- pf(DMN4_F.stat[i,j], DMN4_f1[i,j], DMN4_df1*DMN4_f1[i,j], lower.tail = F)
  }
}

DMN4_pval.ftest.matrix <- as.matrix(DMN4_pval.ftest)
colnames(DMN4_pval.ftest.matrix) <- DMN4_edgeweights.names
rownames(DMN4_pval.ftest.matrix) <- DMN4_SNPs.names

p_DMN4 = as.vector(DMN4_pval.ftest.matrix)

# Distribution of p-values obtained from the null F-dist for FSR for DMN6
hist(p_DMN4, xlab = "P-value", ylab = "Frequency", main = "Distribution of P-value")

# Compute q-values and put them into a matrix
q_DMN4 = p.adjust(p_DMN4, method = "fdr")
q_matrix_DMN4 = matrix(q_DMN4, nrow = 100, ncol = 16)
colnames(q_matrix_DMN4) <- DMN4_edgeweights.names
rownames(q_matrix_DMN4) <- DMN4_SNPs.names

# Create a data.frame containing the p-value and q-vlaue for all Edge Weight-SNP pairs
DMN4_pval.ftest.data <- data.frame(Edgeweights = NA, SNPs = NA, Pvalue = NA, Qvalue = NA)
for (j in 1:16){
  for (i in 1:100){
    DMN4_pval.ftest.data <- rbind(DMN4_pval.ftest.data, 
                                  data.frame(Edgeweights = colnames(DMN4_pval.ftest.matrix)[j], 
                                             SNPs = rownames(DMN4_pval.ftest.matrix)[i], 
                                             Pvalue = DMN4_pval.ftest.matrix[i,j],
                                             Qvalue = q_matrix_DMN4[i,j]))
  }
}

# Get all the pairs with p-value < 0.05 in increasing order
DMN4_export <- DMN4_pval.ftest.data[which(DMN4_pval.ftest.data$Pvalue < 0.05),]
DMN4_export[order(DMN4_export$Pvalue, decreasing = F),]

################################
#####----------------------#####
##### Bootstrap for F test #####
#####----------------------#####
################################

library(MASS)

# Simulate 10000 replicates w/ original Covariates for the top 20 Edge Weight-SNP pairs
DMN4_timepts <- seq(1, 51)
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
DMN4_edge.list <- DMN4_SNP.list <- vector("logical", length = 20)

for (i in 1:20){
  DMN4_edge.list[i] <- which(DMN4_edgeweights.names == DMN4_export[order(DMN4_export$Pvalue, decreasing = F),]$Edgeweights[i])
  DMN4_SNP.list[i] <- which(DMN4_SNPs.names == DMN4_export[order(DMN4_export$Pvalue, decreasing = F),]$SNPs[i])
}

# Set seeds
set.seed(1842)
pair.seeds <- sample(1:500, size = 20, replace = F)

DMN4_beta_rm <- DMN4_sim.F.stat.20pairs <- DMN4_sim.pval.ftest.20pairs <- vector("list", length=20)

# Simulation Start
for (i in 1:20){
  # Set Coviariate Lists
  res.sim <- edge.sim <- sim.est.beta <- sim.est.beta.rm <- sim.y.hat <- sim.y.hat.rm <- sim.res.hat <- sim.res.hat.rm <- sim.rcov.hat <- sim.rcov.hat.rm <- sim.covariate <-  sim.covariate.rm <- vector("list", length = 10000)
  
  sim.SSE1 <- sim.SSE0 <- sim.f1 <- sim.F.stat <- sim.pval.ftest <- vector("logical", length = 10000)
  
  set.seed(pair.seeds[i])
  res.seeds <- sample(x = 1:50000, size = 10000, replace = F)
  
  for (j in 1:10000){
    # simulate residuals
    set.seed(res.seeds[j])
    res.sim[[j]] <- mvrnorm(n=111, mu = DMN4_rmu_hat_rm[[DMN4_edge.list[i]]], 
                            Sigma = DMN4_rcov_hat_rm[,,DMN4_edge.list[i]], 
                            empirical = T)
    
    # Original covariates for the top 20 pairs
    sim.covariate[[j]] <- DMN4_Covariate[,,DMN4_SNP.list[i]]
    sim.covariate.rm[[j]] <- DMN4_Covariate_rm
    
    # Obs data computed using the original covariates with redeuce model (Null Hypothesis)
    edge.sim[[j]] <- matrix(NA, nrow = 111, ncol = 51)
    
    DMN4_beta_rm[[i]] <- DMN4_est_beta_rm[,,DMN4_edge.list[i]]
    
    edge.sim[[j]] <- sim.covariate.rm[[j]]%*%DMN4_beta_rm[[i]] + res.sim[[j]]
    colnames(edge.sim[[j]]) <- time.names
    #  rownames(edge.sim[[j]]) <- sub.names
    
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
    
    sim.timepts <- seq(1, 51)
    sim.rcov.hat[[j]] <- sim.rcov.hat.rm[[j]] <- matrix(NA, nrow = 51, ncol = 51)
    for (t in sim.timepts){
      for (s in sim.timepts){
        sim.rcov.hat[[j]][t,s] <- var(sim.res.hat[[j]][,t], sim.res.hat[[j]][,s])
      }
    }
    
    E.temp <- sim.rcov.hat[[j]]
    tr.E <- sum(diag(E.temp))
    tr.EE <- sum(diag(E.temp%*%E.temp))
    sim.f1[j] <- tr.E^2/tr.EE
    
    # Cov function of the residuals from the reduced model for bootstrap
    sim.res.hat.rm[[j]] <- sim.y.hat.rm[[j]] - edge.sim[[j]]
    for (t in DMN4_timepts){
      for (s in DMN4_timepts){
        sim.rcov.hat.rm[[j]][t,s] <- var(sim.res.hat.rm[[j]][,t], sim.res.hat.rm[[j]][,s])
      }
    }
    
    # Compute F statistics for simualtions
    sim.F.stat[j] <- Ftest.fct(sim.SSE1[j], sim.SSE0[j], sim_df1, sim_df0)
    # Compute p-values for simulations
    sim.pval.ftest[j] <- pf(sim.F.stat[j], sim.f1[j], sim_df1*sim.f1[j], lower.tail = F)
  }
  
  DMN4_sim.F.stat.20pairs[[i]] <- sim.F.stat
  DMN4_sim.pval.ftest.20pairs[[i]] <- sim.pval.ftest
  
  print(paste("pair", i, "done"))
}

# Get the Bootstrap p-values for the top 20 pairs
DMN4_b.value <- vector("logical", length = 20)
for (i in 1:20){
  DMN4_b.value[i] = sum(DMN4_sim.pval.ftest.20pairs[[i]] < DMN4_pval.ftest[DMN4_SNP.list[i],DMN4_edge.list[i]])/10000
}

