# Parametric Bootsrap 

library(pbkrtest)
library(lme4)
# install.packages("/Users/EugMac/Downloads/pbnm_0.3.0.9003.tar.gz", repos=NULL,type="source")


############################## LME and Parametric-Bootstrap for DMN4 #####################################

# read PACE rs-fmri connectivity data

data<-read.csv('/Users/EugMac/Downloads/NEW_SNP/DMN4_PACE_OLD_SPM.csv')

# change network edges names
b=paste0('y',1:16)
names(data)=c("PTID","Age",b,"EXAMDATE","DXCURREN_from_DXCHANGE","PTGENDER","PTHAND","PTEDUCAT")

# read genetic data
data1<-read.csv("/Users/EugMac/Downloads/NEW_SNP/SNPs_top100.csv")
names(data1)
# SNP names
gen_data=data1[,c(3:ncol(data1))]

# change names of SNPs
sp=paste0('x',1:100)
names(data1)= c("PTID","z",sp)

# merge brain connectivity and genetic data
data$Time  <- 1:nrow(data)
dat2 = merge(data,data1,by="PTID",all=F) 
dat2<-dat2[order(dat2$Time), ]

myvars <- names(dat2) %in% c("Time", "EXAMDATE", "DXCURREN_from_DXCHANGE") 
dat2 <- dat2[!myvars]

# Parametric Bootsrap  

Y<-paste0('y',1:16)
X<-paste0('x',1:100)

pval<-matrix(0,nrow=16,ncol=100) # store p-values

for(i in 1:16){
  y<-dat2[,Y[i]]
  for(j in 1:100){
    x<-dat2[,X[j]]
    # model with snp   
    fit_full<-lmer(y~x+z+PTGENDER+PTHAND+PTEDUCAT+Age+ (1|PTID),data=dat2,REML=FALSE)
    # model without snp  
    fit_null<-lmer(y~z+PTGENDER+PTHAND+PTEDUCAT+Age+ (1|PTID),data=dat2,REML=FALSE)
    an<-PBmodcomp(fit_full,fit_null,nsim=10000,cl=1)
    # extract p-value
    pp<-an$test[2,3]
    pval[i,j]=pp
  }
}


# Linear mixed effect model 
for(i in 1:16){
  y<-dat2[,Y[i]]
  for(j in 1:100){
    x<-dat2[,X[j]]
    fit_full<-lme(y~x+z+PTGENDER+PTHAND+PTEDUCAT+Age,random=~+1|PTID,data=dat2,method="ML")
    
    fit_null<-lme(y~z+PTGENDER+PTHAND+PTEDUCAT+Age,random=~+1|PTID,data=dat2,method="ML")
    an<-anova(fit_null,fit_full)
    #cat('Y:',Y[i],'X:',X[j])
    #print(an)
    pp<-an$`p-value`
    pval[i,j]=pp[2]
  }
}

# Histogram plot of p-values
hist(pval,xlab="P-value",main="Distribution of P-values")

# ADJUST P-VALUES
p.vec=as.vector(pval)
p=p.adjust(p.vec, method = "fdr")
fdr_pvalues=matrix(p,nrow=16,ncol=100,byrow=F)

# network edge names
qq<-c("PCC->PCC","MPFC->PCC","LIPC->PCC","RIPC->PCC","PCC->MPFC","MPFC->MPFC","LIPC->MPFC",
"RIPC->MPFC","PCC->LIPC","MPFC->LIPC","LIPC->LIPC","RIPC->LIPC","PCC->RIPC",
"MPFC->RIPC","LIPC->RIPC","RIPC->RIPC")

rownames(fdr_pvalues)<-qq
colnames(fdr_pvalues)<-colnames(gen_data)

# Turn into a 3-column table
d=as.data.frame(as.table(fdr_pvalues))  
d=d[order(abs(d$Freq)),]  # Sort by adjusted p-values (whether +ve or -ve)
#d=subset(d, abs(Freq) < 0.1) # extract adjusted p-values less than .005
m =d%>%top_n(-10,wt=Freq)



########################### LME and Parametric Bootstrap for DMN6 ########################################

data<-read.csv('/Users/EugMac/Downloads/NEW_SNP/DMN6_PACE.csv')

# change network edges names 
b=paste0('y',1:36)
names(data)=c("PTID","Age",b,"EXAMDATE","DXCURREN_from_DXCHANGE","PTGENDER","PTHAND","PTEDUCAT")

# read genetic data
data1<-read.csv("/Users/EugMac/Downloads/NEW_SNP/SNPs_top100.csv")
# SNP names
gen_data=data1[,c(3:ncol(data1))]

# change names of  SNPs

sp=paste0('x',1:100)
names(data1)= c("PTID","z",sp)


# merge brain connectivity and genetic data
data$Time  <- 1:nrow(data)
dat2 = merge(data,data1,by="PTID",all=F) 
dat2<-dat2[order(dat2$Time), ]

myvars <- names(dat2) %in% c("Time", "EXAMDATE", "DXCURREN_from_DXCHANGE") 
dat2 <- dat2[!myvars]


# Fit linear mixed model 

Y<-paste0('y',1:36)
X<-paste0('x',1:100)

pval<-matrix(0,nrow=36,ncol=100) # store p-values

for(i in 1:36){
  y<-dat2[,Y[i]]
  for(j in 1:100){
    x<-dat2[,X[j]]
    # model with snp   
    fit_fulle<-lmer(y~x+z+PTGENDER+PTHAND+PTEDUCAT+Age+ (1|PTID),data=dat2,REML=FALSE)
    # model without snp  
    fit_nulle<-lmer(y~z+PTGENDER+PTHAND+PTEDUCAT+Age+ (1|PTID),data=dat2,REML=FALSE)
    ane<-PBmodcomp(fit_fulle,fit_nulle,nsim=10000,cl=1)
    # extract p-value
    ppe<-ane$test[2,3]
    pvale[i,j]=ppe
  }
}



# Linear mixed effect model 
for(i in 1:36){
  y<-dat2[,Y[i]]
  for(j in 1:100){
    x<-dat2[,X[j]]
    fit_full<-lme(y~x+z+PTGENDER+PTHAND+PTEDUCAT+Age,random=~+1|PTID,data=dat2,method="ML")
    
    fit_null<-lme(y~z+PTGENDER+PTHAND+PTEDUCAT+Age,random=~+1|PTID,data=dat2,method="ML")
    an<-anova(fit_null,fit_full)
    #cat('Y:',Y[i],'X:',X[j])
    #print(an)
    pp<-an$`p-value`
    pval[i,j]=pp[2]
  }
}

# Histogram for p-values
hist(pval,xlab="P-value",main="Distribution of P-values")



