
library(hierGWAS)

data = fread('/Users/EugMac/Downloads/NEW_SNP/plink_recode_genodata.raw',head=T) 
# names of Snps 
data<- data[,c(7:106)]
#impute missing snps_none
dat1 = apply(data, 2, as.numeric) #matrix of snps
# use the median to impute the missing snps
for ( i in which(apply(dat1,2, function(x) {any(is.na(x))}))){
  print(i)
  tempm = dat1[,i]
  tempm[is.na(tempm)]=median(tempm,na.rm=TRUE)
  dat1[,i] = tempm
           }

dat1<-as.data.frame(dat1)
data<-as.data.frame(dat1)
data=data.matrix(data, rownames.force = NA)

# PLOT DENDOGRAM
SNPindex.chrom <- seq(1,100)
chrom <- cluster.snp(data,SNP_index = SNPindex.chrom)
plot(chrom) 


# PCA to extract first principal components of Snps in Cluster.
pca=matrix(0,nrow=112,ncol=12)

G1=c("rs4738020_T", "rs10816805_T", "rs3793566_G", "rs17201409_A", "rs4609263_A")
g1=data[,c(G1)]
g1.pca <- prcomp(g1, center = TRUE,scale. = TRUE)
# extract first PCA score 
pca[,1]=g1.pca$x[,1]

G2=c("kgp6102950_A", "kgp5509939_A", "gp6435439_G", "rs17102906_C")
g2=data[,c(G2)]
g2.pca <- prcomp(g2, center = TRUE,scale. = TRUE)
# extract first PCA score 
pca[,2]=g2.pca$x[,1]

G3=c("kgp3853223_T","kgp4941056_T","rs1203112_A","rs11172895_A","kgp7184403_A","kgp1517824_T","rs16914582_G")
g3=data[,c(G3)]
g3.pca <- prcomp(g3, center = TRUE,scale. = TRUE)
# extract first PCA score 
pca[,3]=g3.pca$x[,1]

G4=c("kgp7750680_C","rs11101707_T","rs6021246_C","rs949200_T","kgp8588069_A","kgp9915984_T","rs2383376_T","kgp5624081_A","kgp1838794_A")
g4=data[,c(G4)]
g4.pca <- prcomp(g4, center = TRUE,scale. = TRUE)
# extract first PCA score 
pca[,4]=g4.pca$x[,1]

G5=c("rs11674577_T", "kgp8575927_C", "rs35061433_T", "kgp10801842_G", "kgp12129398_A")
g5=data[,c(G5)]
g5.pca <- prcomp(g5, center = TRUE,scale. = TRUE)
# extract first PCA score 
pca[,5]=g5.pca$x[,1]


G6=c("rs7935380_T","rs11601321_G","rs11664831_A","kgp3648016_A","kgp11956112_G","rs4510876_T","rs13277723_T","rs11074280_G","kgp4646865_G"); 
g6=data[,c(G6)]
g6.pca <- prcomp(g6, center = TRUE,scale. = TRUE)
# extract first PCA score 
pca[,6]=g6.pca$x[,1]

G7=c("rs2465362_A","kgp9051645_A","kgp5238984_T","kgp9055011_T","kgp3071374_C","rs9388153_T","rs2860410_C","kgp22784175_A","rs7335304_A","rs9317920_G","kgp12216228_G","rs1935110_T")
g7=data[,c(G7)];
g7.pca <- prcomp(g7, center = TRUE,scale. = TRUE)
# extract first PCA score 
pca[,7]=g7.pca$x[,1]



G8=c("kgp5498627_G","kgp1813658_T","kgp7439282_C","kgp9932182_A","rs2646852_G","rs11649752_T","rs1461688_T","rs12283068_A","rs936909_T","kgp12567235_C","rs2030791_C","rs7807857_A")
g8=data[,c(G8)]
g8.pca <- prcomp(g8, center = TRUE,scale. = TRUE)
pca[,8]=g8.pca$x[,1]

G9=c("rs8060934_C","rs1109334_C","kgp7794157_T","kgp11808627_G","rs1625700_A","kgp8158375_C","kgp1147116_A","rs7200842_G","kgp5688649_C","rs6897885_C","kgp2936399_A","kgp9166769_T","rs860876_A","rs7207116_G","kgp10004648_G")
g9=data[,c(G9)]
g9.pca <- prcomp(g9, center = TRUE,scale. = TRUE)
pca[,9]=g9.pca$x[,1]

G10=c("rs6667615_A","rs10512121_G","kgp6027866_C","kgp4931190_C","kgp9433690_G","rs13287994_A")
g10=data[,c(G10)]
g10.pca <- prcomp(g10, center = TRUE,scale. = TRUE)
pca[,10]=g10.pca$x[,1]


G11=c("rs1572421_T","kgp6660468_C","kgp11708175_T","kgp320179_A","rs4699458_T","kgp11152296_A","rs4713396_A","kgp7161203_A")
g11=data[,c(G11)]
g11.pca <- prcomp(g11, center = TRUE,scale. = TRUE)
pca[,11]=g11.pca$x[,1]

G12=c("rs1417416_A","rs12025826_G","kgp1127285_C","kgp6301628_C","rs3862175_C","rs11873190_T","rs9635857_G","kgp2283014_G")
g12=data[,c(G12)]
g12.pca <- prcomp(g12, center = TRUE,scale. = TRUE)
pca[,12]=g12.pca$x[,1]



pca_data<-as.data.frame(pca)
cluster<-paste0('C',1:12)
rownames(pca_data)=rownames(data)
colnames(pca_data)=cluster
write.csv(pca_data,'/Users/EugMac/Downloads/NEW_SNP/pca_12.csv')



  












