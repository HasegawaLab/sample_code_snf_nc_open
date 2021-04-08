#############################################################
#  Example code using simulation data in the SNFtool package (revised on Jan/18/2021)
#############################################################

#Please download SNFtool, amap package for running this example
#For detaied information, please access to the original paper and SNFtool reference
#Wang, B. et al. Similarity network fusion for aggregating data types on a genomic scale. Nat. Methods 11, 333–337 (2014)
#https://cran.r-project.org/web/packages/SNFtool/SNFtool.pdf

if (!require('SNFtool')) install.packages('SNFtool'); library('SNFtool')
if (!require('amap')) install.packages('amap'); library('amap')
library(SNFtool)
library(amap)

#Set all the parameters:
K = 25; #Number of neighbors
alpha = 0.7; #Hyperparameter
t = 25; #Number of Iterations, usually (10~50)

data(Data1)
data(Data2)

#Compute distance matrices
eucli_d1= amap::Dist(Data1,method = "euclidean")%>%as.matrix()
pearson_d2= amap::Dist(Data2,method = "pearson")%>%as.matrix()

W1 = affinityMatrix(eucli_d1, K, alpha)
W2 = affinityMatrix(pearson_d2, K, alpha)

W = SNF(list(W1,W2), K, t)

C =4 # the number of clusters
group_SNF = spectralClustering(W,C)  #Final endotypes information


########################################################################################
# R code to obtain the final endotypes using similarity network fusion in the manuscript
########################################################################################

library(SNFtool)
K=25
alpha=0.7
t=25

#Compute distance matrices
gower_distance <- StatMatch::gower.dist(Data_clinical)%>%as.matrix()
bray_curtis<-vegan::vegdist(Data_microbiome, method="bray")%>%as.matrix()
pearson_transcriptome<-amap::Dist(Data_transcriptome, "pearson")%>%as.matrix()
eucli_metabolo= amap::Dist(Data_metabolome, "euclidean")%>%as.matrix()
#Construct similarity matrices
W1 = affinityMatrix(gower_distance, K, alpha)
W2 = affinityMatrix(bray_curtis, K, alpha)
W3 = affinityMatrix(pearson_transcriptome, K, alpha)
W4 = affinityMatrix(eucli_metabolo, K, alpha)

#These similarity graphs have complementary information about clusters.
#next, we fused all the similarity matrices
#Then the overall fused matrix can be computed by similarity network fusion (SNF):

W = SNF(list(W1,W2,W3,W4), K, t)

C =4 #The number of clusters
group_SNF = spectralClustering(W,C)  #Final endotypes information
