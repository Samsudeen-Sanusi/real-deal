#set working directory

setwd("~/Desktop/Thesis/TCGA_GBM_exp_u133a-2015-02-24")

#Load data
data = read.table("genomicMatrix.txt", sep = "\t", header = TRUE)

#Check data structure
str(data)

#View 
View(data)
class(data)
rownames(data)
colnames(data)

#Number of rows and columns
dim(data)

#[1] 12042   540

#Boxplot for preprocessed data 
library("RColorBrewer")
usr.col=brewer.pal(9, "Set1")
mycols=rep(usr.col,each=3)
boxplot(data, col=usr.col, las=3, names= colnames(data),ylim=c(4,12),cex.axis=0.6,main= "Preprocessd data")


#Assign row name
rownames(data) <- data[,1]
data[,1] <- NULL

#View 
dim(data)
View(data)
class(data)
rownames(data)
colnames(data)

#Convert to data matrix
data_matrix = data.matrix(data)

#Dimentionality
dim(data_matrix) 
View(data_matrix)
class(data_matrix)

#[1] 12042   539


#Samples
samples <- colnames(data_matrix)

#Probe IDs/Gene names
genes <- rownames(data_matrix)

#Filter matrix to 
#  i) elimantes irrelevant and uninteresting DE genes and 
#  ii) reduce no of hypothesis to be tested and uninteresting DE genes

##NOTE: I am keeping genes for which i) at least 270 out of 539 samples (i.e 50%) have an intensity of 100 or above; and ii) ratio of maximal/minimal intensity is at least 1.5.


library(genefilter)
s1<-pOverA(270/539, log2(100))
s2<-function(x) (diff(range(x, na.rm=T))>log2(1.5)) 
ss<-filterfun(s1,s2)
index<-genefilter(data_matrix, ss)
sum(index)

###NOTE: The filtering implemented above reduced the total number of genes from  12042 to 4576.


#filtered expression data matrix
filtered_data<-data_matrix[index,] 

#class and dimension
class(filtered_data)
dim(filtered_data)
View(filtered_data)

#[1] 4576  539

#Find variance 
variance = apply(filtered_data,1, var)

#Select top 100 variance value 
ordered_var = order(variance, decreasing = TRUE)[1:100]

#Filtered data with top variance value
mydata = filtered_data[ordered_var,]

#class and view
class(mydata)
View(mydata)
dim(mydata)

#Remove any missing value that might be present in the data
mydata = na.omit(mydata)

#Check quantile to decide whether to rescale or not
quantile(mydata)


#    0%       25%       50%       75%      100% 
#  2.196606  8.413866  9.636276 10.780020 14.413212 

#Descriptive statistics
des_stat <- data.frame(
  Min = apply(mydata, 2, min), # minimum
  Med = apply(mydata, 2, median), # median
  Mean = apply(mydata, 2, mean), # mean
  SD = apply(mydata, 2, sd), # Standard deviation
  Max = apply(mydata, 2, max) # Maximum
)

des_stat <- round(des_stat, 1)
head(des_stat)

#Since the means and standard deviation values are different, we need to scale the data.

#Scale data
scaled_data =scale(mydata)

#View
View(scaled_data)
class(scaled_data)
dim(scaled_data)
rownames(scaled_data)
colnames(scaled_data)


#Descriptive statistics after scalling
des_stat1 = data.frame(
  Min = apply(scaled_data, 2, min), # minimum
  Med = apply(scaled_data, 2, median), # median
  Mean = apply(scaled_data, 2, mean), # mean
  SD = apply(scaled_data, 2, sd), # Standard deviation
  Max = apply(scaled_data, 2, max) # Maximum
)

des_stat1 = round(des_stat1, 1)
head(des_stat1)

#Scaling gives a mean of 0 and standard deviation of 1

#Define hclust and distance
hclustfunc = function(x) hclust(x, method="ward.D2")
distfunc = function(x) dist(x,method="maximum")

#obtain the clusters
kk <- hclustfunc(distfunc(scaled_data))
clusters <- cutree(kk, 2)
clusters

#HC
library(gplots)
library(RColorBrewer)
pdf(file='my_heatmap.pdf', height=50, width=20)
heatmap.2(scaled_data, trace='none', scale='none',
          hclust=hclustfunc, distfun=distfunc, col=greenred(300), symbreak=F,xlab="Samples", ylab="Genes",
          margins=c(10,15), keysize=0.9, labRow=data$Gene.symbol,density.info='histogram', RowSideColors=as.character(clusters))
dev.off()





#heatmap.2(scaled_data, col=greenred(300), xlab="Samples", ylab="Genes", key = TRUE, keysize = 1.0)



#Correlation with Pearson's method 
correlation = cor(scaled_data,method="pearson")

#Class
class(correlation)
dim(correlation)
View(correlation)

#Correlation heatmap
library(gplots)
library(RColorBrewer)
heatmap(correlation, col = greenred(200), xlab="Samples", ylab="Genes", keysize = 1.0)

#OR

#Define hclust and dist function
hclfunc = function(x) hclust(x, method="ward.D2")
dtfunc = function(x) dist(x,method="maximum")

#obtain the clusters
dd <- hclfunc(dtfunc(correlation))
cls <- cutree(dd, 2)
cls

#Heatmap
library(gplots)
library(RColorBrewer)
pdf(file='sample_heatmap.pdf', height=50, width=20)
heatmap.2(correlation, trace='none', scale='none',
          hclust=hclfunc, disttfun=dtfunc, col=greenred(300), symbreak=F,xlab="Samples", ylab="Samples",
          margins=c(10,15), keysize=0.9,density.info='histogram', RowSideColors=as.character(cls))
dev.off()

#Convert correlation matrix to distance matrix 
dist_matrix = as.dist(1 - correlation)

#class
class(dist_matrix)
View(dist_matrix)

#Samples hierarchical clustering using Ward's method 
hcs = hclust(t(dist_matrix), method = "ward.D2")

#Plot dendogram
plot(hcs, cex = 0.6)

#Cut tree into 2 clusters
groups = cutree(hcs, k=2) 

#draw dendogram with red borders around the 2 clusters 
rect.hclust(hcs, k=2, border="red")

#Number of members in each cluster group
table(groups)

# grp.samp
#  1   2 
# 253 286


#Dendogram for scaled data
#Correlation with Pearson's method (Samples)
corr_samp = cor(scaled_data,method="pearson")

#Class
class(corr_samp)
dim(corr_samp)
View(corr_samp)


#Convert correlation matrix to distance matrix 
dist_matrix.samp = as.dist(1 - corr_samp)

#class
class(dist_matrix.samp)

#Samples hierarchical clustering using Ward's method 
hc.samp = hclust(dist_matrix.samp, method = "ward.D2")

#Plot dendogram
plot(hc.samp, cex = 0.6)

#Cut tree into 2 clusters
grp.samp = cutree(hc.samp, k=2) 

#draw dendogram with red borders around the 2 clusters 
rect.hclust(hc.samp, k=2, border="red")

#Number of members in each cluster group
table(grp.samp)

#
# grp.samp
#  1   2 
# 253 286 

#Define group and get group members
GroupA <- grp.samp==1
GroupB <- grp.samp==2

#Correlation with Pearson's method (Genes)
corr_genes = cor(t(scaled_data),method="pearson")

#Class
class(corr_genes)
dim(corr_genes)
View(corr_genes)


#Convert correlation matrix to distance matrix 
dist_matrix.genes = as.dist(1 - corr_genes)

#class
class(dist_matrix.genes)

#Samples hierarchical clustering using Ward's method 
hc.genes = hclust(dist_matrix.genes, method = "ward.D2")

#Plot dendogram
plot(hc.genes, cex = 0.6)

#Cut tree into 2 clusters
grp.genes = cutree(hc.genes, k=2) 

#draw dendogram with red borders around the 2 clusters 
rect.hclust(hc.genes, k=2, border="red")

#Number of members in each cluster group
table(grp.genes)

#   grp.genes
#    1   2 
#   60  40 

#Define group and get group members
Group_A <- grp.genes==1
Group_B <- grp.genes==2


#Cluster validation
#Install package
install.packages("clValid")
library("clValid")


#Transpose scaled data
t_scaled_data = t(scaled_data)

#Dimensionality and View 
dim(t_scaled_data)
View(t_scaled_data)

# Rows Columns
# 539   100

#Internal validation for samples
inter_samples <- clValid(t_scaled_data, 2:5, clMethods=c("hierarchical","kmeans"),validation="internal")

#Summary
summary(inter_samples)

f_data = t(filtered_data)
View(f_data)

inter <- clValid(f_data, 2:5, clMethods=c("hierarchical","kmeans"),validation="internal")

#Summary
summary(inter)


#Internal validation for genes
inter_genes <- clValid(scaled_data, 2:5, clMethods=c("hierarchical","kmeans"),validation="internal")

#Summary
summary(inter_genes)


#Stability validation
Stability <- clValid(scaled_data, 2:5, clMethods=c("hierarchical","kmeans"),validation="stability")



#B. K-MEANS CLUSTERING
#Determine optimal number of clusters
library("factoextra")
fviz_nbclust(corr_genes, kmeans, method = "silhouette")

#Compute K-means clustering
km <- kmeans(corr_genes, 2, nstart = 25)

#Visualize result
library("factoextra")
fviz_cluster(km, data = corr_genes, frame.type = "convex") + theme_minimal()



#####TESTING FOR DIFFERENTIALLY EXPRESSED GENES (CLASS COMPARISON)
###A. Fold Change

#Define group
GroupA <- grp.samp==1
GroupB <- grp.samp==2

#Fold change mean
FCM = apply(filtered_data[,GroupA],1,mean) / apply(filtered_data[,GroupB],1,mean)

##Top FCM and genes
head(FCM)

#Absolute fold change
AFCM = abs(FCM) 

#Fold change upregulated DE genes
FC_DEGsid = names(FCM)[abs(FCM)>1.0]

#Top 50 FC upregulated DE genes
head(FC_DEGsid, 50)

#Fold change downregulated DE genes
FC_DEGsid_1 = names(FCM)[abs(FCM)<1.0]

#Top 50 FC downregulated DE genes
head(FC_DEGsid_1, 50)



##B. T-TEST
ttestfun = function(x) t.test(x[GroupA],x[GroupB])$p.value 
p.value = apply(filtered_data, 1, ttestfun)

#head
head(p.value)

#p-value plot
hist(p.value, breaks = 30)

#order the genes using the p-values from t-test
tt_DEGsid <- names(p.value)[(p.value<0.05)]
head(tt_DEGsid, 50)

#Order the genes with adjusted p-values from t-test & absolute fold change
fcp.order<-filtered_data[order(which(p.value<0.05 & AFCM>1.0)),]


#dimensionailty and rownames
dim(fcp.order)
rownames(fcp.order)

#Filter row variance to pick the top 100 genes that are differentially expressed
DEGs <-fcp.order[1:100,]
DEGs

#The genes
rownames(DEGs)
colnames(DEGs)

#Class and view
class(DEGs)
View(DEGs)


#CLASS PREDICTION ANALYSIS

#A. FEATURE SELECTION 
#Pick top 100 differentially expressed genes 
DEGs <-fcp.order[1:100,]

#Number of rows and column
dim(DEGs)
View(DEGs)


#Change from matrix to dataframe
DEGs.df <- data.frame(DEGs)

#
View(DEGs.df)
colnames(DEGs.df)


#Name samples
#Define group
GroupA <- grp.samp==1
GroupB <- grp.samp==2

#Name group
names(DEGs.df)[GroupA] = "sample_A"
names(DEGs.df)[GroupB] = "sample_B"

#View
View(DEGs.df)
dim(DEGs.df)


#Transpose DATAFRAME
DEGs.t = as.data.frame(t(DEGs.df))


#View
View(DEGs.t)


#Assign column names
#rownames(DEGs.t) = DEGs.t[1, ] 
#DEGs.t = DEGs.t[-1, ] 
#DEGs.df <- as.data.frame(DEGs.t) 


#Add rownames and column "ID_number"
samplesID <- rownames(DEGs.t)
rownames(DEGs.t) <- NULL
DEGs.t <- cbind(samplesID,DEGs.t)

#craete sample numbering column
DEGs.t["samples_SN"] <- NA
DEGs.t$samples_SN <- 001:nrow(DEGs.t)


#rownames(norm_DEGs.t)
#rownames(norm_DEGs.t) <- NULL
#norm_DEGs.t <- cbind(samples,norm_DEGs.t)

#View
View(DEGs.t)
class(DEGs.t)

#Number of rows and column
dim(DEGs.t)
colnames(DEGs.t)

#data structure and information
str(DEGs.t)

#Count 
table(DEGs.t$samplesID)

# sample_A sample_B 
#   253       286 

#Check Sample ID class
sum(is.na(DEGs.t$samplesID))

#[1] 0 - means samplesID is not a factor

OR

#Class
class(samplesID)

#[1] "character"

#Convert to factor
DEGs.t$samplesID.fac <- factor(DEGs.t$samplesID, levels = c("sample_A", "sample_B"), labels = c("A", "B"))

#Normalize data
norm_fun <- function(x) {
  return ((x - min(x)) / (max(x) - min(x))) 
}

#Call the function 
norm_DEGs.t <- as.data.frame(lapply(DEGs.t[2:101], norm_fun))

#Confirm normalization
summary(norm_DEGs.t$RNF14)
View(norm_DEGs.t)
class(norm_DEGs.t)
str(norm_DEGs.t)
colnames(norm_DEGs.t)

#   Min.   1st Qu.  Median    Mean    3rd Qu.    Max. 
#  0.0000  0.4531   0.5846    0.5633  0.6878    1.0000 

#The result confirms normalization

#Add rowname
rownames(norm_DEGs.t) <- DEGs.t$sample_SN

#Use sample_SN as class labels
DEGs_cl <- DEGs.t[, 1]
names(DEGs_cl) <-DEGs.t$samples_SN
DEGs_cl[1:5]


#Create training and test set randomly
nrow(DEGs.t)   #Get number of rows
rand_permut <- sample(x = 1:539, size = 539)
rand_permut[1:5]

#Save
save(rand_permut, file='rand_permut.RData')

#Import data
load("rand_permut.RData")


#Number of training and testing sets
rand_sample_SN = DEGs.t[rand_permut, "samples_SN"]
length(rand_sample_SN)

#[1] 539, 2/3 for training and 1/3 for testing


#Training and Testing set
testing_sample= as.character(rand_sample_SN[1:180])
training_sample = as.character(rand_sample_SN[181:539])

#Subset data using normilized data (norm_DEGs.t)
training_data = norm_DEGs.t[training_sample, ]
testing_data = norm_DEGs.t[testing_sample, ]

#Training and testing class labels
training_class = DEGs_cl[training_sample]
testing_class = DEGs_cl[testing_sample]

#Display training  set
table(training_class)

#   sample_A sample_B 
#     158       201

#Display testing set
table(testing_class)

#   sample_A sample_B 
#      95       85 

#KNN ALGORITHM
#A. Choose K value

sqrt(nrow(training_data))

#Use an odd number near the square root of size of the training set. This gives k as 18.9473, therefore K is equal to 19

#
library(class)
knn_pred_19 <- knn(train = training_data, test = testing_data,cl = training_class, k=19)
knn_pred_19[1:3]

#table
table(knn_pred_19, testing_class)

#Evaluate model performance
library(gmodels)
CrossTable(x = testing_class, y = knn_pred_19, prop.chisq = FALSE)

#Test other values of K
#K = 17
library(class)
knn_pred_17 <- knn(train = training_data, test = testing_data,cl = training_class, k= 17)
knn_pred_17[1:5]

#table
table(knn_pred_17, testing_class)

#Evaluate model performance
library(gmodels)
CrossTable(x = testing_class, y = knn_pred_17, prop.chisq = FALSE)

#K = 15
library(class)
knn_pred_15 <- knn(train = training_data, test = testing_data,cl = training_class, k= 15)
knn_pred_15[1:5]

#table
table(knn_pred_15, testing_class)

#Evaluate model performance
library(gmodels)
CrossTable(x = testing_class, y = knn_pred_15, prop.chisq = FALSE)


#K = 13
library(class)
knn_pred_13 <- knn(train = training_data, test = testing_data,cl = training_class, k= 13)
knn_pred_13[1:5]

#table
table(knn_pred_13, testing_class)

#Evaluate model performance
library(gmodels)
CrossTable(x = testing_class, y = knn_pred_13, prop.chisq = FALSE)


#K = 11
library(class)
knn_pred_11 <- knn(train = training_data, test = testing_data,cl = training_class, k= 11)
knn_pred_11[1:5]

#table
table(knn_pred_11, testing_class)

#Evaluate model performance
library(gmodels)
CrossTable(x = testing_class, y = knn_pred_11, prop.chisq = FALSE)

#K = 9
library(class)
knn_pred_9 <- knn(train = training_data, test = testing_data,cl = training_class, k= 15)
knn_pred_9[1:5]

#table
table(knn_pred_9, testing_class)

#Evaluate model performance
library(gmodels)
CrossTable(x = testing_class, y = knn_pred_9, prop.chisq = FALSE)

#K = 7
library(class)
knn_pred_7 <- knn(train = training_data, test = testing_data,cl = training_class, k= 7)
knn_pred_7[1:5]

#table
table(knn_pred_7, testing_class)

#Evaluate model performance
library(gmodels)
CrossTable(x = testing_class, y = knn_pred_5, prop.chisq = FALSE)


#K = 5
library(class)
knn_pred_5 <- knn(train = training_data, test = testing_data,cl = training_class, k= 5)
knn_pred_5[1:5]

#table
table(knn_pred_5, testing_class)

#Evaluate model performance
library(gmodels)
CrossTable(x = testing_class, y = knn_pred_5, prop.chisq = FALSE)

#K = 3
library(class)
knn_pred_3 <- knn(train = training_data, test = testing_data,cl = training_class, k= 3)
knn_pred[1:3]

table(knn_pred_3, testing_class)

#Evaluate model performance
library(gmodels)
CrossTable(x = testing_class, y = knn_pred_3, prop.chisq = FALSE)

#Clusreing of filtered data
#Rownames
rownames(filtered_data)

#Heatmap
#Define hclust and distance
hclust_func <- function(x) hclust(x, method="complete")
dist_func <- function(x) dist(x,method="maximum")

#obtain the clusters
mm <- hclust_func(dist_func(filtered_data))
clusters_1 <- cutree(mm, 4)
clusters_1


#
library(gplots)
library(RColorBrewer)
pdf(file='my_heatmap2.pdf', height=50, width=20)
heatmap.2(filtered_data, trace='none', scale='none', 
          hclust=hclust_func, distfun=dist_func, col=greenred(300), symbreak=F,xlab="Samples", ylab="Genes",
          margins=c(10,15), keysize=0.9, labRow=data$Gene.symbol,density.info='histogram', RowSideColors=as.character(clusters_1))
dev.off()


#Survival analysis
#set working directory
setwd("~/Desktop/Thesis/TCGA_GBM_exp_u133a-2015-02-24")

#Install package
install.packages("xlsx")
library("xlsx")


#Clinical data
clinical_data = read.xlsx("clinical_data.xls", sheetIndex= 1, sep = "\t", header = TRUE)

#View 
class(clinical_data)
View(clinical_data)
dim(clinical_data)

#Row and column names
colnames(clinical_data)
rownames(clinical_data)

#Select columns for survival analysis
surv_data <- clinical_data[, c(1,31,43)]

#Dimension and View
dim(surv_data)
View(surv_data)
rownames(surv_data)


#Transpose scaled data matrix
scaled_data.t = t(scaled_data)

#Change from matrix to dataframe
scaled_data.df <- as.data.frame(scaled_data.t)

#Name rowname as SampleID
scaled_data.df <- cbind(SampleID = rownames(scaled_data.df), scaled_data.df)
rownames(scaled_data.df) <- NULL


#View
rownames(scaled_data.df)
class(scaled_data.df)
View(scaled_data.df)
dim(scaled_data.df)

#Merge dataframes 
new_data = merge(surv_data, scaled_data.df, by = c("SampleID"), all = F)

#Row and Column names
rownames(new_data)
colnames(new_data)
dim(new_data)
View(new_data)
class(new_data)

#Set row names
rownames(new_data) <- new_data[,1]
new_data[,1] <- NULL

#View, class and rowname
View(new_data)
class(new_data)
rownames(new_data)

#Transpose new data
new_data.t = as.data.frame(t(new_data))

#class and column name
class(new_data.t)
colnames(new_data.t)
View(new_data.t)


#Name samples
#Define group
GroupA <- grp.samp==1
GroupB <- grp.samp==2

#Assign groups
colnames(new_data.t)[GroupA] = "sample_A"
colnames(new_data.t)[GroupB] = "sample_B"

#View, class and rowname
View(new_data.t)
class(new_data.t)
colnames(new_data.t)

#Retranspose data
my_new_data = as.data.frame(t(new_data.t))

#View, class and rowname
View(my_new_data)
class(my_new_data)
colnames(my_new_data)
rownames(my_new_data)

#Make  SampleID as columnn and change rownames to number
my_new_data$SampleID = rownames(my_new_data)
rownames(my_new_data) = 1:nrow(my_new_data)

#view and dimension
View(my_new_data)
dim(my_new_data)
class(my_new_data)
rownames(my_new_data)
colnames(my_new_data)

#Column sum of NA rows
colSums(is.na(my_new_data))
rowSums(is.na(my_new_data))

#Remove NA rows
my_new_data1 <- my_new_data[complete.cases(my_new_data), ]
str(my_new_data1)

#view and dimension
View(my_new_data1)
dim(my_new_data1)
class(my_new_data1)
rownames(my_new_data1)

#Install Survival package
install.packages("OIsurv")
library("OIsurv")
library("survival")


####
#Create survival object using SampleID
surv_object_ID <- survfit(Surv(as.numeric(my_new_data1$X_TIME_TO_EVENT, my_new_data1$X_EVENT))~my_new_data1$SampleID)

#Get Kaplan-Meier Estimate
summary(surv_object_ID)$surv

#View survival object
summary(surv_object_ID)

#print
print(surv_object_ID)

#Plot Kaplan-Meier Estimate
plot(surv_object_ID, main="Kaplan-Meier estimate with 95% confidence bounds", xlab="time in Days", ylab="survival function", lab=c(10, 10, 7), col=c("blue","red"))
mtext("K-M survival curve for clinical data",3,line=-1,cex=0.9)

plot(surv_object_ID,main="Kaplan-Meier estimate with 95% confidence bounds", xlab="Survival Time in Days", ylab="% Survivorship", col=c("blue","red"))
legend('topright','groups',c("Sample A", "Sample B"), lty=c("solid"), col=c("blue","red"))

#Comparison survival curve
survdiff(Surv(as.numeric(my_new_data1$X_TIME_TO_EVENT,my_new_data1$X_EVENT))~my_new_data1$SampleID, data = my_new_data1, na.action = na.omit, rho = 0)













