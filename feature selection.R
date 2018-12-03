# First Step of viewing normalized spectra count protomics data #

# Data loading 

train<-read.csv('training.csv')
y<-read.csv('y.csv') # Y chromosome proteins


# Here we load libraries for data wrangling and visualisation.
library(tidyverse)
library(limma)
library(ggalluvial)



# peak at the dataset
glimpse(train)

# getting Y chromosome proteins for gender prediction
x <- as.character(y$gene)
y$protein <- gsub( "\\s.*", "", x )
yprotein <- intersect(y$protein,names(train)) # Y chromosome protein found in the training set


#### feature selection with limma using proteins contains no NAs ####

# remove columns contain NAs  
train_NA <- train[,complete.cases(t(train))] 

# create result and expression matrices 
dmsi<-as.data.frame(ifelse(train_NA$msi=='Lo',0,1))
x1<-t(train_NA[,4:(length(train_NA))])

# fit the model
msi_lm<-lmFit(x1, dmsi) #robust linear model

### standard errors have been moderated across proteins (shrunk towards a common value,
### using a Bayesian model)
msi_lm<-eBayes(msi_lm)

# generates top20 with lowest p value
msi_list<-as.data.frame(msi_lm,coefficient=1)
msi_list<-cbind(rownames(msi_list),msi_list)
msi_list<-arrange(msi_list,(p.value))
msi_toplist<-as.character(msi_list[1:20,1]) # top20 with lowest p value

rm(msi_lm); invisible(gc())
##### Check Protein Expression #####

# organize table
lm_top <-select(train, sample, gender, msi, msi_toplist,yprotein)
data <- gather(data = lm_top, key = Class, value = Abundance, -c(1:3))
table<-summarise(group_by(data, Class, msi, gender), expression=mean(Abundance))

table

# box plot 
x_lo<-subset(lm_top, lm_top$msi=='Lo')
x_hi<-subset(lm_top, lm_top$msi=='Hi')
xlo_f<-subset(x_lo, x_lo$gender=='f')
xlo_m<-subset(x_lo, x_lo$gender=='m')
xhi_f<-subset(x_hi, x_hi$gender=='f')
xhi_m<-subset(x_hi, x_hi$gender=='m')

par(mfrow=c(2,2))
for (i in 4:23){
   
   boxplot(x_hi[,i],xhi_f[,i],xhi_m[,i],x_lo[,i],xlo_f[,i],xlo_m[,i], 
           col = c('blue', 'red', 'steelblue','orange','cyan4','brown1'), 
           xlab='msi-H f_Hi m_Hi msi-L f_Lo m_Lo', ylab=colnames(x_hi)[i])
}

