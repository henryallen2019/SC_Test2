## Classification ##

library(caret)
library(e1071)
library(rcpp)
normalize <- function(x, new_min, new_max){
  return ((( x - min(x) ) / ( max(x) - min(x) ))*(new_max - new_min)+  new_min)
}

composition <- function(param_matrix){
  result <- matrix(nrow=nrow(param_matrix), ncol=ncol(param_matrix))
  for(i in c(1:nrow(param_matrix))){
    for(j in c(1:ncol(param_matrix))){
      result[i,j] <- max(pmin(param_matrix[i,], param_matrix[,j]))
    }
  }
  return(result)
}

isEquivalence <- function(param_matrix){
  flag <- TRUE
  for(i in c(1:nrow(param_matrix))){
    if(param_matrix[i,i] != 1){
      flag <- FALSE
    }
  }
  for(i in c(1:nrow(param_matrix))){
    for(j in c(1:ncol(param_matrix))){
      if(param_matrix[i,j] != param_matrix[j,i]){
        flag <- FALSE
      }
    }
  }
  combos <- combn(nrow(param_matrix), 3)
  for(i in c(1:ncol(combos))){
    l <- param_matrix[combos[1,i], combos[2,i]]
    r <- param_matrix[combos[2,i], combos[3,i]]
    mid <- param_matrix[combos[1,i], combos[3,i]]
    if(mid < min(l, r)){
      flag <- FALSE
    }
  }
  return(flag)
}

data <- read.csv("~/iris.csv")
for(i in 1:(ncol(data)-1)){
  data[,i] <- normalize(data[,i],0,1)
}
print(data)
table <- data[1:(ncol(data)-1)]
similarityMatrix <- matrix(c(0), nrow = nrow(table), ncol = nrow(table))
for(i in c(1:nrow(table))){
  for(j in c(1:nrow(table))){
    x <- table[i,]
    y <- table[j,]
    min_res <- c()
    max_res <- c()
    for(f in c(1:length(x))){
      min_res <- c(min_res, min(x[[f]],y[[f]]))
      max_res <- c(max_res, max(x[[f]],y[[f]]))
    }
    similarityMatrix[i, j] <- sum(table[i,] * table[j,]) / (sqrt( sum(table[i,]^2) * sum(table[j,]^2))) #for cosine
    #similarityMatrix[i,j] <- sum(min_res) / sum(max_res)
    similarityMatrix[j,i] <- similarityMatrix[i,j]
  }
}
finished <- FALSE
groups <- list()
output <- rep(0, nrow(data))
while(!finished){
  if( !isEquivalence(similarityMatrix)){
    similarityMatrix <- composition(similarityMatrix)
    print(similarityMatrix)
  }
  else{
    finished <- TRUE
    cut <- 0.98
    similarityMatrix[similarityMatrix < cut] <- 0
    similarityMatrix[similarityMatrix >= cut] <- 1
    for(i in c(1:ncol(similarityMatrix))){
      position <- which(similarityMatrix[,i] == 1)
      if(length(position) > 0)
        groups[[length(groups)+1]] <- position
      label <- length(groups)
      for(j in c(1:length(position))){
        output[position[j]] <- label
      }
      similarityMatrix[position,] <- 0
    }
    break
  }
}
print(groups)

x<- vector()
labels <- vector()
for(i in c(1:length(groups))){
  x<-append(x,length(groups[[i]]))
  labels<-append(labels, paste("group ", i))
}
expectedOutput <- as.numeric(data[, ncol(data)])
y1 <- factor(expectedOutput)
y2 <- factor(output)
cm = confusionMatrix(y1,y2)
print(cm)
par(mfrow=c(1,2))
pie(x,labels,main="Count of records in each Class",col=rainbow(length(x)))
plot(data[,1],data[,2],xlab="Sepal Length",ylab="Sepal Width",col=output,main="Sepal Length vs SepalWidth",pch=2)
legend("topleft",
       legend = labels, col=c("black","red","green"),pch=2)


## Mamdani Model ##

library(frbs)
library(caret)
#use if needed
normalize <- function(x,new_min,new_max) {
  return (((x - min(x)) / (max(x) - min(x)))*(new_max-new_min)+new_min)
}
varinp.mf <- matrix(c(1, 0, 0, 30, NA, 1, 30, 45, 60, NA, 1, 60, 100, 100,NA,
                      1, 0, 0, 30, NA, 1, 30, 45, 60, NA, 1, 60, 100, 100,NA),                  
                    nrow = 5, byrow = FALSE)

num.fvalinput <- matrix(c(3, 3), nrow=1)
varinput.1 <- c("small", "medium_dirtiness", "large")
varinput.2 <- c("not_greasy", "medium_type","greasy")
names.varinput <- c(varinput.1, varinput.2)
range.data <- matrix(c(0,100, 0, 100,0,100), nrow=2)
type.defuz <- "COG"
type.tnorm <- "MIN"
type.snorm <- "MAX"
type.implication.func <- "ZADEH"
name <- "Sim-0"
#newdata<- matrix(c(30, 75,  45, 80, 30, 45), nrow= 3, byrow = TRUE)
newdata<-read.csv("~/dirt.csv")
colnames.var <- c("dirtiness", "type of dirt", "washing time")
num.fvaloutput <- matrix(c(5), nrow=1)
varoutput.1 <- c("very_short", "short", "medium","long","very_long")
names.varoutput <- c(varoutput.1)
varout.mf <- matrix(c(1, 0, 10, 20,NA ,
                      1, 20, 30, 40, NA,
                      1, 40, 50, 60, NA,
                      1,60,70,80,NA,
                      1,80,90,100,NA),
                    nrow = 5, byrow = FALSE)
type.model <- "MAMDANI"
rule <- matrix(c("large","and","greasy","->","very_long",
                 "medium_dirtiness","and","greasy","->","long",
                 "small","and","greasy","->","long",
                 "large","and","medium_type","->","long",
                 "medium_dirtiness","and","medium_type","->","medium",
                 "small","and","medium_type","->","medium",
                 "large","and","not_greasy","->","medium",
                 "medium_dirtiness","and","not_greasy","->","short",
                 "small","and","not_greasy","->","very_short"),
               nrow=9, byrow=TRUE)
object <- frbs.gen(range.data, num.fvalinput, names.varinput, num.fvaloutput, varout.mf,
                   names.varoutput, rule, varinp.mf, type.model, type.defuz, type.tnorm,
                   type.snorm, func.tsk = NULL, colnames.var, type.implication.func, name)
par(mar=c(1,1,1,1))
summary(object)
plotMF(object)
#for(i in 1:(ncol(data)-1)){
# data[,i]=normalize(data[,i],0,1)
#}
res <- predict(object, newdata[,-ncol(newdata)])$predicted.val
print(res)
label<-vector()
labelEncoding<-vector()
for(i in 1:length(res)){
  for(j in 1:ncol(varout.mf)){
    lowerbound=varout.mf[2,j]
    upperbound=varout.mf[4,j]
    if(res[i]<=upperbound && res[i]>lowerbound ){
      label=append(label,varoutput.1[j])
      labelEncoding=append(labelEncoding,j)
    }
  }
}
print(label)
output<-factor(label)
expectedOutput<-factor(newdata[,ncol(newdata)])
levels(expectedOutput)<-varoutput.1
levels(output)<-varoutput.1
cm=confusionMatrix(newdata[,ncol(newdata)],output)
print(cm)
cat("Accuracy",cm$overall["Accuracy"],"\n")
if(is.na(cm$byClass['Pos Pred Value'])){
  cat("Precision",0,"\n" )
}
if(is.na(cm$byClass['Sensitivity'])){
  cat("Recall",0,"\n")
}
if(is.na(cm$byClass['Pos Pred Value'])||is.na(cm$byClass['Sensitivity'])){
  cat("F Measure",0,"\n")
}
plot(newdata[,1],newdata[,2],xlab="dirtiness",ylab="Type of dirt",col=labelEncoding,main="Dirtiness vs DirtType",pch=2)
legend("topleft",
       legend = varoutput.1, col=labelEncoding,pch=2)

dirt.csv
dirtiness,type of dirt,washing_time
30,75,medium
45,80,long
30,45,medium
55,85,very_long
25,23,very_short
65,25,short
15,32,very_short
80,80,long
90,90,long
5,95,medium
67,78,medium
77,11,short
17,24,very_short
6,40,short


## Sugeno Model ##

library(frbs)
#use if needed
normalize <- function(x,new_min,new_max) {
  return (((x - min(x)) / (max(x) - min(x)))*(new_max-new_min)+new_min)
}
varinp.mf <- matrix(c(1, 1000, 1200, 1400, NA, 1, 1400, 1600, 1800, NA,
                      1, 2.5, 5, 7.5, NA, 1, 7.5, 10, 14, NA),                  
                    nrow = 5, byrow = FALSE)

num.fvalinput <- matrix(c(2, 2), nrow=1)
varinput.1 <- c("smallW", "largeW")
varinput.2 <- c("smallT", "largeT")
names.varinput <- c(varinput.1, varinput.2)
range.data <- matrix(c(1000,2000, 0, 15,1000,15000), nrow=2)
type.defuz <- "WAM"
type.tnorm <- "MIN"
type.snorm <- "MAX"
type.implication.func <- "ZADEH"
name <- "Sim-0"
newdata<- read.csv("tsk.csv")
colnames.var <- c("W", "T","AU")
type.model <- "TSK"
num.fvaloutput <- matrix(c(2), nrow=1)
varoutput.1 <- c("small","large")
names.varoutput <- c(varoutput.1)
varout.mf <- matrix(c(1, 2000, 6000, 10000,NA ,
                      1, 10000, 12500, 15000, NA
),
nrow = 5, byrow = FALSE)
func.tsk <- matrix(c(4.6925, -526.2, 2631, 3.4765, -210.5, 2103, 4.6925, -526.2, 2631),
                   nrow = 3, byrow = TRUE)
rule <- matrix(c("largeW", "and", "smallT", "->",
                 "smallW", "or", "largeT", "->",
                 "smallW", "and", "smallT", "->"),
               nrow = 3, byrow = TRUE)
object <- frbs.gen(range.data, num.fvalinput, names.varinput,
                   num.fvaloutput = NULL, varout.mf = NULL, names.varoutput = NULL, rule,
                   varinp.mf, type.model, type.defuz = NULL, type.tnorm, type.snorm,
                   func.tsk, colnames.var, type.implication.func, name)
plotMF(object)
#for(i in 1:(ncol(data)-1)){
# data[,i]=normalize(data[,i],0,1)
#}
res <- predict(object, newdata[,-ncol(newdata)])$predicted.val
print(res)

label<-vector()
labelEncoding<-vector()
for(i in 1:length(res)){
  for(j in 1:ncol(varout.mf)){
    lowerbound=varout.mf[2,j]
    upperbound=varout.mf[4,j]
    if(res[i]<=upperbound && res[i]>lowerbound ){
      label=append(label,varoutput.1[j])
      labelEncoding=append(labelEncoding,j)
    }
  }
}
output<-factor(label)
expectedOutput<-factor(newdata[,ncol(newdata)])
levels(expectedOutput)<-varoutput.1
levels(output)<-varoutput.1
cm=confusionMatrix(expectedOutput,output)
print(cm)
cat("Accuracy",cm$overall["Accuracy"],"\n")
if(is.na(cm$byClass['Pos Pred Value'])){
  cat("Precision",0,"\n" )
}
if(is.na(cm$byClass['Sensitivity'])){
  cat("Recall",0,"\n")
}
if(is.na(cm$byClass['Pos Pred Value'])||is.na(cm$byClass['Sensitivity'])){
  cat("F Measure",0,"\n")
}

plot(newdata[,1],newdata[,2],xlab="W",ylab="T",col=labelEncoding,main="W vs T",pch=2)
legend("topleft",
       legend = varoutput.1, col=labelEncoding,pch=2)



''' w,t,class
1300,6.5,small
1200,12.5,small
1500,13,large
1300,12,large
1600,12.5,small
1100,4.5,small
1800,12,large
1700,10,large
1280,7,small
1100,14,large '''


## C Means ##

library(philentropy)
it=1
eps=3
m=2
no_of_clusters=3
no_of_dimensions=4
no_of_datapoints=150

datapoints<-c()
#for(i in 1:no_of_dimensions){
#  datapoints<-cbind(datapoints,round(runif(no_of_datapoints,min=1,max=10)))
#}
data(iris)
datapoints=cbind(iris[,1],iris[,2],iris[,3],iris[,4])

#step1
clusters<-c()
for(i in 1:no_of_clusters){
  clusters<-cbind(clusters,runif(no_of_datapoints))
}
for(i in 1:nrow(clusters)){
  clusters[i,]=clusters[i,]/sum(clusters[i,])
}
clusters

#datapoints<-cbind(c(1,2,3,4,5,6),c(6,5,8,4,7,9))
#clusters<-cbind(c(0.8,0.9,0.7,0.3,0.5,0.2),c(0.2,0.1,0.3,0.7,0.5,0.8))
centroid<-c()

repeat{
  
  #step2
  prev_centroid=centroid
  centroid<-c()
  for(i in 1:no_of_clusters){
    temp<-c()
    for(j in 1:no_of_dimensions){
      temp<-c(temp,sum((clusters[,i]^m)*datapoints[,j])/sum(clusters[,i]^m))
    }
    centroid<-rbind(centroid,temp)
  }
  rownames(centroid)<-c()
  centroid
  
  #step3
  dissimilarity<-c()
  for(i in 1:no_of_clusters){
    temp<-c()
    for(j in 1:no_of_datapoints){
      temp<-rbind(temp,distance(rbind(centroid[i,],datapoints[j,]),method="euclidean"))
    }
    dissimilarity<-cbind(dissimilarity,temp)
  }
  colnames(dissimilarity)<-c()
  dissimilarity
  
  #step4
  for(i in 1:nrow(dissimilarity)){
    dissimilarity[i,]=((1/dissimilarity[i,])/sum(1/dissimilarity[i,]))
  }
  dissimilarity
  
  cat('Iteration',it,'\n')
  print(cbind(datapoints,clusters,dissimilarity))
  #plot(1:150,rowMax(dissimilarity))
  
  clusters=dissimilarity
  
  #exit condition
  it=it+1
  if(it==3){
    flag=0
    for(i in 1:no_of_clusters){
      if(distance(rbind(centroid[i,],prev_centroid[i,]),method="euclidean")<eps){
        flag=flag+1
      }
    }
    if(flag==no_of_clusters){
      break
    }
  }
}

