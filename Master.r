## cmeans ##

library(clv)
m<-2
data <- read.csv("~/iris.csv")
points<-nrow(data)
dimensions<-ncol(data)-1
membershipMatrix<-matrix(nrow=points,ncol=dimensions)
optimalCluster<-function(){
  dunn1<-vector()
  db1<-vector()
  for(i in 2:10){
    agnesMod<-agnes(data[-ncol(data)])
    vPred<-as.integer(cutree(agnesMod,i))
    clsScatter<-cls.scatt.data(data[-ncol(data)],vPred,"euclidean")
    dunn1<-append(dunn1,clv.Dunn(clsScatter,"centroid","centroid"))
    db1<-append(db1,clv.Davies.Bouldin(clsScatter,"centroid","centroid"))
  }
  plot(c(2:10),dunn1,xlab="No of clusters",ylab="Dunn Index",type="o",main="Dunn Index Evaluation")
  plot(c(2:10),db1,xlab="No of clusters",ylab="Db Index",type="o",main="Db Index Evaluation")
  
}
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
distance<-function(centroids,clusters){
  distanceMatrix<-matrix(nrow=points,ncol=clusters)
  for(i in 1:points){
    for(j in 1:clusters){
      dist <- sqrt( sum( ( data[i,-ncol(data)] - centroids[j,] )^2 ) )
      distanceMatrix[i,j]<-dist
    }
  } 
  return(distanceMatrix)
}
membership<-function(distanceMatrix,clusters){
  
  membershipMatrix<-matrix(0,nrow=points,ncol=clusters)
  for(i in 1:nrow(distanceMatrix)){
    for(j in 1:ncol(distanceMatrix)){
      up=(1/distanceMatrix[i,j])^(1/(m-1))
      down=0
      for(x in 1:ncol(distanceMatrix)){
        down=down+(1/distanceMatrix[i,x])
      }
      down1=down^(1/(m-1))
      membershipMatrix[i,j]<-up/down1
      if(is.nan(membershipMatrix[i,j])){
        membershipMatrix[i,j]=0;
      }
      
    }
  }
  return(membershipMatrix)
}
calculateCentroid<-function(membershipMatrix,clusters){
  centroid<-matrix(nrow=clusters,ncol=dimensions)
  for(i in 1:nrow(centroid)){
    for(j in 1:ncol(centroid)){
      centroid[i,j]<-( sum(data[,j] * membershipMatrix[,i]^m)/sum(membershipMatrix[,i]^m) )
      centroid[i,j]<-round(centroid[i,j],3)
    }
  }
  return(centroid)
}
##for(i in 1:(ncol(data)-1)){
##  data[,i]=normalize(data[,i])
##}
cMeans<-function(clusters){
  prevCentroids=matrix(nrow=clusters,ncol=dimensions)
   for(i in 1:nrow(prevCentroids)){
     for(j in 1:ncol(prevCentroids))
        prevCentroids[i,j]<-round(runif(1,min(data[,j]),max(data[,j])),4)
  }
  centroids=prevCentroids
  flag<-TRUE
  while(flag){
    distanceMatrix<-distance(centroids,clusters)
    memberShipMatrix<-membership(distanceMatrix,clusters)
    centroids<-calculateCentroid(memberShipMatrix,clusters)
    if(identical(prevCentroids,centroids)){
      flag=FALSE
    }
    else{
      prevCentroids=centroids
    }
    i=i+1
  }
  print("FINAL CENTROIDS")
  print(centroids)
  clusterOutput<-vector()
  for(i in 1:nrow(memberShipMatrix)){
    clusterOutput=append(clusterOutput,which.max(memberShipMatrix[i,]))
  }
  plot(iris[c(1, 2)], col = data[,ncol(data)],main="Original clustered output")
  plot(iris[c(1, 2)], col = clusterOutput,main="CMeans clustered output")
  points(centroids[,c(1, 2)], col = 1:length(unique(clusterOutput)), pch = 8, cex = 2)
}
par(mfrow=c(2,2))
optimalCluster()
cMeans(3)


## classification ##


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
cat("Accuracy",cm$overall["Accuracy"],"\n")
if(is.na(cm$byClass['Pos Pred Value'])){
  cat("Precision",0,"\n" )
}else{
  cat("Precision",cm$byClass['Pos Pred Value'],"\n" )
}
if(is.na(cm$byClass['Sensitivity'])){
  cat("Recall",0,"\n")
}else{
  cat("Recall",cm$byClass['Sensitivity'],"\n")
}
if(is.na(cm$byClass['Pos Pred Value'])||is.na(cm$byClass['Sensitivity'])){
  cat("F Measure",0,"\n")
}else{
  Precision=cm$byClass['Pos Pred Value']
  Recall=cm$byClass['Sensitivity']
  Fmeasure=(2 * Precision * Recall) / (Precision + Recall)
  cat("F Measure",Fmeasure,"\n")
}
par(mfrow=c(1,2))
pie(x,labels,main="Count of records in each Class",col=rainbow(length(x)))
plot(data[,1],data[,2],xlab="Sepal Length",ylab="Sepal Width",col=output,main="Sepal Length vs SepalWidth",pch=2)
legend("topleft",
       legend = labels, col=c("black","red","green"),pch=2)


## mamdani ##

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
newdata<-read.csv("~/dirty.csv")
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
print("Defuzzified values")
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
print("Associated Classes")
print(label)
output<-factor(label)
expectedOutput<-factor(newdata[,ncol(newdata)])
levels(expectedOutput)<-varoutput.1
levels(output)<-varoutput.1
cm=confusionMatrix(expectedOutput,output)
print(cm)
cat("Accuracy",cm$overall["Accuracy"],"\n")
if(is.na(cm$byClass['Pos Pred Value'])){
  cat("Precision",0,"\n" )
}else{
  cat("Precision",cm$byClass['Pos Pred Value'],"\n" )
}
if(is.na(cm$byClass['Sensitivity'])){
  cat("Recall",0,"\n")
}else{
  cat("Recall",cm$byClass['Sensitivity'],"\n")
}
if(is.na(cm$byClass['Pos Pred Value'])||is.na(cm$byClass['Sensitivity'])){
  cat("F Measure",0,"\n")
}else{
 
    Precision=cm$byClass['Pos Pred Value']
    Recall=cm$byClass['Sensitivity']
    Fmeasure=(2 * Precision * Recall) / (Precision + Recall)
    cat("F Measure",Fmeasure,"\n")
  
}
Sys.sleep(10)
par(mfrow=c(1,1))
plot(newdata[,1],newdata[,2],xlab="dirtiness",ylab="Type of dirt",col=labelEncoding,main="Dirtiness vs DirtType",pch=2)
legend("topright",
      legend = varoutput.1, col=1:length(varoutput.1),pch=2,inset=c(0,0.35))


## sugeno ##

library(frbs)
library(caret)
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
newdata<- read.csv("~/tsk.csv")
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
par(mfrow=c(2,2))
plotMF(object)
#for(i in 1:(ncol(data)-1)){
 # data[,i]=normalize(data[,i],0,1)
#}
res <- predict(object, newdata[,-ncol(newdata)])$predicted.val

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
print("Defuzzified values ")
print(res)
print("Associated Class")
print(label)
output<-factor(label)
expectedOutput<-factor(newdata[,ncol(newdata)])
levels(expectedOutput)<-varoutput.1
levels(output)<-varoutput.1
cm=confusionMatrix(expectedOutput,output)
print("Confusion Matrix")
print(cm)
cat("Accuracy",cm$overall["Accuracy"],"\n")
if(is.na(cm$byClass['Pos Pred Value'])){
  cat("Precision",0,"\n" )
}else{
  cat("Precision",cm$byClass['Pos Pred Value'],"\n" )
}
if(is.na(cm$byClass['Sensitivity'])){
  cat("Recall",0,"\n")
}else{
  cat("Recall",cm$byClass['Sensitivity'],"\n")
}
if(is.na(cm$byClass['Pos Pred Value'])||is.na(cm$byClass['Sensitivity'])){
  cat("F Measure",0,"\n")
}else{
  Precision=cm$byClass['Pos Pred Value']
  Recall=cm$byClass['Sensitivity']
  Fmeasure=(2 * Precision * Recall) / (Precision + Recall)
  cat("F Measure",Fmeasure,"\n")
}
Sys.sleep(10)
plot(newdata[,1],newdata[,2],xlab="W",ylab="T",col=labelEncoding,main="W vs T",pch=1)
legend("bottomright",
       legend = varoutput.1, col=1:length(varoutput.1),pch=1,inset=c(0,0))

