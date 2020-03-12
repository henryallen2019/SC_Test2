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
