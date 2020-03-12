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
''' 30,75,medium
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
6,40,short '''
