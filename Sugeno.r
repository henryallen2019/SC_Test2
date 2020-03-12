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
