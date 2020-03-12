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
