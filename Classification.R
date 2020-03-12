normalize <- function(x,new_min,new_max) {
  return (((x - min(x)) / (max(x) - min(x)))*(new_max-new_min)+new_min)
}

composition <- function(param_matrix) {
  # Initialize an empty matrix
  result = matrix(nrow = nrow(param_matrix), ncol = ncol(param_matrix))
  for(i in c(1:nrow(param_matrix))) {
    for(j in c(1:ncol(param_matrix))) {
      result[i,j] = max(pmin(param_matrix[i,],param_matrix[,j]))
    }
  }
  return(result)
}

isEquivalenceMatrix <- function(param_matrix) {
  flag = TRUE;
  # Check for reflexive property
  for(i in c(1:ncol(param_matrix))) {
    if(param_matrix[i,i] != 1) {
      flag = FALSE
    }
  }
  # Check for symmetry
  for(i in c(1:nrow(param_matrix))) {
    for(j in c(1:ncol(param_matrix))) {
      if(param_matrix[i,j] != param_matrix[j,i]) {
        flag = FALSE
      }
    }
  }
  # Check for transitive property
  combos = combn(4, 3)
  for(i in c(1: ncol(combos))) {
    l = param_matrix[combos[1,i], combos[2, i]]
    r = param_matrix[combos[2,i], combos[3, i]]
    mid = param_matrix[combos[1,i], combos[3, i]]
    if(mid < min(l, r))
    {
      flag = FALSE;
    }
  }
  return (flag);
}

data=read.csv("~/iris.csv")
for(i in 1:(ncol(data)-1)){
  data[,i]=normalize(data[,i],0,1)
}
print(data)
# Initial Matrix Data
no_of_columns <- 4
table <- data

similarityMatrix = matrix(c(0), nrow = no_of_columns, ncol = no_of_columns)

for(i in c(1:4)){
  for(j in c(i: 4))
  {
    similarityMatrix[i, j] = sum(table[,i] * table[,j])/(sqrt(sum(table[,i]^2) * sum(table[,j]^2)))
    similarityMatrix[j , i] = similarityMatrix[i, j]
  }
}

finished = FALSE

while(!finished) {
  if(!isEquivalenceMatrix(similarityMatrix)){
    similarityMatrix = composition(similarityMatrix)
    print(similarityMatrix)
  }else
  {
    finished = TRUE
    cut = 0.8
    similarityMatrix[similarityMatrix < cut] <- 0
    similarityMatrix[similarityMatrix >= cut] <- 1
    print('GROUPS')
    for(i in c(1: ncol(similarityMatrix))) {
      position = which(similarityMatrix[,i] == 1)
      if(length(position) > 0)
        print(position)
      similarityMatrix[position,] = 0;
    }
    break;
  }
}
