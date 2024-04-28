library(ggplot2)
library(reshape2)
library(e1071)
library(class)
library(randomForest)


isch <-read.delim("isch.txt")
head(isch)
tail(isch)
dim(na.omit(isch))
isch <- na.omit(isch)
summary(isch)
colnames(isch)
isch <- isch[,c(1:5,8,9)]
hist(isch$E)
hist(isch$a)
isch <- subset(isch, a<0.2)
sum(isch$a >0.2)
dim(isch)

x <- isch[, 4:7]
head(x)

xm <- melt(x)
head(xm)
ggplot(xm, aes(variable, value))+geom_boxplot()
ggplot(xm, aes(variable, value))+geom_violin()
ggplot(xm, aes(variable, value))+geom_violin(width=1 ,aes(fill=variable)) + geom_boxplot(width=0.05)

ggplot(isch, aes(A, E , color=Ischemia)) + geom_point() + theme_bw()

# classification
classify2 <-function(x,y){
  if(1.5*x-0.3 > y) return("POS")
  else return("NEG")
}
classify2(1,0.4)

classify <-function(x,y) ifelse(1.5*x-0.3 > y, ("POS") , ("NEG"))


isch$pred <- classify(isch$A , isch$E)
head(isch)

# evaluation
accuracy = sum(isch$pred==isch$Ischemia) / nrow(isch)
accuracy
sensitivity = sum(isch$pred== "POS" & isch$Ischemia=="POS") / sum(isch$Ischemia== "POS")
sensitivity
specificity = sum(isch$pred== "NEG" & isch$Ischemia=="NEG") / sum(isch$Ischemia== "NEG")
specificity

# ROC
classify3 <-function(x,y,arz) ifelse(1.5*x-arz > y, ("POS") , ("NEG"))
evalu <-function(arz){
  pred <- classify3(isch$A,isch$E ,arz)
  sens <- sum(pred== "POS" & isch$Ischemia=="POS") / sum(isch$Ischemia== "POS")
  spec <- sum(pred== "NEG" & isch$Ischemia=="NEG") / sum(isch$Ischemia== "NEG")
  return(c(sens,spec))
}
res <- sapply(seq(-2,2,by=0.01),evalu)
res <- t(res)
dim(res)
head(res)
tail(res)
res[,2] <- 1-res[,2]
res <- data.frame(res)
colnames(res) <- c("sens","spec")
ggplot(res , aes(spec , sens)) +geom_point()





# session 2
head(isch)
isch2 <- isch[,c(1:7)]
N <- nrow(isch2)
isch2 <-isch2[sample(1:N),]
train.set <- isch2[seq(N/2),]
test.set <- isch2[-seq(N/2),]
dim(train.set)
dim(test.set)
head(isch2)
head(train.set)
str(train.set)

# important
isch2["Ischemia"] <- lapply(isch2["Ischemia"] , factor)
isch2["Diastolic.function"] <- lapply(isch2["Diastolic.function"] , factor)
train.set <- isch2[seq(N/2),]
test.set <- isch2[-seq(N/2),]
# svm
model.svm <- svm(Ischemia~.,train.set)
test.pred <- predict(model.svm , test.set[, -2])
head(test.pred)
head(test.set)

# evaluation
accuracy = sum(test.pred==test.set$Ischemia) / nrow(test.set)
accuracy
sensitivity = sum(test.pred== "POS" & test.set$Ischemia=="POS") / sum(test.set$Ischemia== "POS")
sensitivity
specificity = sum(test.pred== "NEG" & test.set$Ischemia=="NEG") / sum(test.set$Ischemia== "NEG")
specificity

class(train.set$Ischemia)
levels(train.set$Ischemia)
table(train.set$Ischemia)

#randomforest
model.RF <- randomForest(train.set[,-2],train.set[,2])
test.pred2 <- predict(model.RF,test.set[, -2])
head(test.pred2)

# evaluation
accuracy = sum(test.pred2==test.set$Ischemia) / nrow(test.set)
accuracy
sensitivity = sum(test.pred2== "POS" & test.set$Ischemia=="POS") / sum(test.set$Ischemia== "POS")
sensitivity
specificity = sum(test.pred2== "NEG" & test.set$Ischemia=="NEG") / sum(test.set$Ischemia== "NEG")
specificity

senspeacc <- function(pred , real){
  acc = sum(pred==real) / length(pred)
  sen = sum(pred== "POS" & real=="POS") / sum(real== "POS")
  spe = sum(pred== "NEG" & real=="NEG") / sum(real== "NEG")
  list(Acc=acc,Sen=sen,Spe=spe)
}
senspeacc(test.pred2,test.set$Ischemia)

#cross validation
N <- nrow(isch2)
isch2 <-isch2[sample(1:N),]
section <- sample(1:5,size=N,replace = T)
head(section,20)

crossvalidation <- function(i){
  test.set= isch2[section==i,]
  train.set= isch2[section!=i,]
  train.set["Ischemia"] <- lapply(train.set["Ischemia"] , factor)
  model.svm <- svm(Ischemia~.,train.set, type="C-classification")
  test.pred <- predict(model.svm , test.set[, -2])
  senspeacc(test.pred,test.set$Ischemia)
}
result <- sapply(1:5, crossvalidation)
result
result <- data.frame(t(result))
class(result$Acc)
result$Acc <- as.numeric(unlist(result$Acc))
result$Sen <- as.numeric(unlist(result$Sen))
result$Spe <- as.numeric(unlist(result$Spe))
is.numeric(result$Acc)
colMeans(result)


# feature selection
N <- nrow(isch2)
isch2 <-isch2[sample(1:N),]
isch2["Ischemia"] <- lapply(isch2["Ischemia"] , factor)
isch2["Diastolic.function"] <- lapply(isch2["Diastolic.function"] , factor)
isch2$Diastolic.function <- as.numeric(isch2$Diastolic.function)
train.set <- isch2[seq(N/2),]
test.set <- isch2[-seq(N/2),]

train.x <- train.set[,-2]
train.y <- train.set[,2]
test.x <- test.set[,-2]
test.y <- test.set[,2]

k=2
feature.select <- function(features){
  model.svm <- svm(train.x[,features],train.y, type="C-classification")
  test.pred <- predict(model.svm , test.x[, features])
  senspeacc(test.pred,test.y)
}
head(train.x)
head(train.y)
head(test.x)
head(test.y)

class(isch2$Ischemia)
class(isch2$Diastolic.function)

feature.select(2)
results <- sapply(1:ncol(train.x),feature.select)
results
colnames(train.x)

feature.set <- combn(seq(ncol(train.x)),k)
feature.set
results2 <- sapply(1:ncol(feature.set), function(i) feature.select(feature.set[,i]))
results2

# feature extraction
fs <- combn(4:7 , 2) # E A e a
fs
multiply <- sapply(seq(ncol(fs)), function(i) isch2[,fs[1,i]]*isch2[,fs[2,i]])
head(multiply)
colnames(multiply) <- sapply(seq(ncol(fs)), function(i) paste0(colnames(isch2)[fs[1,i]],"*",colnames(isch2)[fs[2,i]]))

div1 <- sapply(seq(ncol(fs)), function(i) isch2[,fs[1,i]]/isch2[,fs[2,i]])
head(div1)
colnames(div1) <- sapply(seq(ncol(fs)), function(i) paste0(colnames(isch2)[fs[1,i]],"/",colnames(isch2)[fs[2,i]]))

div2 <- sapply(seq(ncol(fs)), function(i) isch2[,fs[2,i]]/isch2[,fs[1,i]])
head(div2)
colnames(div2) <- sapply(seq(ncol(fs)), function(i) paste0(colnames(isch2)[fs[2,i]],"/",colnames(isch2)[fs[1,i]]))

newisch <- cbind(isch,multiply,div1,div2)
head(newisch)

N <- nrow(newisch)
newisch <-newisch[sample(1:N),]
newisch["Ischemia"] <- lapply(newisch["Ischemia"] , factor)
newisch["Diastolic.function"] <- lapply(newisch["Diastolic.function"] , factor)
newisch$Diastolic.function <- as.numeric(newisch$Diastolic.function)
train.set <- newisch[seq(N/2),]
test.set <- newisch[-seq(N/2),]

train.x <- train.set[,-2]
train.y <- train.set[,2]
test.x <- test.set[,-2]
test.y <- test.set[,2]

feature.set <- combn(seq(ncol(train.x)),k)
feature.set
ggplot(newisch,aes(A,a,color=Ischemia)) + geom_point()
nrow(newisch)