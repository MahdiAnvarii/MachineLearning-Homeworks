# session 3
# preprocessing data #

data("iris")
head(iris)
dim(iris)
library(ggplot2)
iris <- iris[sample(1:nrow(iris)),]
ggplot(iris,aes(Sepal.Length,Sepal.Width,color=Species))+geom_point()
ggplot(iris,aes(Petal.Length,Petal.Width,color=Species))+geom_point()
w <- runif(5)
w
data.x <- cbind(iris[,1:4],1)
data.y <- ifelse(iris$Species=="versicolor",1,0)
head(data.x)
head(data.y)
N <- nrow(data.x)

train.x <-data.x[seq(N/2),]
train.y <-data.y[seq(N/2)]
test.x <-data.x[-seq(N/2),]
test.y <-data.y[-seq(N/2)]

# perceptron

perceptron.predict <- function(x,w){
  if (sum(w * x) <= 0)
    return(0)
  else return(1)
}

update2 <- function(x,w,y,lr){
  p <- perceptron.predict(x,w)
  e= y-p # error
  w <<- w + lr * e * x # update w in global namespace with lr as a rate
  NULL
}


lr <- 0.01
acct <- sapply(1:400,function(i){
  invisible(sapply(1:nrow(train.x), function(i) update2(as.numeric(train.x[i,]),w,train.y[i],lr)))

  pred <- sapply(1:nrow(test.x),function(i) perceptron.predict(as.numeric(test.x[i,]),w))

  acc <- sum(pred==test.y) / length(test.y)
  acc
})
plot(acct)

# batch update

update.minibatch <- function(x,w,y,lr){
  p <- sapply(1:nrow(x),function(i) perceptron.predict(as.numeric(x[i,]),w))
  e <- y-p
  w <<- w + lr * rowMeans(e*t(x))
}

mb.size <- 10
lr <- 0.01
acct <- sapply(1:400,function(i){
  mb <- sample(1:nrow(train.x),mb.size)
  
  update.minibatch(train.x[mb,],w,train.y[mb],lr)
  
  pred <- sapply(1:nrow(test.x),function(i) perceptron.predict(as.numeric(test.x[i,]),w))
  
  acc <- sum(pred==test.y) / length(test.y)
  acc
})
plot(acct)