#' GPD Classifier - testing
#'
#' This function is used to evaluate a test set for a pre-trained GPD classifier. It can be used to perform open set classification based on the generalized Pareto distribution.
#' @param train  data matrix containing the train data. Class labels should not be included.
#' @param test a data matrix containing the test data.
#' @param pre a list obtained with the function \code{\link{gpdcTrain}}.
#' @param prob logical indicating whether p-values should be returned.
#' @param alpha threshold to be used if \code{prob} is equal to \code{FALSE}. It must be between 0 and 1.
#' @details For details on the method and parameters see Vignotto and Engelke (2018).
#' @return If \code{prob} is equal to \code{TRUE}, a vector containing the p-values for each point is returned. A high p-value results in the classification of the corresponding test data as a known point, since this hypothesis cannot be rejected. If the p-value is small, the corresponding test data is classified as an unknown point. If \code{prob} is equal to \code{TRUE}, a vector of predicted values is returned.
#' @author Edoardo Vignotto \cr
#' \email{edoardo.vignotto@unige.ch}
#' @references Vignotto, E., & Engelke, S. (2018). Extreme Value Theory for Open Set Classification-GPD and GEV Classifiers. \emph{arXiv preprint arXiv:1808.09902}.
#' @seealso \code{\link{gpdcTrain}}
#' @examples
#' trainset <- LETTER[1:15000,]
#' testset <- LETTER[-(1:15000), -1]
#' knowns <- trainset[trainset$class==1, -1]
#' gpdClassifier <- gpdcTrain(train = knowns, k = 10)
#' predicted <- gpdcTest(train = knowns, test = testset, pre = gpdClassifier)
#' @export

gpdcTest <- function(train,test,pre,prob = TRUE,alpha = 0.01) {

  stopifnot(is.data.frame(train))
  for(i in 1:ncol(train)) {
    stopifnot(is.numeric(train[,i]))
  }
  stopifnot(is.data.frame(test))
  for(i in 1:ncol(test)) {
    stopifnot(is.numeric(test[,i]))
  }
  stopifnot(ncol(test)==ncol(train))
  stopifnot(alpha>0 & alpha<1)
  stopifnot(length(pre)==3)

  p <- ncol(train)
  k <- pre$k
  pshapes <- pre$pshapes
  balls <- pre$balls
  q <- 1-1/(1*k)
  output <- rep(0,dim(test)[1])
  distances <- -RANN::nn2(train,test,k=k+1)$nn.dist
  threshold <- distances[,k+1]
  distances <- distances[,-(k+1)]
  R <- distances/threshold
  shape <- apply(R,1,function(x) mean(log(x[x != 0])))
  ball <- apply(cbind(shape,threshold),1,function(x)
    -evd::qgpd(q,x[2],x[1]*x[2],x[1]))
  pshape <- p*shape


  #for(i in 1:dim(test)[1]) output[i] <- min(mean(balls>ball[i]),mean(pshapes>pshape[i]))
  for(i in 1:dim(test)[1]) output[i] <- mean(pshapes>pshape[i] | (pshapes<pshape[i] & (balls>ball[i])))

  if(!prob) {

    quant <- seq(1,1-alpha,length.out=length(pshapes)-round(stats::quantile(1:length(pshapes),1-alpha)))
    quant <- quant[quant != 1]

    z <- NA
    minz <- Inf

    for(i in quant) {
      #minTemp <- abs(min(mean(balls>stats::quantile(balls, i)),mean(pshapes>stats::quantile(pshapes, i)))-alpha)
      minTemp <- abs(mean(pshapes>stats::quantile(pshapes, i) | (pshapes<stats::quantile(pshapes, i) & (balls>stats::quantile(balls, i))) )-alpha)
      if(minTemp<minz) {
        minz <- minTemp
        z <- i
      }
    }

    shape_threshold <- stats::quantile(pshapes,z)
    ball_threshold <- stats::quantile(balls,z)

    output <- rep("unknown",dim(test)[1])
    for(i in 1:length(pshape)) {
      if((pshape[i]<shape_threshold) & ball[i]<ball_threshold) {
        output[i] <- "known"
      }
    }
    output <- factor(output, levels=c("known","unknown"))
  }

  output
}
