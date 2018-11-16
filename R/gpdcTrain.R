#' GPD Classifier - training
#'
#' This function is used to train a GPD classifier. It can be used to perform open set classification based on the generalized Pareto distribution.
#' @param train a data matrix containing the train data. Class labels should not be included.
#' @param k the number of upper order statistics to be used.
#' @details For details on the method and parameters see Vignotto and Engelke (2018).
#' @return A list of three elements.
#' \item{pshapes}{the estimated rescaled shape parameters for each point in the training dataset.}
#' \item{balls}{the estimated radius for each point in the training dataset.}
#' \item{k}{the number of upper order statistics used.}
#' @note Data are not scaled internally; any preprocessing has to be done externally.
#' @author Edoardo Vignotto \cr
#' \email{edoardo.vignotto@unige.ch}
#' @references Vignotto, E., & Engelke, S. (2018). Extreme Value Theory for Open Set Classification-GPD and GEV Classifiers. \emph{arXiv preprint arXiv:1808.09902}.
#' @seealso \code{\link{gpdcTest}}
#' @examples
#' trainset <- LETTER[1:15000,]
#' knowns <- trainset[trainset$class==1, -1]
#' gpdClassifier <- gpdcTrain(train = knowns, k = 10)
#' @export

gpdcTrain <- function(train,k) {

  stopifnot(is.data.frame(train))
  for(i in 1:ncol(train)) {
    stopifnot(is.numeric(train[,i]))
  }
  stopifnot(is.numeric(k) & k==(k%/%1))

  q <- 1-1/(1*k)
  p <- ncol(train)

  distances <- -RANN::nn2(train,train,k=k+2)$nn.dist[,-1]
  threshold <- distances[,k+1]
  distances <- distances[,-(k+1)]

  R <- (distances)/threshold


  shape <- apply(R,1,function(x) mean(log(x[x != 0])))
  ball <- apply(cbind(shape,threshold),1,function(x)
    -evd::qgpd(q,x[2],x[1]*x[2],x[1]))
  pshape <- p*shape


  return(list(pshapes=pshape, balls=ball, k=k))

}
