#' GEV Classifier - training
#'
#' This function is used to train a GEV classifier. It can be used to perform open set classification based on the generalized extreme value distribution.
#' @param train a data matrix containing the train data. Class labels should not be included.
#' @details For details on the method and parameters see Vignotto and Engelke (2018).
#' @return A numeric vector of two elements containing the estimated parameters of the fitted reversed Weibull.
#' @note Data are not scaled internally; any preprocessing has to be done externally.
#' @author Edoardo Vignotto \cr
#' \email{edoardo.vignotto@unige.ch}
#' @references Vignotto, E., & Engelke, S. (2018). Extreme Value Theory for Open Set Classification - GPD and GEV Classifiers. \emph{arXiv preprint arXiv:1808.09902}.
#' @seealso \code{\link{gevcTest}}
#' @examples
#' trainset <- LETTER[1:15000,]
#' knowns <- trainset[trainset$class==1, -1]
#' gevClassifier <- gevcTrain(train = knowns)
#' @export

gevcTrain <- function(train) {

  stopifnot(is.data.frame(train))
  for(i in 1:ncol(train)) {
    stopifnot(is.numeric(train[,i]))
  }

  distances <- 100*RANN::nn2(train,train,k=2)$nn.dist[,-1]
  if(sum(distances == 0) != 0) {
    distances <- distances[distances != 0]
  }
  fitdistrplus::mledist(distances, "weibull",optim.method="Nelder-Mead")$estimate
}
