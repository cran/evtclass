#' GEV Classifier - testing
#'
#' This function is used to evaluate a test set for a pre-trained GEV classifier. It can be used to perform open set classification based on the generalized Pareto distribution.
#' @usage gevcTest(train, test, pre, prob = TRUE, alpha)
#' @param train a data matrix containing the train data. Class labels should not be included.
#' @param test a data matrix containing the test data.
#' @param pre a numeric vector of parameters obtained with the function \code{\link{gevcTrain}}.
#' @param prob logical indicating whether p-values should be returned.
#' @param alpha threshold to be used if \code{prob} is equal to \code{FALSE}. It must be between 0 and 1.
#' @details For details on the method and parameters see Vignotto and Engelke (2018).
#' @return If \code{prob} is equal to \code{TRUE}, a vector containing the p-values for each point is returned. A high p-value results in the classification of the corresponding test data as a known point, since this hypothesis cannot be rejected. If the p-value is small, the corresponding test data is classified as an unknown point. If \code{prob} is equal to \code{TRUE}, a vector of predicted values is returned.
#' @author Edoardo Vignotto \cr
#' \email{edoardo.vignotto@unige.ch}
#' @references Vignotto, E., & Engelke, S. (2018). Extreme Value Theory for Open Set Classification-GPD and GEV Classifiers. \emph{arXiv preprint arXiv:1808.09902}.
#' @seealso \code{\link{gevcTrain}}
#' @examples
#' trainset <- LETTER[1:15000,]
#' testset <- LETTER[-(1:15000), -1]
#' knowns <- trainset[trainset$class==1, -1]
#' gevClassifier <- gevcTrain(train = knowns)
#' predicted <- gevcTest(train = knowns, test = testset, pre = gevClassifier)
#' @export


gevcTest <- function(train,test,pre,prob=TRUE,alpha = 0.01) {

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
  stopifnot(length(pre)==2)

  distances <- 100*RANN::nn2(train,test,k=2)$nn.dist[,-1]
  out <- sapply(distances,function(x) 1-stats::pweibull(x,pre[1],pre[2]))
  if(!prob) {
    output <- rep("unknown", length(out))
    output[out>alpha] <- "known"
    out <- as.factor(output)
  }
  out
}
