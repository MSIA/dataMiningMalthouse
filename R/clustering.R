#' Summarizing K-Means Fits
#'
#' summary method for class "kmeans".
#'
#' @param fit an object of class "kmeans", usually, a result of a call to kmeans
#' @export
summary.kmeans = function(fit)
{
  p = ncol(fit$centers)
  K = nrow(fit$centers)
  n = sum(fit$size)
  xbar = t(fit$centers)%*%fit$size/n
  print(data.frame(
    n=c(fit$size, n),
    Pct=(round(c(fit$size, n)/n,2)),
    round(rbind(fit$centers, t(xbar)), 2),
    RMSE = round(sqrt(c(fit$withinss/(p*(fit$size-1)), fit$tot.withinss/(p*(n-K)))), 4)
  ))
  cat("SSE=", fit$tot.withinss, "; SSB=", fit$betweenss, "; SST=", fit$totss, "\n")
  cat("R-Squared = ", fit$betweenss/fit$totss, "\n")
  cat("Pseudo F = ", (fit$betweenss/(K-1))/(fit$tot.withinss/(n-K)), "\n\n");
  invisible(list(Rsqr=fit$betweenss/fit$totss,
                 F=(fit$betweenss/(K-1))/(fit$tot.withinss/(n-K))) )
}

#' Plot Method for K-Means Fits
#'
#' summary method for class "kmeans".
#'
#' @param fit an object of class "kmeans", usually, a result of a call to kmeans
#' @param boxplot whether to plot box-and-whisker around data points
#' @export
plot.kmeans = function(fit, boxplot=F)
{
  require(lattice)
  p = ncol(fit$centers)
  k = nrow(fit$centers)
  plotdat = data.frame(
    mu=as.vector(fit$centers),
    clus=factor(rep(1:k, p)),
    var=factor( 0:(p*k-1) %/% k, labels=colnames(fit$centers))
  )
  print(dotplot(var~mu|clus, data=plotdat,
                panel=function(...){
                  panel.dotplot(...)
                  panel.abline(v=0, lwd=.1)
                },
                layout=c(k,1),
                xlab="Cluster Mean"
  ))
  invisible(plotdat)
}

#' Find best k for K-Means Clustering
#'
#' Search for the best k to perform k-means clustering on a data matrix.
#'
#' @param x numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param Ks a vector of candidate K values to perform K-Means Clustering
#' @param nstart number of times to repeat with different random seeds for each k
#' @export
search.k <- function(x, Ks, nstart = 100) {
  F <- double(length(Ks))
  offset <- Ks[1] - 1
  for (K in Ks) {
    fit <- kmeans(x, K, nstart = nstart)
    F[K - offset] <- summary(fit)$F
    sse[K - offset] = fit$tot.withinss

    # silhouette
    si2 = silhouette(fit$cluster, dist(x, "euclidean"))
    si[K - offset] = summary(si2)$avg.width
  }
  par(mfrow=c(1, 3))
  plot(Ks, F, type = "b", xlab = "Number Clusters K")
  plot(Ks, si, type="b")
  plot(Ks, sse, type="b")
}
