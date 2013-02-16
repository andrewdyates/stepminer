## Fit step up function as used in boolean implication mining.
## Andrew Yates 2013
sse <- function(v) sum((v-mean(v))^2)

## compute sum squared error per step location "i"
## 1:i -> low step
## i+1:n -> high step
fit.upstep <- function(v, do.plots=FALSE) {
  R <- list()
  v <- sort(v)
  n <- length(v)
  R$v.sse <- sapply(1:(n-1), function(i) sse(v[1:i])+sse(v[(i+1):n]))
  R$idx <- which.min(R$v.sse)
  R$sum.sse <- R$v.sse[R$idx]
  R$th <- (v[R$idx]+v[R$idx+1])/2
  R$low <- mean(v[1:R$idx])
  R$high <- mean(v[(R$idx+1):n])
  R$mean <- mean(v)
  R$mean.idx <- tail(which(v<=mean(v)), 1)
  R$median <- median(v)
  R$median.idx <- max(floor(length(v)/2),1)
  R
}

plot.sse <- function(R, add.mean.median=FALSE, ...) {
  plot(R$v.sse, xlab="Highest Index In Low Step", ylab="Sum of Sum Squared Errors", ...)
  abline(h=R$sum.sse, col="#cc0000", lty=1)
  abline(v=R$idx, col="#cc0000", lty=3)
  points(R$idx, R$v.sse[R$idx], col="#ff0000", pch=15, cex=3)
  if(add.mean.median) {
    # mean
    abline(v=R$mean.idx, col="#0000cc", lty=3)
    abline(h=R$v.sse[R$mean.idx], col="#0000cc", lty=1)
    points(R$mean.idx, R$v.sse[R$mean.idx], col="#0000cc", pch=16, cex=3)
    # median
    abline(v=R$median.idx, col="#009900", lty=3)
    abline(h=R$v.sse[R$median.idx], col="#009900", lty=1)
    points(R$median.idx, R$v.sse[R$median.idx], col="#009900", pch=17, cex=3)    
  }
}


plot.stepfit <- function(R, v, add.mean.median=FALSE, ...) {
  v <- sort(v)
  plot(v, xlab="Rank", ylab="Log Expression", ...)
  abline(h=R$th, col="#cc0000", lty=2, lwd=2)
  abline(v=R$idx, col="#990000", lty=3)
  abline(h=R$high, col="#990000", lwd=1, lty=3)
  abline(h=R$low, col="#990000", lwd=1, lty=3)
  if(add.mean.median) {
    abline(h=R$mean, col="#0000cc", lty=2, lwd=2)
    #abline(v=R$mean.idx, col="#0000cc", lty=3, lwd=1)
    points(R$mean.idx, R$mean, col="#0000cc", pch=16, cex=3)
    abline(h=R$median, col="#009900", lty=2, lwd=2)
    #abline(v=R$median.idx, col="#009900", lty=3, lwd=1)
    points(R$median.idx, R$median, col="#009900", pch=17, cex=3)
  }
  segments(x0=0, y0=R$low, x1=R$idx, y1=R$low, col="#ff0000", lwd=6)
  segments(x0=R$idx, y0=R$low, x1=R$idx, y1=R$high, col="#ff0000", lwd=6)
  segments(x0=R$idx, y0=R$high, x1=length(row1)+1, y1=R$high, col="#ff0000", lwd=6)
  points(R$idx+0.5, R$th, col="#ff0000", pch=15, cex=3)
}


