## R functions


# By Nora Mitchell 
standardize <- function(x) {
  mu <- mean(x, na.rm=TRUE)
  sigma <- sd(x, na.rm=TRUE)
  y <- (x - mu)/sigma
  return(y)
}


drop.levels <- function(dat) {
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}


get.complete.cases <- function(x){
  ok <- complete.cases(x)
  x <- x[ok,]
  x <- drop.levels(x)
  return(x)
}


plot.estimates <- function(x) {
  if (class(x) != "summary.mcmc")
    x <- summary(x)
  n <- dim(x$statistics)[1]
  par(mar=c(4, 14, 3, 1))
  plot(x$statistics[,1], n:1,
       yaxt="n", ylab="",
       xlim=range(x$quantiles)*1.2,
       pch=19,
       cex.lab=1.2,
       xlab="Posterior means and 95% credible intervals")
  grid()
  axis(2, at=n:1, rownames(x$statistics), las=2)
  arrows(x$quantiles[,1], n:1, x$quantiles[,5], n:1, code=0)
  abline(v=0, lty=2)
}


# Ethan Linck's coefficient of variation function
cv <- function(data){
  mean <- mean(data)
  sd <- sd(data)
  coef <- sd/mean
  return(coef)
}