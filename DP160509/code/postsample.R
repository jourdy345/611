postsample <- function(n, M) {
  ##data
  x <- rnorm(n,mean = 2, sd = 2)
  p <- rep(1/(M+n), n)
  base <- rnorm(1)
  p <- c( M/(M+n), p)
  x <- c(base, x) 

  ##posterior sample
  post <- sample(size = 1000, x = x, prob = p, replace = TRUE)
  return(post)
}

##plot cdf
par(mfrow=c(1,2))

plot(ecdf(postsample(8,5)), verticals = TRUE, do.points = FALSE, 
     lwd = 3, col = 'blue', ylab = expression(F(theta)), xlab = expression(theta),
     main = 'n=8', xlim = c(-4,8), ylim = c(0,1))
curve(pnorm(x,mean=2, sd=2), from =-4, to = 8, add = TRUE,lwd=2)
curve(pnorm, from =-4, to = 8, add = TRUE, lty = 2)
legend(2,0.2, legend=c("true G", "baseline G0", "E(G | y)"),
         lty=c(1,2,1), lwd=c(2,1,3), col=c('black', 'black', 'blue'), bty="n")

plot(ecdf(postsample(50,5)), verticals = TRUE, do.points = FALSE, 
     lwd = 2, col = 'blue', ylab = expression(F(theta)), xlab = expression(theta),
     main = 'n=50', xlim = c(-4,8), ylim = c(0,1))
curve(pnorm(x,mean=2, sd=2), from =-4, to = 8, add = TRUE,lwd=2)
curve(pnorm, from =-4, to = 8, add = TRUE, lty = 2)
legend(2,0.2, legend=c("true G", "baseline G0", "E(G | y)"),
         lty=c(1,2,1), lwd=c(2,1,3), col=c('black', 'black', 'blue'), bty="n")



