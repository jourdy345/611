StickBreaking = function(num,n,M, fun) {
  theta = fun(n)
  beta = rbeta(n,1,M)
  p = c()
  p[1] = beta[1]
  p = c(beta[1], beta[2:n] * cumprod(1-beta[1:(n-1)]))
  cluster = sample(size=num, x=theta, replace=T, prob=p)
  return(cluster)
}

par(mfrow=c(1,2))

plot(ecdf(StickBreaking(1000, 100, 5, rnorm)), verticals = TRUE, do.points = FALSE,col='blue',
     ylab=expression(Fn (theta)),xlab=expression(theta),main='',xlim=c(-3,3),ylim=c(0,1))
#lines(ecdf(StickBreaking(1000,100,0.2,rnorm)), col='red',type='l')
curve(pnorm, from =-3, to = 3, add = TRUE,lwd=2)
legend('topleft', legend='M=5')

plot(ecdf(StickBreaking(1000, 100, 8, rnorm)), verticals = TRUE, do.points = FALSE,col='blue',
     ylab=expression(Fn (theta)),xlab=expression(theta),main='',xlim=c(-3,3),ylim=c(0,1))
#lines(ecdf(StickBreaking(1000,100,4,rnorm)), col='red',type='l')
curve(pnorm, from =-3, to = 3, add = TRUE,lwd=2)
legend('topleft', legend='M=8')

plot(ecdf(StickBreaking(1000, 100, 10, rnorm)), verticals = TRUE, do.points = FALSE,col='blue',
     ylab=expression(Fn (theta)),xlab=expression(theta),main='',xlim=c(-3,3),ylim=c(0,1))
#lines(ecdf(StickBreaking(1000,100,20,rnorm)), col='red',type='l')
curve(pnorm, from =-3, to = 3, add = TRUE,lwd=2)
legend('topleft', legend='M=10')

plot(ecdf(StickBreaking(1000, 100, 100, rnorm)), verticals = TRUE, do.points = FALSE,col='blue',
     ylab=expression(Fn (theta)),xlab=expression(theta),main='',xlim=c(-3,3),ylim=c(0,1))
#lines(ecdf(StickBreaking(1000,100,100,rnorm)), col='red',type='s')
curve(pnorm, from =-3, to = 3, add = TRUE,lwd=2)
legend('topleft', legend='M=100')