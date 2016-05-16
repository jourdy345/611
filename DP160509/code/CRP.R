crp = function(num, M, fun) {
	cluster = c(1)	#cluster number
	n = c(1)	#number of elements in cluster
	theta=fun(1)
	
	for (index in 1:(num-1)){
		p = c( n / (sum(n) + M), M / (sum(n) + M) )
		#cat("p_",index, ": ", p, "\n")
		i = sample(x = 1:(length(cluster)+1), size = 1, prob = p, replace = TRUE)
		if (i == length(cluster) + 1){
			newTheta = fun(1)
			theta=c(theta,newTheta)
			n=c(n,1)
			cluster=1:(length(cluster) + 1)
		} else {
			n[i]=n[i]+1
		}
	}
	return(list(n=n,cluster=cluster,theta=theta))
}
a = crp(10,5,runif); n <- a$n
t(as.data.frame(n,row.names=as.character(round(a$theta,2))))
barplot(a$n, xlab=expression(theta), ylab='n', main='Predictive Distribution', 
        names.arg=round(a$theta,2))
