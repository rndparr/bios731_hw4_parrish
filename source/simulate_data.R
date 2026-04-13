

###############################################################
## SIMULATE DATA
###############################################################
simulate_data <- function(n, mu_true = c(0, 5, 10, 20), pi_true = rep(0.25, 4), sigma2 = 1) {
	K <- length(mu_true)
	sigma <- sqrt(sigma2)

	# sample clusters
	c_true <- sample(1:K, size = n, prob = pi_true, replace = TRUE)
	c_tab <- table(c_true)

	# initialize y
	y <- rep(NA, n)

	# get y values by cluster
	sapply(1:K, 
		function(i) { y[which(c_true == i)] <<- rnorm(c_tab[[i]], mu_true[i], sd = sigma); }) |>
		invisible()

	return(data.frame(y = y, c = c_true))
}