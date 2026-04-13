

###############################################################
## SAMPLE CLUSTER
###############################################################
sample_cluster <- function(y, mu) {
	K <- length(mu)
	n <- length(y)

	# probability vector pi_i for cluster sample
	pi_i <- sapply(1:K, 
		function(k){
			sapply(y, function(y_i) {
				exp(-0.5 * (y_i - mu[k])^2) 
			})
		})
	
	# normalize
	pi_i <- pi_i / rowSums(pi_i)

	# sample clusters
	c_sample <- sapply(1:n, 
		function(i) { sample(1:K, size = 1, prob = pi_i[i, ]) })

	return(c_sample)
}

###############################################################
## SAMPLE MU
###############################################################
sample_mu <- function(y, c, K, sigma2) {
	# matrix of n_k, sum(y_k) values
	nk_sum_yk <- sapply(1:K, 
		function(k) {
			y_k <- y[c == k]; 
			return( rbind(length(y_k), sum(y_k)) ) 
		})

	# mu_star, sigma2_star values
	sigma2_star <- 1 / ( nk_sum_yk[1, 1:K] + (1 / sigma2) )
	mu_star <- sigma2_star * nk_sum_yk[2, 1:K]

	# sample mus
	mu_sample <- sapply(1:K, 
		function(k){ rnorm(1, mu_star[k], sqrt(sigma2_star[k])) } )

	mu_sample <- sort(mu_sample)

	return(mu_sample)
}

###############################################################
## GIBBS SAMPLER - 1 chain
###############################################################
gibbs_chain <- function(y, K, sigma2, n_iter = 10000, burnin = 2000) {
	n <- length(y)

	# fix index
	n_iter <- n_iter + 1
	burnin <- burnin + 1

	# initialize matrices
	c <- matrix(0, nrow = n_iter, ncol = n)
	mu <- matrix(0, nrow = n_iter, ncol = K)

	# initial parameters
	c[1, ] <- sample(1:K, n, replace = TRUE)
	mu[1, ] <- rep(mean(y), K) #((max(y) - min(y)) / (K - 1) * 0:(K - 1)) + min(y)

	# update
	for (iter in 2:n_iter) {
		c[iter, ] <- sample_cluster(y, mu[iter - 1, ])
		mu[iter, ] <- sample_mu(y, c[iter - 1, ], K, sigma2)
	}

	# get estimates (sans burnin samples)
	mu_est <- colMeans(mu[(burnin + 2):n_iter, ])
	c_est <- apply(c[(burnin + 2):n_iter, ], 2, function(x){ which.max(tabulate(x)) })

	# confidence intervals
	mu_ci_l <- apply(mu[(burnin + 2):n_iter, ], 2, quantile, probs = 0.025)
	mu_ci_u <- apply(mu[(burnin + 2):n_iter, ], 2, quantile, probs = 0.975)

	# fix order (mu small to large)
	# out <- reorder_mu_c(mu_est, c_est)

	# return(list(mu = mu, c = c, mu_est = out$mu_est, c_est = out$c_est, ind_fix = out$ind, burnin = burnin - 1))

	return(list(mu = mu, c = c, mu_est = mu_est, mu_ci_l = mu_ci_l, mu_ci_u = mu_ci_u, c_est = c_est, burnin = burnin - 1))
}


###############################################################
## GIBBS SAMPLER-MULTIPLE CHAINS
###############################################################

# add time output to an individual chain (for use with gibbs_multi_chain)
timed_gibbs_chain <- function(y, K, sigma2, n_iter = 10000, burnin = 2000) {
	# run chain, timed
	gibbs_time <- func_time(
		chain_res <- gibbs_chain(y, K, sigma2, n_iter = n_iter, burnin = burnin)
		)

	# add time to output
	chain_res$time <- gibbs_time

	return(chain_res)
}


# run gibbs sampler with multiple chains
gibbs_multi_chain <- function(y, K, sigma2, n_chains = 4, n_iter = 10000, burnin = 2000){

	# replicate gibbs_chain n_chains times
	res <- setNames(
		replicate(n_chains, timed_gibbs_chain(y, K, sigma2, n_iter = n_iter, burnin = burnin), simplify=FALSE), 
		paste0('chain_', 1:n_chains))

	# summary
	res[['burnin']] <- burnin
	res[['mu_est']] <- sapply(1:n_chains, function(j) { res[[paste0('chain_', j)]]$mu_est }) |> rowMeans()

	return(res)
}


