

## GET ELBO
get_elbo <- function(y, pi_i, mu_star, sigma2_star, sigma2) {
	n <- length(y)
	K <- length(mu_star)

	prior_term <- -0.5 * sum( log(20 * pi * sigma2) + ( (mu_star^2 + sigma2_star) / sigma2) ) - n * log(K)
	lik_term <- sum(pi_i * (-0.5 * log(2 * pi) - 0.5 * (-2 * sum(outer(y, mu_star)) + n * sum(mu_star^2 + sigma2_star))))
	entropy_c <- - sum(pi_i * log(pi_i + 1e-12)) # with numerical stability fix
	entropy_mu <- 0.5 * sum(log( (2 * pi_i * sigma2_star) + 1e-12 ) + 1) # with numerical stability fix

	return(prior_term + lik_term + entropy_c + entropy_mu)
}


## GET PI_I
get_pi_i <- function(y, mu_star, sigma2_star) {
	# initial log pi_i value
	log_pi_i <- ( outer(y, mu_star) - 0.5 * outer(rep(1, n), mu_star^2 + sigma2_star))

	# fix recommended by Blei, et al for stability
	row_max_log_pi_i <- apply(log_pi_i, 1, max)
	pi_i <- exp(log_pi_i - row_max_log_pi_i)
	pi_i <- pi_i / rowSums(pi_i) # normalize

	return(pi_i)
}


## CAVI MODEL
cavi_model <- function(y, K, sigma2, mu_star = NULL, sigma2_star = NULL, tol = 1e-6, max_iter = 100) {
	n <- length(y)

	# fix index
	max_iter <- max_iter + 1

	# initialize parameters
	if (is.null(mu_star)){
		mu_star <- ((max(y) - min(y)) / (K - 1) * 0:(K - 1)) + min(y)
	}

	if (is.null(sigma2_star)){
		sigma2_star <- rep(sigma2, K)
	}

	# initiate pi_i
	pi_i <- get_pi_i(y, mu_star, sigma2_star)

	# initiate elbo vector
	elbo <- rep(NA, max_iter)
	elbo[1] <- get_elbo(y, pi_i, mu_star, sigma2_star, sigma2)

	# iterations
	converged <- FALSE
	for (iter in 2:max_iter) {
		nk <- colSums(pi_i)

		# update mu_star, sigma2_star
		sigma2_star <- 1 / (nk + 1 / sigma2)
		mu_star <- sigma2_star * colSums(pi_i * y)

		# update pi_i
		pi_i <- get_pi_i(y, mu_star, sigma2_star)
		
		# get elbo
		elbo[iter] <- get_elbo(y, pi_i, mu_star, sigma2_star, sigma2)

		# test tolerance
		if ( abs(elbo[iter] - elbo[iter - 1]) <= tol ) {
			converged <- TRUE
			break
		}
	}

	# remove NA vals from elbo
	elbo <- elbo[!is.na(elbo)]

	# cluster estimate set to cluster with highest prob
	c_est <- apply(pi_i, 1,  which.max)


	# fix order (mu small to large)
	out <- reorder_mu_c(mu_star, c_est)

	return(list(mu_est = out$mu_est, c_est = out$c_est, elbo = elbo, n_iter = length(elbo) - 1, converged = converged))
}




