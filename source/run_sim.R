

###############################################################
## PARSE ARGUMENTS
###############################################################
# Rscript run_sim_i.R ${i} ${ncores}

# parse arguments
args=(commandArgs(TRUE))
cat(paste0('\nargs:\n'))
print(args)
cat()
if(length(args) == 0) {
	stop('Error: No arguments supplied!')
} else  {
	i <- as.integer(args[[1]])
	n <- as.integer(args[[2]])
}

# print session info to log
print(sessionInfo())

###############################################################
## SOURCE FUNCTIONS
###############################################################
print('Importing source scripts.')
source(here::here('source', 'utility.R'))
source(here::here('source', 'simulate_data.R'))
source(here::here('source', 'gibbs_functions.R'))
source(here::here('source', 'cavi_functions.R'))

###############################################################
# PARAMETERS
###############################################################
K <- 4
mu_true <- c(0, 5, 10, 20)
pi_true <- rep(0.25, 4)
sigma2 <- 1
n_chain <- 4

mu_cols <- c(paste0('mu_', 1:4))
mu_ci_cols <- c(paste0('mu_', 1:4, '_ci_l'), paste0('mu_', 1:4, '_ci_u'))
mu_c_cols <- c(mu_cols, paste0('c_', 1:n))

# set seed
seed <- 041126 + n + i
set.seed(seed)

###############################################################
# SIMULATE DATA
###############################################################
print('Simulating data.')

# dataframe includes y and true c 
y_dat <- simulate_data(n, mu_true = mu_true, pi_true = pi_true, sigma2 = sigma2)

# true dat
true_dat <- data.frame(rbind(c(mu_true, y_dat$c)))
colnames(true_dat) <- mu_c_cols


# y-vector
y <- y_dat$y

###############################################################
# FIT
###############################################################
# Gibbs sampler
print('Running Gibbs Sampler.')
gibbs_time <- func_time(
	gibbs_res <- gibbs_chain(y, K, sigma2, n_iter = 10000, burnin = 2000)
	)
gibbs_res$time <- gibbs_time

# CAVI
print('Running CAVI.')
cavi_time <- func_time(
	cavi_res <- cavi_model(y, K, sigma2, mu_star = NULL, sigma2_star = NULL, tol = 1e-6, max_iter = 500)
	)
cavi_res$time <- cavi_time

###############################################################
# OUTPUT
###############################################################
print('Saving output.')

# output all method output
save(true_dat, gibbs_res, cavi_res, 
	file = here::here('sim', 'sim_data', paste0('sim_n', n, '_', i, '.Rds')))

# output results data
out <- rbind(
	# rbind(c(i, n, 'true', NA, mu_true)), # true values
	rbind(c(i, n, 'gibbs', gibbs_time, gibbs_res$mu_est, gibbs_res$mu_ci_l, gibbs_res$mu_ci_u)), # Gibbs values
	rbind(c(i, n, 'cavi', cavi_time, cavi_res$mu_est, cavi_res$mu_ci_l, cavi_res$mu_ci_u)) # CAVI values
	)
out <- as.data.frame(out)
colnames(out) <- c('i', 'n', 'method', 'time', mu_cols, mu_ci_cols)

out[, c('time', mu_cols, mu_ci_cols)] <- mapply(as.numeric, out[, c('time', mu_cols, mu_ci_cols)])

# bias
for (k in 1:K){
	out[[paste0('mu_', k, '_bias')]] <-  out[[paste0('mu_', k)]] - mu_true[k]
}

# coverage
for (k in 1:K){
	out[[paste0('mu_', k, '_coverage')]] <- (mu_true[k] >= out[[paste0('mu_', k, '_ci_l')]]) & (mu_true[k] <= out[[paste0('mu_', k, '_ci_u')]])
}

# entire mu vector covered
out[['coverage']] <- apply(out[, paste0('mu_', 1:4, '_coverage')], 1, all)

# append to dataframe
write.table(
	out,
	here::here('sim', 'sim_out', paste0('sim_n', n, '_i', i, '.txt')),
	quote = FALSE,
	append = FALSE,
	row.names = FALSE,
	col.names = FALSE,
	sep = '\t')

print('Done!')


