
###############################################################
## REORDER MU, C
###############################################################
# sort mu smallest to largest, correct cluster values 
reorder_mu_c <- function(mu_est, c_est){
	# re-order mu, get index for reordering
	fix_order <- setNames(sort(mu_est, index.return=TRUE), c('mu_est', 'ind'))

	# re-order cs to correspond to corrected mu vector
	fix_order$c_est <- match(c_est, fix_order$ind)

	# export ind because may need to see that in gibbs results
	return(fix_order)
}

###############################################################
## FUNCTION TIMING
###############################################################
# function time
func_time <- function(x){
	system.time(x)[['elapsed']]
}
