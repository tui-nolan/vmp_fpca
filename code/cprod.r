########## R function: cprod ##########

# For constructing a vector version of the crossprod function

# Created: 18 MAR 2020
# Last Updated: 13 APR 2020

cprod <- function(x, y) {
	
	if(missing(y)) {
		
		if(!is.vector(x)) {
			
			stop("Use the crossprod function for matrix inner products")
		}
		
		y <- x
	}
	
	if(!is.vector(y) & !is.vector(x)) {
			
		stop("Use the crossprod function for matrix inner products")
	}
	
	ans <- as.vector(crossprod(x, y))
	return(ans)
}

########## End of cprod ##########