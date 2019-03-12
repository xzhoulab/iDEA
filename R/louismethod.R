######################################################################
## iDEA project
## Editor: Shiquan Sun
## Date: 2019-2-4 14:23:13
## Modified: 2019-3-11 18:15:06
## Affiliation: University of Michigan
##				Department of Biostatistics				
######################################################################
# 
#' Correct the p-values using Louis method
#'
#' A part of variance is missed when we use EM algorithm, since we should use the complete-data gradient vector or second derivative matrix rather than the incomplete data likelihood. 
#' See details: "Thomas A. Louis. Finding the Observed Information Matrix when Using the EM Algorithm, Journal of the Royal Statistical Society. 
#' Series B (Methodological) Vol. 44, No. 2 (1982), pp. 226-233".
#' 
#' @param object iDEA object
#' @param ... Ignored
#'
#' @return Returns a iDEA object with corrected p-values stored in object@louis.
#' 
#' @import doSNOW
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import stats
#' @import utils
#' 
#' @export
#'
#'
LouisCorrect <- function(object){
	num_core <- object@num_core
	# using parallel to correct
	#library(doSNOW)
	if(num_core > 1){
		if(num_core > detectCores()){warning("LOUIS:: the number of cores you're setting is larger than detected cores!");num_core = detectCores()}
	}#end fi
	cl <- makeCluster(num_core)
	registerDoSNOW(cl)
	num_annot <- length(object@emmcmc)
	
	cat(paste("## ===== LOUIS METHOD ==== ##\n"))
	cat("## correcting p-values ...\n")
	
	pb <- txtProgressBar(max = num_annot, style = 3)
	progress <- function(n) setTxtProgressBar(pb, n)
	opts <- list(progress = progress)
	
	iannot<-1
	# parallel
	res_all <- foreach(iannot=1:num_annot, .combine=rbind, .options.snow=opts) %dopar% {
		res <- object@emmcmc[[iannot]]
		##################
		## Louis function
		LouisMethod <- function(res){

			numer <- (1-res$pip)*res$pip*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta1*res$sigma2_e))*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta2*res$sigma2_e))

			denom <- res$pip*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta1*res$sigma2_e)) + 
			(1-res$pip)*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta2*res$sigma2_e))
	
			# define missing information matrix
			T <- matrix(0, ncol=2, nrow=2)
			T[1,1] <- sum( as.numeric((numer/denom^2) * res$A[,1]*res$A[,1] ) )
			T[1,2] <- sum( as.numeric((numer/denom^2) * res$A[,1]*res$A[,2] ) )
			T[2,1] <- sum( as.numeric((numer/denom^2) * res$A[,1]*res$A[,2] ) )
			T[2,2] <- sum( as.numeric((numer/denom^2) * res$A[,2]*res$A[,2] ) )
	
			annot_var <- diag( solve(res$info_mat - T) )
			# return the results
			return(annot_var)
		}# end function
		#################
		
		if(abs(res$annot_coef[2])<5 & !is.na(res$annot_coef[2]) & !is.null(res$annot_coef) ){
			try_test <- try( louis_var <- LouisMethod(res), silent=T)
			if(class(try_test)=="try-error"){
				print("try-error!")
			}else{
				reseach <- data.frame(annot_id=names(object@emmcmc)[iannot], annot_coef=res$annot_coef[2], annot_var=res$annot_var[2], annot_var_louis=louis_var[2], sigma2_b=res$sigma2_beta1)
			}# end fi
		}# end fi
	}# end louis
	close(pb)
	stopCluster(cl) 

	# remove the negative variance genes
	res_all <- res_all[which(res_all$annot_var_louis>0), ]
	zscore <- res_all$annot_coef/sqrt(res_all$annot_var)
	zscore_louis <- res_all$annot_coef/sqrt( res_all$annot_var_louis )
	res_all$pvalue_louis <- 2*pnorm(-abs(zscore_louis))
	res_all$pvalue <- 2*pnorm(-abs(zscore))
	
	object@louis <- res_all
	rm(res_all)
	# return the results
	return(object)
}#end funcs




