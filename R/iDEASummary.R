####################################################################################################
## Package : iDEA
## Version : 1.0.1
## Date    : 2019-2-4 14:31:30
## Modified: 2019-9-25 21:17:54
## Title   : Integrative Differential Expression and Gene Set Enrichment Analysis
## Authors : Shiquan Sun
## Contacts: shiquans@umich.edu
##           University of Michigan, Department of Biostatistics
####################################################################################################

#' Integrative DE and GSEA based summary statistics
#'
#' Integrative DE and GSEA into a unfied framwork. The parameters of the model were infered via EM-MCMC algorithm.
#'
#' @param object iDEA object
#' @param fit_null Bool variable to indicate whether fitting the model without the annoation
#' @param init_beta Initial value for gene effect size, beta
#' @param init_tau Initial value for annotations, including the intercept
#' @param min_degene The threshold for the number of detected DE genes. For some of extremely cases, 
#' the method does not work when the number of detected DE genes is 0.
#' @param em_iter Maximum iteration for EM algorithm
#' @param mcmc_iter Maximum iteration for MCMC algorithm
#' @param fit.tol Tol for fitting the model
#' @param verbose Print the progresses
#' @param ... Ignored
#' 
#' @import pbmcapply
#' @import stats
#' @import utils
#' @return Returns a iDEA object with EM-MCMC results stored in object@emmcmc.
#'
#' @export
#'
#'
#'
iDEA.fit <- function(object,
					fit_null=FALSE,
				    init_beta=NULL, 
					init_tau=c(-2,0.5),
					min_degene=5,
					em_iter=15, 
					mcmc_iter=1000, 
					fit.tol=1e-5,  
					verbose=TRUE, ...) {
  # number of cores
	num_core <- object@num_core
	
	if(num_core > 4){
		if(num_core>detectCores()){warning("iDEA:: the number of cores you're setting is larger than detected cores!");num_core = detectCores()}
	}#end fi
	
	# On Windows, set cores to be 1 because of mclapply
	if(.Platform$OS.type == "windows"){num_core <- 1}
	#cl <- makeCluster(num_core)
	#registerDoSNOW(cl)
	#registerDoParallel(cores=num_core)
	# number of genes
	num_gene <- object@num_gene
	num_annot <- length(object@annotation)

	cat(paste("## ===== iDEA INPUT SUMMARY ==== ##\n"))
	cat(paste("## number of annotations: ", num_annot,"\n"))
	cat(paste("## number of genes: ", num_gene,"\n"))
	cat(paste("## number of cores: ", num_core,"\n"))
	
	if(is.null(init_beta)){
		#init_beta <- rep(0, num_gene)
		init_beta <- object@summary[,1] # to speed up, initialized with estimated effect size
	}# end fi
	#************************#
	#       EM - MCMC        #
	#************************#
	# iDEA summary data, gsea
	res_idea <- NULL
	if(fit_null){
		# fit the model under the null
		if(verbose){cat(paste("## fitting the null model ... \n"))}
		Annot <- as.matrix(data.frame(rep(1, num_gene)))
		t1 <- system.time(model1 <- try( res <- EMMCMCStepSummary(object@summary[,1], object@summary[,2], as.matrix(Annot), init_beta, init_tau[1], em_iter, mcmc_iter, min_degene) ))
		if(class(model1) != "try-error"){
			rownames(res$pip) <- object@gene_id
			res$converged   <- TRUE
			res$ctime   <- t1[3]
		}else{res <- NULL}
		object@null_model <- res
	}#end fi
	# fit the model under the alternative
	if(verbose){cat(paste("## fitting the alternative model ... \n"))}
	res_idea <- pbmclapply(1:num_annot, FUN = function(x) { 
		Annot <- rep(0, object@num_gene)
		Annot[object@annotation[[x]]] <- 1
		Annot <- Annot - mean(Annot)
		Annot <- as.matrix(data.frame(rep(1, num_gene), Annot) )
		t1 <- system.time(model1 <- try( res <- EMMCMCStepSummary(object@summary[,1], object@summary[,2], as.matrix(Annot), init_beta, init_tau, em_iter, mcmc_iter, min_degene) ))
		if(class(model1) != "try-error"){
			rownames(res$pip) <- object@gene_id
			colnames(res$pip) <- "PIP"
			rownames(res$beta) <- object@gene_id
			colnames(res$beta) <- "beta"
			rownames(res$annot_coef) <- c("tau_1", "tau_2")
			colnames(res$annot_coef) <- "annot_coef"
			rownames(res$annot_var) <- c("tau_1", "tau_2")
			colnames(res$annot_var) <- "annot_var"
			res$converged   <- TRUE
			res$ctime   <- t1[3]
		}else{res <- NULL}
		return(res)}, mc.cores = getOption("mc.cores", num_core)
	)# end parallel
	
	names(res_idea) <- object@annot_id
	object@de <- res_idea
	rm(res_idea)
	# return results
	#object@louis <- LouisMethod(object, num_core=num_core)
	return(object)
}# end function 


######################################################
#             iDEA with weight FIT FUNCTION		     #
######################################################
#' Integrative DE and GSEA based summary statistics and weight for each gene
#'
#' Integrative DE and GSEA into a unfied framwork. The parameters of the model were infered via EM-MCMC algorithm.
#'
#' @param object iDEA object
#' @param weight The weight for each gene to construct the gene-gene interation
#' @param fit_null Bool variable to indicate whether fitting the model without the annoation
#' @param init_beta Initial value for gene effect size, beta
#' @param init_tau Initial value for annotations, including the intercept
#' @param min_degene The threshold for the number of detected DE genes. For some of extremely cases, 
#' the method does not work when the number of detected DE genes is 0.
#' @param em_iter Maximum iteration for EM algorithm
#' @param mcmc_iter Maximum iteration for MCMC algorithm
#' @param fit.tol Tol for fitting the model
#' @param verbose Print the progresses
#' @param ... Ignored
#'
#' @import pbmcapply
#'
#' @return Returns a iDEA object with EM-MCMC results stored in object@de.
#'
#' @export
#'
#'
iDEAWeight.fit <- function(object,
					weight=NULL,
					fit_null=FALSE,
					init_beta=NULL, 
					init_tau=c(-2,0.5),
					min_degene=5,
					em_iter=20, 
					mcmc_iter=1000, 
					fit.tol=1e-5,  
					verbose=TRUE, ...) {

	##number of cores
	num_core <- object@num_core
	if(num_core > 4){
		if(num_core>detectCores()){warning("iDEAw:: the number of cores you're setting is larger than detected cores!");num_core = detectCores()}
	}#end fi
	## On Windows, set cores to be 1
	if(.Platform$OS.type == "windows"){num_core <- 1}

	## number of genes
	num_gene <- object@num_gene
	num_annot <- length(object@annotation)

	cat(paste("## ===== iDEA-WEIGHT INPUT SUMMARY ==== ##\n"))
	cat(paste("## number of annotations: ", num_annot,"\n"))
	cat(paste("## number of genes: ", num_gene,"\n"))
	cat(paste("## number of cores: ", num_core,"\n"))
	
	if(is.null(init_beta)){
		init_beta <- rep(0, num_gene)
	}# end fi
	## check weights -- a vector
	if(is.null(weight)){
		stop("iDEAw:: the weight is equired to pre-compute before fitting the idea.weight function! Please call 'weight <- ComputeWeight(counts)' ")
	}# end fi
	object@weight <- weight
	#************************#
	#      EM - MCMC        #
	#************************#
	##iDEA summary data
	res_idea <- NULL
	if(fit_null){
		## fit the model under the null
		Annot <- as.matrix(data.frame(rep(1, num_gene)))
		t1 <- system.time(model1 <- try( res <- EMMCMCStepSummaryWeight(object@summary[,1], object@summary[,2], weight, as.matrix(Annot), init_beta, init_tau[1], em_iter, mcmc_iter, min_degene) ))
		if(class(model1) != "try-error"){
			rownames(res$pip) <- object@gene_id
			res$converged   <- TRUE
			res$ctime   <- t1[3]
		}else{res <- NULL}
		object@null_model <- res
	}# end fi
	## fit the model under the alternative
	res_idea <- pbmclapply(1:num_annot, FUN = function(x){ 
		Annot <- rep(0, object@num_gene)
		Annot[object@annotation[[x]]] <- 1
		Annot <- Annot - mean(Annot)
		Annot <- as.matrix(data.frame(rep(1, num_gene), Annot) )
		t1 <- system.time(model1 <- try( res <- EMMCMCStepSummaryWeight(object@summary[,1], object@summary[,2], weight, as.matrix(Annot), init_beta, init_tau, em_iter, mcmc_iter, min_degene) ))
		if(class(model1) != "try-error"){
			rownames(res$pip) <- object@gene_id
			colnames(res$pip) <- "PIP"
			rownames(res$beta) <- object@gene_id
			colnames(res$beta) <- "beta"
			rownames(res$annot_coef) <- c("tau_1", "tau_2")
			colnames(res$annot_coef) <- "annot_coef"
			rownames(res$annot_var) <- c("tau_1", "tau_2")
			colnames(res$annot_var) <- "annot_var"
			res$converged <- TRUE
			res$ctime <- t1[3]
		}else{res <- NULL}
		return(res)}, mc.cores = getOption("mc.cores", num_core)
	)# end parallel
	
	names(res_idea) <- object@annot_id
	object@de <- res_idea
	rm(res_idea)
	## return results
	return(object)
}# end function 

######################################################
#             LOUIS METHOD TO CORRECT P-VALUES       #
######################################################
#' Correct the p-values using Louis method
#'
#' A part of variance is missed when we use EM algorithm, since we should use the complete-data gradient vector or second derivative matrix rather than the incomplete data likelihood. 
#' See details: "Thomas A. Louis. Finding the Observed Information Matrix when Using the EM Algorithm, Journal of the Royal Statistical Society. 
#' Series B (Methodological) Vol. 44, No. 2 (1982), pp. 226-233".
#' 
#' @param object iDEA object
#' @param ... Ignored
#'
#' @return Returns a iDEA object with corrected p-values stored in object@gsea.
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
iDEA.louis <- function(object){
	num_core <- object@num_core
	
	# using parallel to correct
	#library(doSNOW)
	if(num_core > 1){
		if(num_core > detectCores()){warning("LOUIS:: the number of cores you're setting is larger than detected cores!");num_core = detectCores()}
	}#end fi
	cl <- makeCluster(num_core)
	registerDoSNOW(cl)
	num_annot <- length(object@de)
	
	pb <- txtProgressBar(max = num_annot, style = 3)
	progress <- function(n) setTxtProgressBar(pb, n)
	opts <- list(progress = progress)
	
	
	# parallel
	res_all <- foreach(i=1:num_annot, .combine=rbind, .options.snow=opts) %dopar% {
		res <- object@de[[i]]
		Annot <- rep(0, object@num_gene)
		Annot[object@annotation[[i]]] <- 1
		Annot <- Annot - mean(Annot)
		Annot <- as.matrix(data.frame(rep(1, object@num_gene), Annot) )
		##################
		## Louis function
		LouisMethod <- function(res, A){
			numer <- (1-res$pip)*res$pip*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta1*res$sigma2_e))*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta2*res$sigma2_e))

			denom <- res$pip*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta1*res$sigma2_e)) + 
			(1-res$pip)*pnorm(res$beta, mean=0, sd=sqrt(res$sigma2_beta2*res$sigma2_e))
	
			# define missing information matrix
			T <- matrix(0, ncol=2, nrow=2)
			T[1,1] <- sum( as.numeric((numer/denom^2) * A[,1]*A[,1] ) )
			T[1,2] <- sum( as.numeric((numer/denom^2) * A[,1]*A[,2] ) )
			T[2,1] <- sum( as.numeric((numer/denom^2) * A[,1]*A[,2] ) )
			T[2,2] <- sum( as.numeric((numer/denom^2) * A[,2]*A[,2] ) )
	
			annot_var <- diag( solve(res$info_mat - T) )
			# return the results
			return(annot_var)
		}# end function
		#################
		
		if(abs(res$annot_coef[2])<5 & !is.na(res$annot_coef[2]) ){
			try_test <- try( louis_var <- LouisMethod(res, Annot), silent=T)
			if(class(try_test)=="try-error"){
				print("try-error!")
			}else{
				reseach <- data.frame(annot_id=object@annot_id[i], annot_coef=res$annot_coef[2], annot_var=res$annot_var[2], annot_var_louis=louis_var[2], sigma2_b=res$sigma2_beta1)
			}# end fi
		}# end fi
	}# end louis
	close(pb)
	stopCluster(cl) 

	# remove the negative variance genes
	pos_index <- which(res_all$annot_var_louis>0)
	zscore <- res_all$annot_coef[pos_index]/sqrt(res_all$annot_var[pos_index])
	zscore_louis <- res_all$annot_coef[pos_index]/sqrt( res_all$annot_var_louis[pos_index] )
	res_all$pvalue_louis <- 2*pnorm(-abs(zscore_louis))
	res_all$pvalue <- 2*pnorm(-abs(zscore))
	
	object@gsea <- res_all
	rm(res_all)
	# return the results
	return(object)
}#end funcs

#########################################
#             CODE END                  #
#########################################
