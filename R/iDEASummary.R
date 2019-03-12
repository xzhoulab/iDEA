####################################################################################################
## Package : iDEA
## Version : 1.0.0
## Date    : 2019-2-4 14:31:30
## Modified: 2019-3-11 18:12:49
## Title   : Integrative Differential Expression and Gene Set Enrichment Analysis
## Authors : Shiquan Sun, Ying Ma and Xiang Zhou
## Contacts: shiquans@umich.edu and yingma@umich.edu 
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
	# iDEA summary data
	res_idea <- NULL
	if(fit_null){
		# fit the model under the null
		if(verbose){cat(paste("## fitting the model without annotation ... \n"))}
		Annot <- as.matrix(data.frame(rep(1, num_gene)))
		t1 <- system.time(model1 <- try( res <- EMMCMCStepSummary(object@summary[,1], object@summary[,2], as.matrix(Annot), init_beta, init_tau[1], em_iter, mcmc_iter, min_degene) ))
		if(class(model1) != "try-error"){
			res$A <- Annot
			colnames(res$A) <- c("intercept")
			rownames(res$pip) <- object@gene_id
			res$converged   <- TRUE
			res$ctime   <- t1[3]
		}else{res <- NULL}
		object@null_model <- res
	}#end fi
	# fit the model under the alternative
	if(verbose){cat(paste("## fitting the model with annotation ... \n"))}
	res_idea <- pbmclapply(1:num_annot, FUN = function(x) { 
		Annot <- rep(0, object@num_gene)
		Annot[object@annotation[[x]]] <- 1
		Annot <- Annot - mean(Annot)
		Annot <- as.matrix(data.frame(rep(1, num_gene), Annot) )
		t1 <- system.time(model1 <- try( res <- EMMCMCStepSummary(object@summary[,1], object@summary[,2], as.matrix(Annot), init_beta, init_tau, em_iter, mcmc_iter, min_degene) ))
		if(class(model1) != "try-error"){
			res$A <- Annot
			colnames(res$A) <- c("intercept", names(object@annotation)[x])
			rownames(res$pip) <- object@gene_id
			res$converged   <- TRUE
			res$ctime   <- t1[3]
		}else{res <- NULL}
		return(res)}, mc.cores = getOption("mc.cores", num_core)
	)# end parallel
	
	names(res_idea) <- object@annot_id
	object@emmcmc <- res_idea
	rm(res_idea)
	# return results
	#object@louis <- LouisMethod(object, num_core=num_core)
	return(object)
}# end function 

######################################################
#             Obtain summary FIT FUNCTION			 #
######################################################
#' Get summary data from the existing DE method
#' @param counts scRNAseq raw counts matrix,  genes (row) by cells (column) 
#' @param celltype cell types for counts matrix, two cell types are allowed
#' @param method DE method to obtain the summary data
#' 
#' @import Matrix
#' @import edgeR
#' @import DESeq2
#' @import zinbwave
#' @import zingeR
#' @return Returns summary data for iDEA
#'
#' @export
GetSummaryData <- function(counts, celltype, method="DESeq2"){
	# convert into integer type
	iCounts <- apply(counts, 2, function(x){
		#x <- as.numeric(x)
		storage.mode(x) <- 'integer'
		return(x) }
	)
	rm(counts)
	if( tolower(method) == "edger"){
		d <- DGEList(iCounts)
		d <- suppressWarnings(edgeR::calcNormFactors(d))
		design <- model.matrix(~celltype)
		weights <- zingeR::zeroWeightsLS(counts = d$counts, design = design, maxit = 500, normalization = "TMM", verbose = TRUE)
		d$weights <- weights
		d <- estimateDisp(d, design)
		fit <- glmFit(d,design)
		res <- zinbwave::glmWeightedF(fit, ZI=TRUE, independentFiltering=TRUE)
		zscore <- sqrt( qchisq(1-res$table$PValue,df=1) )
		beta_se <- res$table$logFC/zscore
		summaryD <- data.frame(beta=res$table$logFC, beta_var=beta_se*beta_se)
		rownames(summaryD) <- rownames(res$table)
		return(summaryD)
	}else if( tolower(method) == "deseq2"){
		require(DESeq2)
		colData <- data.frame(traits = celltype)
		design <- model.matrix(~celltype)
		dse <- DESeq2::DESeqDataSetFromMatrix(countData = iCounts, colData = colData, design = ~traits)
		weights <- zingeR::zeroWeightsLS(counts = iCounts, design = design, maxit = 500, normalization = "DESeq2_poscounts", colData = colData, designFormula = ~traits, verbose = TRUE)
		assays(dse)[["weights"]] <- weights
		dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
		dse <- DESeq2::estimateDispersions(dse)
		dse <- DESeq2::nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-2)
		res <- as.data.frame( DESeq2::results(dse) )
		summaryD <- data.frame(beta=res$log2FoldChange, beta_var=res$lfcSE*res$lfcSE)
		rownames(summaryD) <- rownames(res)
		return(summaryD)
	}else{
		warning("Only two DE methods (DESeq2 and edgeR) are provided to obtain the summary data.")
	}# end fi
}# end func

######################################################
#             iDEA with weight FIT FUNCTION				 #
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
#' @return Returns a iDEA object with EM-MCMC results stored in object@emmcmc.
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
			res$converged <- TRUE
			res$ctime <- t1[3]
		}else{res <- NULL}
		return(res)}, mc.cores = getOption("mc.cores", num_core)
	)# end parallel
	
	names(res_idea) <- object@annot_id
	object@emmcmc <- res_idea
	rm(res_idea)
	## return results
	return(object)
}# end function 

#########################################
#             CODE END                  #
#########################################
