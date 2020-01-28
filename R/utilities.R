########################################
## iDEA Project
## Editor : Shiquan Sun, Ying Ma
## Date   : 2019-2-4 13:43:16
## Modified: 2020-01-21 08:55:27
########################################

## to creat a new object for the project
setClass("iDEA", slots=list(
	summary = "data.frame",
	annotation = "list",
	project = "character",
	gene_id = "character",
	annot_id = "character",
	weight = "data.frame",
	num_gene = "numeric",
	num_core = "numeric",
	noGS = "list",
	de = "list",
	gsea = "data.frame",
    BMA_pip = "data.frame"
) )


#' Setup the iDEA object, 
#'
#' @param summary Summary statistics (the estimated gene effect size and its variance, g1 x 2 matrix) from single-cell RNAseq differential expression tools, i.g. zingeR, MAST, etc. The summary statistics file should be in data.frame data format with gene name as row names, while the column names are not required but the order of the column matters: the first column should be the coefficient and the second column should be the variance of the coefficient for each gene.
#' @param annotation Pre-definted gene specific annotations (gene sets, g2 x m matrix), i.e., GO term, pathways etc. The gene specific annotation file is required data.frame data format with gene name as row names, while the header is the gene specific annotation name. The row names of annotation file should match exactly with row names of summary statistics and should be in the same type, i.e. gene symbol or transcription id etc. If not, one solution is to use biomaRt R package to convert the gene name to make sure they are consistent each other.
#' @param project Project name (string).
#' @param max_var_beta Include genes where the variances which are smaller than 'max_var_beta' are maintained.
#' @param min_precent_annot The threshold of coverage rate (CR), i.e., the number of annotated genes (gene set size) divided by the number of tested genes.
#' @param num_core Number of cores for parallel implementation.
#' @import parallel
#' @import methods
#' 
#' @return Returns a iDEA object with the summary statistics stored in object@summary and object@annotation.
#' object@project, object@num_core, object@num_gene, object@num_annot, object@gene_id, object@annot_id, are also initialized.
#'
#' @export
#'

CreateiDEAObject <- function(summary, annotation, project = "iDEA", max_var_beta = 100, min_precent_annot = 0.0025, num_core = 10){
	
	## check data fromat
	if(!is.data.frame(summary)){
		summary <- as.data.frame(summary)
	}# end fi
	colnames(summary) <- c("beta", "beta_var")
	keep_index <- summary$beta_var<max_var_beta & !is.na(summary$beta_var)
	summary <- as.data.frame(summary[keep_index,])
	
	if(!is.data.frame(annotation)){
		annotation <- as.data.frame(annotation)
	}# end fi
	
	# On Windows, set cores to be 1 because of mclapply
	if(.Platform$OS.type == "windows"){num_core <- 1}
	
	## merge two data with gene names
	if(nrow(summary) != nrow(annotation)){
		summary$GeneID <- rownames(summary)
		annotation$GeneID <- rownames(annotation)
		combined_data <- merge(summary, annotation, by="GeneID")
		summary <- combined_data[, 2:3]
		rownames(summary) <- combined_data$GeneID
		annotation <- combined_data[, 4:ncol(combined_data)]
		rownames(annotation) <- combined_data$GeneID
		rm(combined_data)
	}# end fi
	
	## number of genes
	num_gene <- nrow(annotation)
	precent_annot <- apply(annotation, 2, sum)/num_gene
	## filtering out at specificed precentage genes are annotated, 
	annotation <- annotation[, which(precent_annot>min_precent_annot)] #  
    annotation = as.data.frame(annotation)
	## convert annotation into list
	annot_list <- mclapply(1:ncol(annotation), FUN = function(x) { which(annotation[,x]==1) }, mc.cores = getOption("mc.cores", num_core)
	)# end parallel
	names(annot_list) <- colnames(annotation)
	## inheriting
	object <- new(
		Class = "iDEA",
		gene_id = rownames(annotation),
		annot_id = colnames(annotation),
		summary = summary,
		annotation = annot_list,
		num_gene = num_gene,
		num_core = num_core,
		project = project
	)
	return(object)
}# end function


#' Compute the weight for each gene to construct gene-gene correlation
#'
#' @param counts Raw count data, each column is a cell and each row is a gene.
#' @param normalization Bool variable to denote whether to do normalization before computing the weight, The default is TRUE.
#' @param lib_size Read deepth for each cell, For default, summarize all read counts across whole genes.
#' @param method The method to compute correlation matrix, pearson (default).
#' 
#' @return Return the weight for each gene, p x 1 vector
#'
#' @export
#'
ComputeWeight <- function(counts, normalization=TRUE, lib_size=NULL, method="pearson"){
	## compute r square
	if(normalization){
		data_fpkm <- QQNorm(counts, lib_size=lib_size)
	}else{
		data_fpkm <- counts
	}#end fi
	rm(counts)
	rsq <- cor(t(data_fpkm), method=method)^2
	## alternatively we can remove the jth gene itself
	weight <- apply(rsq, 2, sum)
	## return weight for each gene
	return(1.0/weight)
}# end func


#' Quantile-Quantile normalization for count data
#'
#' @param counts Raw count data, each column is a cell and each row is a gene.
#' @param lib_size Read deepth for each cell, For default, summarize all read counts across whole genes.
#' 
#' @return Return normalized data
#'
#' @export
#'

QQNorm <- function(counts, lib_size=NULL){
	
	## compute read deepth
	if(is.null(lib_size)){
		lib_size <- apply(counts, 2, sum)
	}# end fi
	
	num_gene <- nrow(counts)
	num_cell <- ncol(counts)

	norm_counts <- apply(counts, 1, function(x){
		m <- x/lib_size
		n <- as.matrix(m)
		s <- sample(num_cell, num_cell, replace=F)
		n[s] <- qqnorm(n[s], plot.it=F)$x
		return(n)
	} )
	colnames(norm_counts) <- colnames(counts)
	rownames(norm_counts) <- rownames(counts)
	## return cell length parameter
	return(norm_counts)
}# end func


######################################################
#             Calculate FDR by permuted null         #
######################################################
#' We constructed a permuted null distribution by permuting gene labels. Specifically, we permuted the gene labels of annotation matrix to construced the null pvakue for gene sets. Then we calculate FDR for each gene sets given the permuted null distribution
#'
#' @param object.alt iDEA object under the alternative. iDEA analysis in this object was performed using real gene sets
#' @param object.null iDEA object under the permuted null. iDEA analysis in this object was performed by permuting the gene labels in gene sets.
#' @param numPermute Number of permutation used in object.null. Default is 10.
#'
#' @import stats
#' @import utils
#' @return Returns a data frame which stores the results of iDEA with FDR values for each gene sets
#'
#' @export
#'

iDEA.FDR <- function(object.alt,object.null,numPermute = 10){
    alt = object.alt@gsea$pvalue_louis
    names(alt) = object.alt@gsea$annot_id
    null = object.null@gsea$pvalue_louis
    order_alt = alt[order(alt,decreasing = F)]
    order_null = null[order(null,decreasing = F)]
    total_false <- 0
    fdr_values <-c()
    for (ifdr in 1:length(order_alt)){
        total_false <- sum(order_null < order_alt[ifdr])
        e<-total_false/(numPermute*ifdr)
        fdr_values<-c(fdr_values,e)
    }
    names(fdr_values) = names(alt)
    FDR = data.frame(annot_id = names(fdr_values),fdr = fdr_values)
    df = merge(object.alt@gsea,FDR,by = "annot_id")
    df = df[match(names(alt),df$annot_id),]
    return(df)
}

    
    
    
    
    
    



