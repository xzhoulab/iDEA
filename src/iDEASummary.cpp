#include <iostream>
#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>
#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;
using namespace Rcpp;

#define ARMA_DONT_PRINT_ERRORS

// A template version of the above
void swap(double &a, double &b)
{
	double t = a;
	a = b;
	b = t;
} // end func
//**********************************************************************//
//                   Fitting Logistic Regression                        //
//**********************************************************************//
void Mstep(mat A, vec POST_gamma, vec &tau, vec &var_coef, mat &J, vec sigma2_e, vec beta, vec beta2, double &sigma2_beta, double &sigma2_beta2)
{
	double a_0beta = 3.0;
	double b_0beta = 20.0;
	//========================================================
	// *update sigma2_b, modidifed by sun, 2018-10-29 06:53:12
	sigma2_beta = (sum((beta % beta % POST_gamma) / (2.0 * sigma2_e)) + b_0beta) / (sum(POST_gamma) / 2.0 + a_0beta + 1.0);
	vec tau_old = tau;
	size_t true_iter = 1;
	//mat J = zeros<mat>( tau.n_elem, tau.n_elem );
	size_t nr_iter = 10; // maximum iteration of newton raphson
	while (true_iter)
	{
		vec pi = exp(A * tau_old) / (1 + exp(A * tau_old));
		J = A.t() * diagmat(pi % (1 - pi)) * A;
		tau = tau + inv_sympd(J) * (A.t() * (POST_gamma - pi));
		if (arma::norm(tau_old - tau) < 1e-03 || (++true_iter) > nr_iter)
		{
			break;
		}
		else
		{
			tau_old = tau;
		}
	}
	//exit(1);
	// the variance of annotation
	mat inv_J = inv_sympd(J);
	var_coef = inv_J.diag();

} // end function

//**********************************************************************//
//                   Log Likelihood Function                            //
//**********************************************************************//
double logL(vec SSR, double sigma2b, double sigma2e, vec beta, vec gamma, vec pi)
{

	double log_density = 0.5 * dot(SSR, SSR) / sigma2e;
	log_density += 0.5 * SSR.n_elem * log(sigma2e);
	log_density += 0.5 * sum(gamma % (log(sigma2b) + log(sigma2e) + beta % beta / (sigma2b * sigma2e) - 2 * log(pi)));
	log_density -= sum((1 - gamma) % log(1 - pi));
	log_density += (0.1 + 1) * log(sigma2b) + 0.1 / sigma2b;
	return -log_density;
}

double logL(vec SSR, double sigma2b, vec sigma2e, vec beta, vec gamma, vec pi)
{
	// here sigma2e is a vector, input by users
	double log_density = 0.5 * dot(SSR, diagmat(sigma2e) * SSR);
	log_density += 0.5 * sum(log(sigma2e));

	log_density += 0.5 * sum(gamma % (log(sigma2b) + log(sigma2e) + beta % beta / (sigma2b * sigma2e) - 2 * log(pi)));
	log_density -= sum((1 - gamma) % log(1 - pi));
	log_density += (0.1 + 1) * log(sigma2b) + 0.1 / sigma2b;
	return -log_density;
}

double logL(vec SSR, double sigma2b, double sigma2b2, vec sigma2e, vec beta, vec beta2, vec gamma, vec pi)
{
	// here sigma2e is a vector, input by users
	double log_density = 0.5 * dot(SSR, diagmat(sigma2e) * SSR);
	log_density += 0.5 * sum(log(sigma2e));

	log_density += 0.5 * sum(gamma % (log(sigma2b) + log(sigma2e) + beta % beta / (sigma2b * sigma2e) - 2 * log(pi)));
	log_density += 0.5 * sum((1 - gamma) % (log(sigma2b2) + log(sigma2e) + beta2 % beta2 / (sigma2b2 * sigma2e) - 2 * log(1 - pi)));
	log_density += (0.1 + 1) * log(sigma2b) + 0.1 / sigma2b;
	return -log_density;
}

double MlogL(mat A, vec tau, vec prior_pi)
{

	double a = as_scalar(dot(prior_pi, A * tau));
	double b = as_scalar(sum(log(1 + exp(A * tau))));

	return a - b;
} // end func

void CheckGamma(arma::vec &gamma, int min_gene)
{
	// at least one gene is selected
	if (sum(gamma) < min_gene)
	{
		arma::uvec index_vec(gamma.n_elem);
		index_vec.fill(1);
		index_vec = cumsum(index_vec) - 1;
		arma::uvec shut_index_vec = shuffle(index_vec);
		gamma.zeros();
		// change elements of A greater than 0 to 1
		for (size_t i = 0; i <= min_gene; i++)
		{
			gamma(shut_index_vec(i)) = 1;
		} // end for
	}
	else if (sum(gamma) > 0.6 * gamma.n_elem)
	{
		// at most half genes are selected
		// cout<<"index_vec length"<<endl;
		arma::uvec index_vec = find(gamma > 0);
		arma::uvec shut_index_vec = shuffle(index_vec);
		gamma.zeros();
		// change elements of A greater than 0 to 1
		for (size_t i = 0; i < ceil(0.5 * gamma.n_elem); i++)
		{
			gamma(shut_index_vec(i)) = 1;
		} // end for
	}	 // end fi
} // end func

//*************************************************//
//                    MCMC Step                    //
//*************************************************//
void Estep(vec y, vec sigma2_e, vec &beta, vec &POST_gamma, vec &POST_beta, vec &POST_beta2, double &sigma2_beta, double &sigma2_beta2, vec pi, double &log_l_hat, int mcmc_iter, int min_degene)
{
	int num_gene = sigma2_e.n_elem;
	//int num_cvt = W.n_cols;

	// *set the initial values of parameters
	// *burn in times
	int s_step = ceil(0.2 * mcmc_iter);
	int w_step = mcmc_iter - s_step;

	vec mean_beta = zeros<vec>(num_gene);
	vec var_beta = zeros<vec>(num_gene);
	mean_beta = beta;

	vec mean_beta2 = zeros<vec>(num_gene);
	vec var_beta2 = zeros<vec>(num_gene);

	vec Ebeta = beta;
	// *initial value for sampling sigma2_e
	vec XEbeta = Ebeta;
	double WEalpha = 0.0;
	vec SSR = y - XEbeta - WEalpha;

	// *the length of beta and gamma sampling equals to the number of genes
	//vec beta_s =zeros<vec>(num_gene);
	vec beta_s = beta;
	vec beta_s2 = beta;

	vec gamma_s = zeros<vec>(num_gene);
	for (size_t i = 0; i < num_gene; i++)
	{
		if (pi(i) > as_scalar(randu(1)))
		{
			gamma_s(i) = 1;
		} // *pi is the probability of gamma = 1
		else
		{
			gamma_s(i) = 0;
		}
	} // end
	vec POST_Ebeta = zeros<vec>(num_gene);

	// probability of gamma, related to pi
	mat log_pi = zeros<mat>(num_gene, 2);
	vec pi_b1 = zeros<vec>(num_gene);

	double xHSSR = 0.0, xHx = 0.0;
	//double xHx2 = 0.0;

	// to avoid label switching issue
	if (sigma2_beta < sigma2_beta2)
	{
		swap(sigma2_beta, sigma2_beta2);
	} // end fi
	// *the main loop of MCMC
	for (size_t iter = 0; iter < (s_step + w_step); iter++)
	{

		// *update the beta in turn,beta from spike and slab
		// Ebeta = beta_s%gamma_s; // pair-wise multipl
		Ebeta = beta_s % gamma_s + beta_s2 % (1.0 - gamma_s); // pair-wise multipl
		if (iter >= w_step)
		{
			POST_Ebeta += Ebeta;
			POST_beta += beta_s % gamma_s;
			POST_beta2 += beta_s2 % (1 - gamma_s);
		} // end fi

		//XEbeta = X*Ebeta;
		XEbeta = Ebeta;
		SSR.zeros();
		SSR = y - WEalpha;
		// n x 1

		// *for each gene
		for (int i = 0; i < num_gene; i++)
		{
			xHSSR = y(i); // *here, Xbeta has removed ith row

			// *for case, gamma = 1
			// *modified by sun 2018.7.26.18:32:51
			//xHx = dot(X.col(i), X.col(i) ) + 1/sigma2_beta;
			xHx = 1.0 + 1.0 / sigma2_beta;
			mean_beta(i) = xHSSR / xHx;
			//var_beta[i] = sigma2_e/xHx;
			var_beta(i) = sigma2_e(i) / xHx;
			// *sampling beta with mean and variance
			beta_s(i) = as_scalar(randn(1)) * sqrt(var_beta(i)) + mean_beta(i);
			log_pi(i, 1) = mean_beta(i) * mean_beta(i) / (2 * var_beta(i)) + log(sqrt(var_beta(i))) - log(sqrt(sigma2_beta)) + log(pi(i));
			log_pi(i, 0) = log(1.0 - pi(i));
			//beta_s2(i) = as_scalar(randn(1))*sqrt(var_beta2(i)) + mean_beta2(i);
			beta_s2(i) = 0;
			//cout<<"pi(1) = " << exp(log_pi(i, 1)) <<"pi(0) = "<<exp(log_pi(i, 0)) <<endl;
			// *new pi
			log_pi.row(i) -= max(log_pi.row(i));

			// * the variable used for sampling
			pi_b1(i) = exp(log_pi(i, 1)) / (exp(log_pi(i, 0)) + exp(log_pi(i, 1)));
			//cout<<"pi_b1(i) = " << pi_b1(i) <<endl;

			if (pi_b1(i) > 0.9999)
			{
				pi_b1(i) = 0.9999;
			} // end fi
			if (pi_b1(i) < 0.0001)
			{
				pi_b1(i) = 0.0001;
			} // end fi
			//cout<<"Estep::pi_b1 = "<<pi_b1[i]<<endl;

			// Bernoulli distribution
			if (pi_b1(i) > as_scalar(randu(1)))
			{
				gamma_s(i) = 1; // bool variable
			}
			else
			{
				gamma_s(i) = 0;
			} // end if

			Ebeta(i) = beta_s(i) * gamma_s(i) + beta_s2(i) * (1.0 - gamma_s(i));
		} // end for each gene

		if ((sum(gamma_s) < min_degene) | (sum(gamma_s) > 0.6 * gamma_s.n_elem))
		{
			CheckGamma(gamma_s, min_degene);
		}

		XEbeta = Ebeta;
		if (iter >= w_step)
		{
			POST_gamma += gamma_s;
		} // end fi

		//================================================================
		// *sampling alpha for covariates
		//WEalpha.zeros();
		WEalpha = 0.0;
	} // end MCMC iteration

	POST_gamma /= (double)s_step;
	POST_beta /= (double)s_step;
	POST_beta2 /= (double)s_step;
	POST_Ebeta /= (double)s_step;
	XEbeta.zeros();
	XEbeta = POST_Ebeta;
	//WEalpha = W*(POST_alpha);

	SSR.zeros();
	SSR = y - XEbeta - WEalpha;
	log_l_hat = logL(SSR, sigma2_beta, sigma2_beta2, sigma2_e, POST_beta, POST_beta2, POST_gamma, pi);

	//beta = POST_beta;
	beta = XEbeta;
} // *end function

//*******************************************************************//
//                Expectation Maximazation                           //
//*******************************************************************//
//' EM-MCMC function based on summary statistics, effect size and variance of effect size
//' @param yin Estimated beta
//' @param varbetain Variance
//' @param Ain Annotation
//' @param betain Initial beta
//' @param tauin Initial annotation coefficient
//' @param em_iterin EM iteration
//' @param mcmc_iterin MCMC iteration
//' @param min_degenein Minimum DE genes
//'
//' @return A list
//'
//' @export
// [[Rcpp::export]]
SEXP EMMCMCStepSummary(SEXP yin, SEXP varbetain, SEXP Ain, SEXP betain, SEXP tauin, SEXP em_iterin, SEXP mcmc_iterin, SEXP min_degenein)
{
	try
	{						  // *EM Algorithm
		vec y = as<vec>(yin); // *dim = num_sample x 1
		mat A = as<mat>(Ain); // *dim = num_gene x num_annot

		vec var_beta = as<vec>(varbetain);
		vec beta = as<vec>(betain); // *dim = num_gene x 1, initial value can be calculated by linear model
		vec tau = as<vec>(tauin);   // *annotation coefficient

		int em_iter = Rcpp::as<int>(em_iterin);
		int mcmc_iter = Rcpp::as<int>(mcmc_iterin);
		int min_degene = Rcpp::as<int>(min_degenein);

		int num_gene = A.n_rows;
		vec POST_beta = zeros<vec>(num_gene);
		vec POST_beta2 = zeros<vec>(num_gene);

		vec POST_gamma;
		double sigma2_beta = 2.0;
		double sigma2_beta2 = 0.05;
		// *compute pi
		vec prior_pi = exp(A * tau) / (1.0 + exp(A * tau));
		vec tau_new = tau;
		double E_flogL_old = 999999.0, E_flogL = 0.0, E_ElogL = 0.0;
		vec tau_var = zeros<vec>(tau.n_elem);
		// variance matrix for annotation, 2 x 2
		mat info_mat = zeros<mat>(tau.n_elem, tau.n_elem);
		//mat hist_tau = zeros<mat>(tau.n_elem, em_iter);
		for (size_t iter = 0; iter < em_iter; iter++)
		{
			// * only store the last iteration results
			POST_gamma = zeros<vec>(var_beta.n_elem);

			//======================================
			// *MCMC step to estimate the gamma
			//double tstart1 = clock();
			Estep(y, var_beta, beta, POST_gamma, POST_beta, POST_beta2, sigma2_beta, sigma2_beta2, prior_pi, E_ElogL, mcmc_iter, min_degene);
			//double time_mcmc = (clock() - tstart1) / (double(CLOCKS_PER_SEC));
			Mstep(A, POST_gamma, tau_new, tau_var, info_mat, var_beta, POST_beta, POST_beta2, sigma2_beta, sigma2_beta2);
			//double time_nr = (clock() - tstart2) / (double(CLOCKS_PER_SEC));

			// *compute the pi
			prior_pi = exp(A * tau_new) / (1.0 + exp(A * tau_new));

			// *compute the full likelihood here
			double E_MlogL = MlogL(A, tau_new, prior_pi);
			E_flogL = E_MlogL + E_ElogL;
			//hist_tau.col(iter) = tau_new;

			if ((abs(E_flogL - E_flogL_old) < 1e-02) || max(abs(tau_new - tau) / abs(tau_new + tau)) < 1e-02)
			{
				break;
			}
			else
			{
				E_flogL_old = E_flogL;
				tau = tau_new;
			} // end fi
		}	 // end for EM iteration

		return List::create(Named("beta") = beta, Named("sigma2_beta1") = sigma2_beta,
							Named("sigma2_beta2") = sigma2_beta2, Named("annot_coef") = tau_new,
							Named("annot_var") = tau_var, Named("pip") = POST_gamma,
							Named("sigma2_e") = var_beta, Named("info_mat") = info_mat);
	}
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end func

//**********************************************************************//
//                   Fitting Logistic Regression                        //
//**********************************************************************//
void MstepWeight(mat A, vec weight, vec POST_gamma, vec &tau, vec &var_coef, mat &J, vec sigma2_e, vec beta, vec beta2, double &sigma2_beta, double &sigma2_beta2)
{
	// *tau is coefficent of annotation
	double a_0beta = 3.0;
	double b_0beta = 20.0;
	// double b_0beta = 10*(a_0beta-1)/num_gene;
	//========================================================
	// *update sigma2_b, modidifed by sun, 2018-10-29 06:53:12
	sigma2_beta = (sum((beta % beta % POST_gamma % weight) / (2.0 * sigma2_e)) + b_0beta) / (sum(POST_gamma % weight) / 2.0 + a_0beta + 1.0);
	//sigma2_beta=10;

	vec tau_old = tau;
	size_t true_iter = 1;
	//mat J = zeros<mat>( tau.n_elem, tau.n_elem );
	int nr_iter = 10; // maximum iteration of newton raphson
	while (true_iter)
	{
		vec pi = exp(A * tau_old) / (1 + exp(A * tau_old));

		J = A.t() * diagmat(pi % (1 - pi)) * A;

		tau = tau + inv_sympd(J) * (A.t() * (POST_gamma - pi));

		//cout<<"Mstep:: tau = "<<tau<<endl;
		if (arma::norm(tau_old - tau) < 1e-03 || (++true_iter) > nr_iter)
		{
			break;
		}
		else
		{
			tau_old = tau;
		}
	}
	//exit(1);
	// the variance of annotation
	mat inv_J = inv_sympd(J);
	var_coef = inv_J.diag();

} // end function

//*************************************************//
//                    MCMC Step                    //
//*************************************************//
void EstepWeight(vec y, vec sigma2_e, vec weight, vec &beta, vec &POST_gamma, vec &POST_beta, vec &POST_beta2, double &sigma2_beta, double &sigma2_beta2, vec pi, double &log_l_hat, int mcmc_iter, int min_degene)
{
	int num_gene = sigma2_e.n_elem;
	//int num_cvt = W.n_cols;

	// *set the initial values of parameters
	// *burn in times
	int s_step = ceil(0.2 * mcmc_iter);
	int w_step = mcmc_iter - s_step;

	//==================
	// *for beta, sigma2_beta: inverse gamma distribution(a_0beta,b_0beta)
	double a_0beta = 3.0;
	double b_0beta = 20.0;
	// double b_0beta = 10*(a_0beta-1)/num_gene;
	double a_0beta2 = 2000.0;
	double b_0beta2 = 2.0;

	vec mean_beta = zeros<vec>(num_gene);
	vec var_beta = zeros<vec>(num_gene);
	mean_beta = beta;

	vec mean_beta2 = zeros<vec>(num_gene);
	vec var_beta2 = zeros<vec>(num_gene);

	vec Ebeta = beta;
	// *initial value for sampling sigma2_e
	vec XEbeta = Ebeta;
	double WEalpha = 0.0;
	vec SSR = y - XEbeta - WEalpha;

	// *the length of beta and gamma sampling equals to the number of genes
	vec beta_s = beta;
	vec beta_s2 = beta;

	vec gamma_s = zeros<vec>(num_gene);
	for (size_t i = 0; i < num_gene; i++)
	{
		if (pi(i) > as_scalar(randu(1)))
		{
			gamma_s(i) = 1;
		} // *pi is the probability of gamma = 1
		else
		{
			gamma_s(i) = 0;
		}
	} // end

	vec POST_Ebeta = zeros<vec>(num_gene);

	// *modified by sun
	//double POST_sigma2_e = 0.0;
	double POST_sigma2_b = 0.0;
	double POST_sigma2_b2 = 0.0;

	// probability of gamma, related to pi
	mat log_pi = zeros<mat>(num_gene, 2);
	vec pi_b1 = zeros<vec>(num_gene);

	double xHSSR = 0.0, xHx = 0.0, xHx0 = 0.0, wHSSR = 0.0, wHw = 0.0;
	double xHx2 = 0.0;

	// to avoid label switching issue
	if (sigma2_beta < sigma2_beta2)
	{
		swap(sigma2_beta, sigma2_beta2);
	} // end fi
	// *the main loop of MCMC
	for (size_t iter = 0; iter < (s_step + w_step); iter++)
	{

		Ebeta = beta_s % gamma_s + beta_s2 % (1.0 - gamma_s); // pair-wise multipl
		if (iter >= w_step)
		{
			POST_Ebeta += Ebeta;
			POST_beta += beta_s % gamma_s;
			POST_beta2 += beta_s2 % (1 - gamma_s);
		} // end fi

		//XEbeta = X*Ebeta;
		XEbeta = Ebeta;
		SSR.zeros();
		SSR = y - WEalpha;
		// n x 1

		// *for each gene
		for (int i = 0; i < num_gene; i++)
		{

			xHSSR = y(i); // *here, Xbeta has removed ith row

			// *for case, gamma = 1
			// *modified by sun 2018.7.26.18:32:51
			xHx = 1.0 + 1.0 / sigma2_beta;
			mean_beta(i) = xHSSR / xHx;
			// modified by sun, 2019-2-5 10:20:14
			var_beta(i) = weight(i) * sigma2_e(i) / xHx;
			// *sampling beta with mean and variance
			beta_s(i) = as_scalar(randn(1)) * sqrt(var_beta(i)) + mean_beta(i);
			// modified by sun, 2019-2-5 10:20:55
			log_pi(i, 1) = mean_beta(i) * mean_beta(i) / (2 * var_beta(i)) + log(sqrt(var_beta(i))) - weight(i) * log(sqrt(sigma2_beta)) + weight(i) * log(pi(i));
			log_pi(i, 0) = log(1.0 - pi(i));

			// second component
			beta_s2(i) = 0;
			// *new pi
			log_pi.row(i) -= max(log_pi.row(i));

			// * the variable used for sampling
			pi_b1(i) = exp(log_pi(i, 1)) / (exp(log_pi(i, 0)) + exp(log_pi(i, 1)));
			//cout<<"pi_b1(i) = " << pi_b1(i) <<endl;

			if (pi_b1(i) > 0.9999)
			{
				pi_b1(i) = 0.9999;
			}
			if (pi_b1(i) < 0.0001)
			{
				pi_b1(i) = 0.0001;
			}
			//cout<<"Estep::pi_b1 = "<<pi_b1[i]<<endl;

			// *Bernoulli distribution
			if (pi_b1(i) > as_scalar(randu(1)))
			{
				gamma_s(i) = 1; // bool variable
			}
			else
			{
				gamma_s(i) = 0;
			} // *end if

			Ebeta(i) = beta_s(i) * gamma_s(i) + beta_s2(i) * (1.0 - gamma_s(i));
		} // *end for each gene
		if (sum(gamma_s) < min_degene | sum(gamma_s) > 0.6 * gamma_s.n_elem)
		{
			CheckGamma(gamma_s, min_degene);
		}
		XEbeta = Ebeta;
		if (iter >= w_step)
		{
			POST_gamma += gamma_s;
		}

		WEalpha = 0.0;
	} // *end MCMC iteration

	POST_gamma /= (double)s_step;
	POST_beta /= (double)s_step;
	POST_beta2 /= (double)s_step;
	POST_Ebeta /= (double)s_step;

	XEbeta.zeros();
	XEbeta = POST_Ebeta;
	//WEalpha = W*(POST_alpha);

	SSR.zeros();
	SSR = y - XEbeta - WEalpha;
	log_l_hat = logL(SSR, sigma2_beta, sigma2_beta2, sigma2_e, POST_beta, POST_beta2, POST_gamma, pi);
	beta = XEbeta;
} // *end function

//*******************************************************************//
//                Expectation Maximazation                           //
//*******************************************************************//
//' EM-MCMC with weight function based on summary statistics, effect size and variance of effect size
//' @param yin Estimated beta
//' @param varbetain Variance
//' @param weightin Weigth vector
//' @param Ain Annotation
//' @param betain Initial beta
//' @param tauin Initial annotation coefficient
//' @param em_iterin EM iteration
//' @param mcmc_iterin MCMC iteration
//' @param min_degenein Minimum DE genes
//'
//' @return A list
//'
//' @export
// [[Rcpp::export]]
SEXP EMMCMCStepSummaryWeight(SEXP yin, SEXP varbetain, SEXP weightin, SEXP Ain, SEXP betain, SEXP tauin, SEXP em_iterin, SEXP mcmc_iterin, SEXP min_degenein)
{
	try
	{						  // *EM Algorithm
		vec y = as<vec>(yin); // *dim = num_sample x 1
		mat A = as<mat>(Ain); // *dim = num_gene x num_annot
		vec weight = as<vec>(weightin);

		vec var_beta = as<vec>(varbetain);
		vec beta = as<vec>(betain); // *dim = num_gene x 1, initial value can be calculated by linear model
		//vec alpha = as<vec>(alphain); // *fix effect coefficient
		vec tau = as<vec>(tauin); // *annotation coefficient

		int em_iter = Rcpp::as<int>(em_iterin);
		int mcmc_iter = Rcpp::as<int>(mcmc_iterin);
		int min_degene = Rcpp::as<int>(min_degenein);

		int num_gene = A.n_rows;
		vec POST_beta = zeros<vec>(num_gene);
		vec POST_beta2 = zeros<vec>(num_gene);

		vec POST_gamma;
		double sigma2_beta = 2.0;
		double sigma2_beta2 = 0.05;
		// *compute pi
		vec prior_pi = exp(A * tau) / (1.0 + exp(A * tau));
		vec tau_new = tau;
		double E_flogL_old = 999999.0, E_flogL = 0.0, E_ElogL = 0.0;
		vec tau_var = zeros<vec>(tau.n_elem);
		// variance matrix for annotation, 2 x 2
		mat info_mat = zeros<mat>(tau.n_elem, tau.n_elem);
		//mat hist_tau = zeros<mat>(tau.n_elem, em_iter);
		// *the main loop of EM algorithm
		/////////////////////////////////////////
		//==================
		// *for beta, sigma2_beta: inverse gamma distribution(a_0beta,b_0beta)
		double a_0beta = 3.0;
		double b_0beta = 20.0;
		// double b_0beta = 10*(a_0beta-1)/num_gene;
		double a_0beta2 = 2000.0;
		double b_0beta2 = 2.0;
		// a_beta, b_beta, a_beta2, b_beta2
		double a_beta = sum(prior_pi) / 2.0 + a_0beta;
		double b_beta = sum(beta / (2.0 * var_beta)) + b_0beta;

		double a_beta2 = sum(1.0 - prior_pi) / 2.0 + a_0beta2;
		double b_beta2 = sum(beta / (2.0 * var_beta)) + b_0beta2;
		/////////////////////////////////////////////////////////////////
		for (size_t iter = 0; iter < em_iter; iter++)
		{
			// * only store the last iteration results
			POST_gamma = zeros<vec>(var_beta.n_elem);

			//======================================
			// *MCMC step to estimate the gamma
			//cout<<"E-step: "<< iter+1 <<endl;
			double tstart1 = clock();
			//Estep(y, se_beta, W, beta, alpha, POST_gamma, prior_pi, E_ElogL, mcmc_iter);
			EstepWeight(y, var_beta, beta, weight, POST_gamma, POST_beta, POST_beta2, sigma2_beta, sigma2_beta2, prior_pi, E_ElogL, mcmc_iter, min_degene);
			double time_mcmc = (clock() - tstart1) / (double(CLOCKS_PER_SEC));
			//cout<<"LogL (mcmc) = "<< E_ElogL << endl;
			//=======================================
			// *logist regression step to estimate annotation coefficients, tau
			//cout<<"M-step: "<< iter+1 <<endl;
			double tstart2 = clock();
			//cout<<"POST_gamma = "<<POST_gamma<<endl;
			MstepWeight(A, weight, POST_gamma, tau_new, tau_var, info_mat, var_beta, POST_beta, POST_beta2, sigma2_beta, sigma2_beta2);
			double time_nr = (clock() - tstart2) / (double(CLOCKS_PER_SEC));
			//cout<<"est annot coef, tau = "<< tau_new << endl;

			// *compute the pi
			//cout<<"EMstep:: compute the pi"<<endl;
			prior_pi = exp(A * tau_new) / (1.0 + exp(A * tau_new));

			// *compute the full likelihood here
			double E_MlogL = MlogL(A, tau_new, prior_pi);
			E_flogL = E_MlogL + E_ElogL;
			//hist_tau.col(iter) = tau_new;

			if ((abs(E_flogL - E_flogL_old) < 1e-02) | max(abs(tau_new - tau) / abs(tau_new + tau)) < 1e-02)
			{
				break;
			}
			else
			{
				E_flogL_old = E_flogL;
				tau = tau_new;
			}
			//if(max(abs(tau_new - tau)/abs(tau_new+tau))<1e-03 ){ break; }
			//else{ tau = tau_new; }

		} // end for EM iteration

		return List::create(Named("beta") = beta, Named("sigma2_beta1") = sigma2_beta,
							Named("sigma2_beta2") = sigma2_beta2, Named("annot_coef") = tau_new,
							Named("annot_var") = tau_var, Named("pip") = POST_gamma,
							Named("sigma2_e") = var_beta, Named("info_mat") = info_mat);
	}
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end func with weight

/////////////////////////////////////////////////////////////////////////////////////////
//                             CODE END HERE                                           //
/////////////////////////////////////////////////////////////////////////////////////////
