#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
Rcpp::List updateMuCpp(const arma::vec& idx,const arma::mat& s,const arma::mat& mu,
	const arma::mat& mu_star,const arma::vec& sigma,const arma::vec& a,const arma::vec& rho, 
	const arma::vec& nu,arma::vec& int_,arma::vec& int_star,arma::vec& Xbeta_mut,
	arma::vec& Xbeta_mut_star,arma::vec& delta_d,arma::mat& delta_d_star,arma::vec& delta_t,arma::vec& phi){
			
	Rcpp::RNGScope scope;
	
	int N=delta_t.n_elem;	
	int n=idx.n_elem;
	arma::vec out(N,fill::zeros);
	int j;
	arma::vec mh_star_mu(1);
	arma::vec mh_0_mu(1);
	Rcpp::NumericVector mh(1);
	Rcpp::NumericVector rn(1);
	arma::rowvec dif(2);	
	arma::vec kernel_1(1);
	arma::vec kernel_2(1);
	arma::vec dmvt(1);
	arma::vec dmvt_star(1);
	arma::mat Sigma(2,2);
	arma::mat Sigma_tilda(2,2);

	for(int i=0;i<n;i++){
		j = idx(i)-1;
		Sigma(0,0) = sigma(j)*sigma(j);
		Sigma(0,1) = Sigma(1,0) = sigma(j)*sigma(j)*sqrt(a(j))*rho(j);
		Sigma(1,1) = sigma(j)*sigma(j)*a(j);
		Sigma_tilda = Sigma;
		Sigma_tilda(0,1) = Sigma_tilda(1,0) = -Sigma(1,0);
		dif = s.row(j)-mu.row(j);
		kernel_1 = pow(arma::as_scalar(1+1/nu(j)*sum((dif*inv(Sigma))%dif)),-(nu(j)+2)/2);
		kernel_2 = pow(arma::as_scalar(1+1/nu(j)*sum((dif*inv(Sigma_tilda))%dif)),-(nu(j)+2)/2);
		dmvt = arma::as_scalar(log(kernel_1+kernel_2));
		dif = s.row(j)-mu_star.row(j);
		kernel_1 = pow(arma::as_scalar(1+1/nu(j)*sum((dif*inv(Sigma))%dif)),-(nu(j)+2)/2);
		kernel_2 = pow(arma::as_scalar(1+1/nu(j)*sum((dif*inv(Sigma_tilda))%dif)),-(nu(j)+2)/2);
		dmvt_star = arma::as_scalar(log(kernel_1+kernel_2));
		mh_star_mu = dmvt_star + Xbeta_mut_star(j) - (delta_d_star(j,0)/(delta_t(j)*phi)) 
				-(delta_d_star(j,1)/(delta_t(j+1)*phi)) - log(int_star(j+1));
	   	mh_0_mu = dmvt + Xbeta_mut(j) - (delta_d(j)/(delta_t(j)*phi))  
	    		-(delta_d(j+1)/(delta_t(j+1)*phi)) - log(int_(j+1));
		mh = exp(mh_star_mu-mh_0_mu);
		rn = runif(1,0,1);
		bool res = Rcpp::is_true(any(mh>rn));
		if(res){
			out(j) = 1;
			delta_d(j) = delta_d_star(j,0);
			delta_d(j+1) = delta_d_star(j,1);
		}
	}
	return Rcpp::List::create(Rcpp::Named("out")=out,Rcpp::Named("delta.d")=delta_d);
}
