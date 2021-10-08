// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(mvtnorm)]]
#include <mvtnormAPI.h>
#include <RcppDist.h>
using namespace Rcpp;




// [[Rcpp::export]]
List DirNegCurvRcpp(arma::vec which_par, arma::vec gr, arma::mat hess) 
{
  int npar = which_par.size();
  arma::vec null=arma::zeros<arma::vec>(2);
  arma::mat Delta_fz=arma::zeros<arma::mat>(npar,npar);
  List eigenvects(npar);
  for (int i=0; i<(npar-1); i++)
  {
    List eigenvects_i(npar);
    for (int j=(i+1); j<npar; j++)
    {
      int i_par=which_par(i)-1;
      int j_par=which_par(j)-1;
      arma::ivec ind = { i_par,j_par };
      arma::uvec indices = arma::conv_to<arma::uvec>::from(ind);

      arma::mat hij=hess(indices,indices);
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, hij);
      double eijmin=eigval(0);
      if (eijmin<0) 
      {
        arma::vec eigvec0=eigvec.col(0);
        arma::vec grij=gr(indices);
        arma::vec grd = grij.t()*eigvec0;
        double grd0 = grd(0);
        arma::vec sign_grd = sign(grd);
        double sign_grd0 = sign_grd(0);
        if (sign_grd0==0) sign_grd0 = 1;
        arma::vec d = -sign_grd0*eigvec0;
        Delta_fz(i,j) = -sign_grd0*grd0+0.5*eijmin;
        eigenvects_i(j)=d;
      }
    }
    eigenvects(i)=eigenvects_i;
  }
  double Delta_fz_min;
  int i_par, j_par;
  arma::vec eigenvector(2);
  Delta_fz_min = Delta_fz.min();
  if (Delta_fz_min<0)
  {
    arma::uword sel1v = Delta_fz.index_min(); // in case of more than one value, returns the index of the first one
    arma::uvec sel1 = ind2sub(size(Delta_fz), sel1v);
    arma::uword i=sel1(0);
    arma::uword j=sel1(1);
    List eigenvector_i = eigenvects(i);
    arma::vec eigenvector_ij=eigenvector_i(j);
    eigenvector=eigenvector_ij;
    i_par=which_par(i);
    j_par=which_par(j);
  }
  
  else
  {
    Delta_fz_min=0;
    i_par=1;
    j_par=1;
    eigenvector=null;
  }
  return List::create(
    _["i_par"] = i_par,
    _["j_par"] = j_par,
    _["Delta_fz_min"] = Delta_fz_min,
    _["eigenvector"] = eigenvector
  );
}


// [[Rcpp::export]]
List DirNetwRcpp(arma::vec gr, arma::mat hess, arma::vec escl) 
{
  int npar = gr.size();
  arma::vec null=arma::zeros<arma::vec>(2);
  arma::mat Delta_fz=arma::mat(npar,npar);
  Delta_fz.fill(9999);
  List direction(npar);
  for (int i=0; i<(npar-1); i++)
  {
    if(all(escl!=i))
    {
      List direction_i(npar);
      for (int j=(i+1); j<npar; j++)
      {
        if(all(escl!=j))
        {
          arma::ivec ind = { i,j };
          arma::uvec indices = arma::conv_to<arma::uvec>::from(ind);
          arma::mat hij=hess(indices,indices);
          arma::vec grij=gr(indices);
          double det_hij=hij(0,0)*hij(1,1)-hij(1,0)*hij(0,1);
          bool notsing= abs(det_hij>0.000001);
          if (notsing)
          {
            arma::mat inv_hij(2,2);
            inv_hij(0,0)=hij(1,1);
            inv_hij(1,1)=hij(0,0);
            inv_hij(0,1)=-hij(0,1);
            inv_hij(1,0)=-hij(1,0);
            inv_hij/=det_hij;
            arma::vec drc = -inv_hij*grij;
            direction_i(j) = drc;
            double fn_drc = pow(accu(pow(drc,2)),0.5);
            arma::vec tmp = grij.t()*drc;
            Delta_fz(i,j) = 0.5*tmp(0) / fn_drc;
          }
          else
          {
            direction_i(j) = null;
          }
        }
        direction(i)=direction_i;
      }
    }
  }
  double Delta_fz_min = Delta_fz.min();
  arma::uword sel1v = Delta_fz.index_min();
  arma::uvec sel1 = ind2sub(size(Delta_fz), sel1v);
  arma::uword i=sel1(0);
  arma::uword j=sel1(1);
  List direction_i = direction(i);
  arma::vec direction_ij=direction_i(j);
  
  return List::create(
    _["i_par"] = i+1,
    _["j_par"] = j+1,
    _["Delta_fz_min"] = Delta_fz_min,
    _["direction"] = direction_ij
  );
}

// [[Rcpp::export]]
arma::mat bi_prob_Rcpp(double tau_i, double tau_j, arma::mat vary_ij)
{
  int n=2;
  int nu=0;
  double lower[2]= {0,1};
  double upper[2]= {tau_i,tau_j};
  int infin[2]= {0,0};
  double corrF[1];
  corrF[0]=vary_ij(0,1);
  double delta[2]= {0,0};
  int maxpts=25000;
  double abseps = 0.001;
  double releps = 0;
  double error = 0;
  double value = 0;
  int inform=0;
  int one=1;
  
  mvtnorm_C_mvtdst(&n, &nu, lower, upper, infin, corrF, delta, 
                   &maxpts, &abseps, &releps, &error, &value, &inform, &one);
  
  arma::mat pi_ij(n,n);
  pi_ij(0,0)=value;
  pi_ij(1,1)=1- ::Rf_pnorm5(tau_i, 0.0, 1.0, 1, 0)-::Rf_pnorm5(tau_j, 0.0, 1.0, 1, 0)+value;
  pi_ij(0,1)=::Rf_pnorm5(tau_i, 0.0, 1.0, 1, 0)-value;
  pi_ij(1,0)=::Rf_pnorm5(tau_j, 0.0, 1.0, 1, 0)-value;
  
  return(pi_ij);
}



// [[Rcpp::export]]
double lik_fa2Rcpp(arma::vec par, List bifreq, int nitems, int D, double eta=0.0, bool fixrotat=false)
{
  arma::mat parmat = arma::conv_to<arma::mat>::from(par);
  parmat.reshape(D+1,nitems);
  arma::mat thresholds=parmat.row(D);
  arma::mat Lambda = parmat.rows(0,D-1);
  if (fixrotat)
  {
    arma::uvec upper_indices = trimatu_ind(size(Lambda),nitems-D+1);
    Lambda(upper_indices)*=0.0;
  }
  arma::mat varf;
  varf = Lambda.t() * Lambda;
  arma::mat vary=varf;
  vary.diag().ones();
  double loglik=0;
  for (int i=0; i<(nitems-1); i++)
  {
    List bifreq_i = bifreq[i];
    for (int j=(i+1); j<nitems; j++) 
    {
      arma::mat bifreq_ij=bifreq_i[j];
      arma::ivec ind = { i,j };
      arma::uvec indices = arma::conv_to<arma::uvec>::from(ind);
      arma::mat vary_ij=vary.submat(indices,indices);
      double tau_i = thresholds(i);
      double tau_j = thresholds(j);
      arma::mat pi_ij=bi_prob_Rcpp(tau_i,tau_j,vary_ij);
      double lij = accu(bifreq_ij % log(pi_ij));
      loglik += lij;
    }
  }
  if (eta>0)
  {
    arma::mat Lambda2 = square(Lambda);
    arma::mat Lambda2sum = sum(Lambda2,1);//rowSums->sum_j
    Lambda2sum += 0.01;
    arma::mat Lambda2sum_sqrt = sqrt(Lambda2sum);
    arma::mat pen = sum(Lambda2sum_sqrt,0);//colSums->sum_d
    double penalty=pen(0,0);
    loglik-=eta*penalty;
  }
  return(-loglik);
}


// [[Rcpp::export]]
arma::mat grad_lik_fa2Rcpp(arma::vec par,List bifreq,int nitems,int D, double eta=0,bool fixrotat=false)
{
  arma::mat parmat = arma::conv_to<arma::mat>::from(par);
  parmat.reshape(D+1,nitems);
  arma::mat thresholds=parmat.row(D);
  arma::mat Lambda = parmat.rows(0,D-1);
  if (fixrotat)
  {
    arma::uvec upper_indices = trimatu_ind(size(Lambda),nitems-D+1);
    Lambda(upper_indices)*=0.0;
  }
  arma::mat varf;
  varf = Lambda.t() * Lambda;
  arma::mat vary=varf;
  vary.diag().ones();

  arma::mat der1=arma::zeros<arma::mat>(D+1,nitems);
  arma::mat der_rho=arma::zeros<arma::mat>(nitems,nitems);
  for (int i=0; i<(nitems-1); i++)
  {
    List bifreq_i = bifreq[i];
    double tau_i = thresholds(i);
    double dnorm_tau_i = ::Rf_dnorm4(tau_i,0.0,1.0,0);
    for (int j=(i+1); j<nitems; j++)
    {
      arma::mat bifreq_ij=bifreq_i[j];
      arma::ivec ind = { i,j };
      arma::uvec indices = arma::conv_to<arma::uvec>::from(ind);
      arma::mat vary_ij=vary.submat(indices,indices);
      double rho_ij = vary_ij(0,1);
      double tau_j = thresholds(j);
      arma::mat pi_ij=bi_prob_Rcpp(tau_i,tau_j,vary_ij);

      arma::mat der_pi_ij(2,2);

      // der of pi_ij with respect to tau_i
      der_pi_ij(0,0) = dnorm_tau_i * ::Rf_pnorm5((tau_j-rho_ij*tau_i)/pow(1-pow(rho_ij,2),0.5), 0.0, 1.0, 1, 0);
      der_pi_ij(0,1) = dnorm_tau_i - der_pi_ij(0,0);
      der_pi_ij(1,0) = - der_pi_ij(0,0);
      der_pi_ij(1,1) = - dnorm_tau_i + der_pi_ij(0,0);

      der1(D,i) += accu(bifreq_ij / pi_ij % der_pi_ij);

      // der of pi_ij with respect to tau_j
      double dnorm_tau_j = ::Rf_dnorm4(tau_j,0.0,1.0,0);

      der_pi_ij(0,0) = dnorm_tau_j * ::Rf_pnorm5((tau_i-rho_ij*tau_j)/pow(1-pow(rho_ij,2),0.5), 0.0, 1.0, 1, 0);
      der_pi_ij(0,1) =  - der_pi_ij(0,0);
      der_pi_ij(1,0) = dnorm_tau_j - der_pi_ij(0,0);
      der_pi_ij(1,1) = - dnorm_tau_j + der_pi_ij(0,0);

      der1(D,j) += accu(bifreq_ij / pi_ij % der_pi_ij);

      // der of p_ij with respect to rho

      arma::mat tau_ij = thresholds.cols(indices);
      arma::vec m=arma::zeros<arma::vec>(2);
      arma::vec density_tau_i_tau_j = dmvnorm(tau_ij,m,vary_ij);
      der_pi_ij(0,0) = density_tau_i_tau_j(0);
      der_pi_ij(0,1) = -density_tau_i_tau_j(0);
      der_pi_ij(1,0) = -density_tau_i_tau_j(0);
      der_pi_ij(1,1) = density_tau_i_tau_j(0);

      der_rho(i,j) += accu(bifreq_ij / pi_ij % der_pi_ij);
    }
  }
  der_rho += der_rho.t();
  der1.rows(0,D-1) = Lambda * der_rho;
  if (eta>0)
  {
    // first derivatives
    arma::mat Lambda2 = square(Lambda);
    arma::mat Lambda2sum = sum(Lambda2,1);//rowSums
    Lambda2sum += 0.01;
    arma::mat Lambda2sum_sqrt = sqrt(Lambda2sum);
    arma::mat penaltyder=Lambda;
    for (int d=0; d<D; d++)
      penaltyder.row(d)/=Lambda2sum_sqrt(d);
    der1.rows(0,D-1) -= eta*penaltyder;
  }
  if (fixrotat)
  {
    arma::uvec upper_indices = trimatu_ind(size(der1),nitems-D+1);
    der1(upper_indices)*=0.0;
  }
  
  return(-der1);
}





// [[Rcpp::export]]
List der_lik_fa2Rcpp(arma::vec par,List bifreq,int nitems,int D, double eta=0.0)
{
  arma::mat parmat = arma::conv_to<arma::mat>::from(par);
  parmat.reshape(D+1,nitems);
  arma::mat thresholds=parmat.row(D);
  arma::mat Lambda = parmat.rows(0,D-1).t();
  arma::mat varf;
  varf = Lambda * Lambda.t();
  arma::mat vary=varf;
  vary.diag().ones();
  
  arma::mat der2=arma::zeros<arma::mat>((D+1)*nitems,(D+1)*nitems);
  arma::mat der1=arma::zeros<arma::mat>((D+1)*nitems,1);
  arma::mat der1_lambda=arma::zeros<arma::mat>(D*nitems,1);
  arma::mat der1_tau=arma::zeros<arma::mat>(nitems,1);
  arma::mat der2_tau=arma::zeros<arma::mat>(nitems,nitems);
  arma::mat der2_lambda=arma::zeros<arma::mat>(D*nitems,D*nitems);
  arma::mat der_rho=arma::zeros<arma::mat>(nitems,nitems);
  arma::mat der2_rho=arma::zeros<arma::mat>(nitems,nitems);
  arma::mat der2_lambda_tau=arma::zeros<arma::mat>(nitems,D*nitems);
  arma::mat desn(nitems*D,nitems*D);
  desn.eye();

  for (int i=0; i<(nitems-1); i++)
  {
    List bifreq_i = bifreq[i];
    double tau_i = thresholds(i);
    double dnorm_tau_i = ::Rf_dnorm4(tau_i,0.0,1.0,0);
    for (int j=(i+1); j<nitems; j++) 
    {
      arma::mat bifreq_ij=bifreq_i[j];
      arma::ivec ind = { i,j };
      arma::uvec indices = arma::conv_to<arma::uvec>::from(ind);
      arma::mat vary_ij=vary.submat(indices,indices);
      double rho_ij = vary_ij(0,1);
      double tau_j = thresholds(j);
      arma::mat pi_ij=bi_prob_Rcpp(tau_i,tau_j,vary_ij);
      
      // der of pi_ij with respect to tau_i
      arma::mat der_pi_ij_tau_i(2,2);
      double fn_rho = pow(1-pow(rho_ij,2),0.5);
      double val1 = (tau_j-rho_ij*tau_i)/fn_rho;
      der_pi_ij_tau_i(0,0) = dnorm_tau_i * ::Rf_pnorm5(val1, 0.0, 1.0, 1, 0);
      der_pi_ij_tau_i(0,1) = dnorm_tau_i - der_pi_ij_tau_i(0,0);
      der_pi_ij_tau_i(1,0) = - der_pi_ij_tau_i(0,0);
      der_pi_ij_tau_i(1,1) = - dnorm_tau_i + der_pi_ij_tau_i(0,0);
      der1_tau(i) += accu(bifreq_ij / pi_ij % der_pi_ij_tau_i);
      
      // second der of pi_ij with respect to tau_i^2
      arma::mat der2_pi_ij_tau_i(2,2);
      der2_pi_ij_tau_i(0,0) = -tau_i * dnorm_tau_i* ::Rf_pnorm5(val1, 0.0, 1.0, 1, 0) - dnorm_tau_i * ::Rf_dnorm4(val1,0.0,1.0,0) * rho_ij/fn_rho;
      der2_pi_ij_tau_i(0,1) = -tau_i*dnorm_tau_i - der2_pi_ij_tau_i(0,0);
      der2_pi_ij_tau_i(1,0) = - der2_pi_ij_tau_i(0,0);
      der2_pi_ij_tau_i(1,1) = tau_i * dnorm_tau_i + der2_pi_ij_tau_i(0,0);
      
      arma::mat tmp = - pow(der_pi_ij_tau_i,2) / pow(pi_ij,2) + der2_pi_ij_tau_i / pi_ij;
      der2_tau(i,i) += accu(bifreq_ij % tmp);

      // der of pi_ij with respect to tau_j
      arma::mat der_pi_ij_tau_j(2,2);
      double dnorm_tau_j = ::Rf_dnorm4(tau_j,0.0,1.0,0);
      double val2 = (tau_i-rho_ij*tau_j)/fn_rho;
      
      der_pi_ij_tau_j(0,0) = dnorm_tau_j * ::Rf_pnorm5(val2, 0.0, 1.0, 1, 0);
      der_pi_ij_tau_j(0,1) =  - der_pi_ij_tau_j(0,0);
      der_pi_ij_tau_j(1,0) = dnorm_tau_j - der_pi_ij_tau_j(0,0);
      der_pi_ij_tau_j(1,1) = - dnorm_tau_j + der_pi_ij_tau_j(0,0);

      der1_tau(j) += accu(bifreq_ij / pi_ij % der_pi_ij_tau_j);
      
      // second der of pi_ij with respect to tau_j^2
      arma::mat der2_pi_ij_tau_j(2,2);
      der2_pi_ij_tau_j(0,0) = -tau_j * dnorm_tau_j* ::Rf_pnorm5(val2, 0.0, 1.0, 1, 0) - dnorm_tau_j * ::Rf_dnorm4(val2,0.0,1.0,0) * rho_ij/fn_rho;
      der2_pi_ij_tau_j(0,1) = - der2_pi_ij_tau_j(0,0);
      der2_pi_ij_tau_j(1,0) = -tau_j*dnorm_tau_j - der2_pi_ij_tau_j(0,0);
      der2_pi_ij_tau_j(1,1) = tau_j * dnorm_tau_j + der2_pi_ij_tau_j(0,0);
      
      tmp = - pow(der_pi_ij_tau_j,2) / pow(pi_ij,2) + der2_pi_ij_tau_j / pi_ij;
      der2_tau(j,j) += accu(bifreq_ij % tmp);
      
      // second der of pi_ij with respect to tau_i and tau_j
      arma::mat der2_pi_ij_tau_i_tau_j(2,2);
      der2_pi_ij_tau_i_tau_j(0,0) = dnorm_tau_i * ::Rf_dnorm4(val1,0.0,1.0,0) /fn_rho;
      der2_pi_ij_tau_i_tau_j(0,1) = - der2_pi_ij_tau_i_tau_j(0,0);
      der2_pi_ij_tau_i_tau_j(1,0) = - der2_pi_ij_tau_i_tau_j(0,0);
      der2_pi_ij_tau_i_tau_j(1,1) = der2_pi_ij_tau_i_tau_j(0,0);
      tmp = - der_pi_ij_tau_i % der_pi_ij_tau_j / pow(pi_ij,2) + der2_pi_ij_tau_i_tau_j / pi_ij;
      der2_tau(i,j) += accu(bifreq_ij % tmp);

      // der of p_ij with respect to rho
      arma::mat der_pi_ij_rho(2,2);
      arma::mat tau_ij = thresholds.cols(indices);
      arma::vec m=arma::zeros<arma::vec>(2);
      arma::vec density_tau_i_tau_j = dmvnorm(tau_ij,m,vary_ij);
      der_pi_ij_rho(0,0) = density_tau_i_tau_j(0);
      der_pi_ij_rho(0,1) = -density_tau_i_tau_j(0);
      der_pi_ij_rho(1,0) = -density_tau_i_tau_j(0);
      der_pi_ij_rho(1,1) = density_tau_i_tau_j(0);
      der_rho(i,j) += accu(bifreq_ij / pi_ij % der_pi_ij_rho);
      
      // second der of pi_ij with respect to rho^2
      arma::mat der2_pi_ij_rho(2,2);
      double z = pow(tau_i,2) - 2*rho_ij * tau_i * tau_j + pow(tau_j,2);
      double fn1_rho = 1-pow(rho_ij,2);
      der2_pi_ij_rho(0,0) = density_tau_i_tau_j(0)/fn1_rho*rho_ij+density_tau_i_tau_j(0)*(tau_i*tau_j*fn1_rho-z*rho_ij)/pow(fn1_rho,2);
      der2_pi_ij_rho(0,1) = -der2_pi_ij_rho(0,0);
      der2_pi_ij_rho(1,0) = -der2_pi_ij_rho(0,0);
      der2_pi_ij_rho(1,1) = der2_pi_ij_rho(0,0);
      tmp = - pow(der_pi_ij_rho,2) / pow(pi_ij,2) + der2_pi_ij_rho / pi_ij;
      der2_rho(i,j) += accu(bifreq_ij % tmp);
      
      // second der of loglik with respect to lambda
      arma::mat desn_i = desn.rows(i*D,i*D+D-1);
      arma::mat desn_j = desn.rows(j*D,j*D+D-1);
      arma::mat der_rho_lambda = Lambda.row(j) * desn_i + Lambda.row(i) * desn_j;
      der1_lambda += der_rho(i,j) * der_rho_lambda.t();
      der2_lambda += der_rho(i,j) * (desn_j.t() * desn_i+desn_i.t() * desn_j);
      der2_lambda += der2_rho(i,j) * der_rho_lambda.t() * der_rho_lambda;

      // second der of pi_ij with respect to rho and tau_i
      arma::mat der2_pi_ij_rho_tau(2,2);
      der2_pi_ij_rho_tau(0,0) = density_tau_i_tau_j(0)*(-2*tau_i+2*rho_ij*tau_j)/(2*fn1_rho);
      der2_pi_ij_rho_tau(0,1) = -der2_pi_ij_rho_tau(0,0);
      der2_pi_ij_rho_tau(1,0) = -der2_pi_ij_rho_tau(0,0);
      der2_pi_ij_rho_tau(1,1) = der2_pi_ij_rho_tau(0,0);
      
      // second der of loglik with respect to rho and tau_i
      tmp = -der_pi_ij_tau_i % der_pi_ij_rho / pow(pi_ij,2) + der2_pi_ij_rho_tau / pi_ij;
      
      double der2_rho_tau = accu(bifreq_ij % tmp);
      // second der of loglik with respect to lambda and tau_i
      der2_lambda_tau.row(i) += der2_rho_tau * der_rho_lambda;
      
      // second der of pi_ij with respect to rho and tau_j
      der2_pi_ij_rho_tau(0,0) = density_tau_i_tau_j(0)*(-2*tau_j+2*rho_ij*tau_i)/(2*fn1_rho);
      der2_pi_ij_rho_tau(0,1) = -der2_pi_ij_rho_tau(0,0);
      der2_pi_ij_rho_tau(1,0) = -der2_pi_ij_rho_tau(0,0);
      der2_pi_ij_rho_tau(1,1) = der2_pi_ij_rho_tau(0,0);
      // second der of loglik with respect to rho and tau_j
      tmp = -der_pi_ij_tau_j % der_pi_ij_rho / pow(pi_ij,2) + der2_pi_ij_rho_tau / pi_ij;
      der2_rho_tau = accu(bifreq_ij % tmp);
      // second der of loglik with respect to lambda and tau_i
      der2_lambda_tau.row(j) += der2_rho_tau * der_rho_lambda;
    }
  }
  arma::uvec intercepts = arma::regspace<arma::uvec>(D,D+1,(D+1)*nitems-1);
  arma::uvec nointer = arma::regspace<arma::uvec>(0,(D+1)*nitems-1);
  nointer.shed_rows(intercepts);
  der1.rows(intercepts) = der1_tau;
  der1.rows(nointer) = der1_lambda;
  arma::mat der2_tau_tmp = der2_tau;
  der2_tau_tmp.diag().zeros();
  der2_tau += der2_tau_tmp.t();
  der2(intercepts,intercepts) = der2_tau;
  der2(intercepts,nointer) = der2_lambda_tau;
  der2(nointer,intercepts) = der2_lambda_tau.t();
  der2(nointer,nointer) = der2_lambda;
  
  if (eta>0)
  {
    // first derivatives
    arma::mat Lambda2 = square(Lambda);
    arma::mat Lambda2sum = sum(Lambda2+0.001,0);//colSums
    arma::mat Lambda2sum_sqrt = sqrt(Lambda2sum);
    arma::mat penaltyder=Lambda;
    for (int d=0; d<D; d++)
      penaltyder.col(d)/=Lambda2sum_sqrt(d);
    arma::mat penaltydervec = reshape(penaltyder.t(),D*nitems,1);
    der1.rows(nointer) -= eta*penaltydervec;
    
    // second derivatives
    arma::mat Lambda2sum_m1p5 = pow(Lambda2sum,-1.5);
    arma::mat Lambdavec = reshape(Lambda.t(),D*nitems,1);
    arma::mat LambdaLambda = Lambdavec * Lambdavec.t();
    arma::mat diagD(D, D, arma::fill::eye);
    arma::mat desnmat = repmat(diagD,nitems,nitems);
    LambdaLambda %= desnmat;
    arma::mat Lambda2sum_m1p5_desn = repmat(diagmat(Lambda2sum_m1p5),nitems,nitems);
    arma::mat penaltyder2 = - LambdaLambda % Lambda2sum_m1p5_desn;
    penaltyder2.diag() += 1/ repmat(Lambda2sum_sqrt,1,nitems);
    der2(nointer,nointer) -= eta * penaltyder2;
  }

  return List::create(
    _["grad"] = -der1,
    _["hess"] = -der2
  );
}






// [[Rcpp::export]]
arma::mat NormalOgiveLikRcpp(arma::vec par, arma::mat data, arma::mat nodes, 
                             arma::vec weights, arma::mat numpatt, int D=1, bool fixrotat=false)
{
  int nitems = data.n_cols;
  int nq = nodes.n_rows;
  
  arma::mat parmat = arma::conv_to<arma::mat>::from(par);
  parmat.reshape(D+1,nitems);
  
  if (fixrotat)
  {
    arma::uvec upper_indices = trimatu_ind(size(parmat),nitems-D+1);
    parmat(upper_indices)*=0.0;
  }

  arma::mat onesmat = arma::ones<arma::mat>(nq,1);
  arma::mat uno_nodes = join_rows(nodes, onesmat);//cbind
  arma::mat prob = uno_nodes*parmat;
  prob.for_each( [](arma::mat::elem_type& val) { val = ::Rf_pnorm5(val, 0.0, 1.0, 1, 0); } );
  
  arma::mat probt = prob.t();
  arma::uvec ids = find(probt >= 1-1.4e-07);
  probt.elem(ids).fill(1-1.4e-07);
  ids = find(probt <= 1.4e-07);
  probt.elem(ids).fill(1.4e-07);
  
  arma::mat prodj= exp(data*log(probt)+(1-data)*log(1-probt));
  
  arma::mat logsumq = log(prodj * weights) ;
  arma::mat mlik =  - sum(logsumq%numpatt);
  
  return mlik;
}















