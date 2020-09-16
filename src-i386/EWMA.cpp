// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; using namespace arma;

double MaxTlower(double tauMax, double sigma, double lambda, int n) {
  double w = 2 * tauMax;
  return ceil(log((12 * n * w * (2 - lambda) * sigma * sqrt(lambda / n / (2 - lambda)) - w) / 36 / lambda / pow(sigma, 2)) / 2 / log(1 - lambda));
}


double thetaFunc(double lambda, int h){
  return 1 - pow(1 - lambda, 2 * h);
}

double tauFunc(double L, double sigma, double lambda, int h, int t) {
  double theta = thetaFunc(lambda, h);
  double res = 0;
  if (h == 0) {
    res = L * sigma * sqrt(lambda / (2 - lambda)) / (2 * t + 1);
  } else if (h > 0) {
    res = L * sigma * sqrt(lambda / (2 - lambda) * theta) / (2 * t + 1);
  }
  return res;
}

double upperFunc(double mu, double lambda, double tau, int k, int l) {
  return mu + (2 * l - (1 - lambda) * 2 * k + 1) / lambda * tau;
}

double lowerFunc(double mu, double lambda, double tau, int k, int l) {
  return mu + (2 * l - (1 - lambda) * 2 * k - 1) / lambda * tau;
}

double muFunc(double U, int m) {
  return R::qnorm(U, 0.0, 1.0, 1, 0) / sqrt(m); 
}

double sigmaFunc(double V, double nu, double ubc) {
  return sqrt(R::qchisq(V, nu, 1, 0) / nu) / ubc;
}

double qFunc(double lower, double upper) {
  return R::pnorm(upper, 0.0, 1.0, 1, 0) - R::pnorm(lower, 0.0, 1.0, 1, 0);
}

arma::mat QFunc(double U, double V, double L, double lambda, double nu, double ubc, int h, int t, int m) {
  
  arma::mat Q(2 * t + 1, 2 * t + 1);
  double mu = muFunc(U, m);
  double sigma = sigmaFunc(V, nu, ubc);
  double tau = tauFunc(L, sigma, lambda, h, t);
  double lower = 0;
  double upper = 0;
  int k = 0;
  int l = 0;
  
  int i;
  int j;
  
  for (i = 0; i < 2 * t + 1; i++) {
    k = i - t;
    for (j = 0; j < 2 * t + 1; j++) {
      l = j - t;
      
      lower = lowerFunc(mu, lambda, tau, k, l);
      upper = upperFunc(mu, lambda, tau, k, l);
      Q(i, j) = qFunc(lower, upper);
    }
  }
  
  return Q;
  
}

int checkInvert(arma::mat Q) {
  arma::colvec QrowSums = arma::sum(Q, 1);
  int m = Q.n_cols;
  int Flg1 = 0;
  int Flg2 = 0;
  int out = 0;
  
  int i;
  
  for (i = 0; i < m; i++) {
    if (0 < QrowSums(i) and QrowSums(i) <= 1)  Flg1++;
    if (0 < QrowSums(i) or QrowSums(i) < 1)  Flg2++;
  }
  
  if (Flg1 == m or Flg2 > 0) out = 1;
  
  return out;
}

List QMaxFunc(double U, double V, double L, double lambda, double nu, double ubc, int tmin, int tmax, int m) {
  
  arma::mat Q = QFunc(U, V, L, lambda, nu, ubc, 0, tmax, m);
  int invert = 0;
  int t_ = 0;
  int i_ = 0;
  
  do {
    i_++;
    
    invert = checkInvert(Q);
    
    if (invert == 0) {
      t_ = tmax - i_;
      Q = QFunc(U, V, L, lambda, nu, ubc, 0, t_, m); 
    }
    
    if (t_ == tmin) {
      stop("Matrix Q is not invertible");
    }
    
  } while (invert == 0);
  
  return List::create(_["Q"] = Q, _["t"] = (Q.n_cols - 1) / 2);
  
}

arma::mat XiFunc(int t, int loc) {
  int g = 2 * t + 1;
  arma::mat out = zeros(1, g);
  
  out(0, loc) = 1;
  
  return out;
  
}

double CARLZeroFunc(double U, double V, double L, double lambda, double nu, double ubc, int hmax, int t, int m, arma::mat Xi, arma::mat QMax) {
  
  int g = 2 * t + 1;
  arma::mat Qc = eye(g, g);
  arma::mat Nmat = eye(g, g);
  arma::mat o1 = ones(g, 1);
  arma::mat ey = eye(g, g);
  arma::mat out1(1, 1);
  double out;
  
  int hh;
  
  for (hh = 0; hh < hmax - 1; hh++) {
    
    Qc = Qc * QFunc(U, V, L, lambda, nu, ubc, hh + 1, t, m);
    
    Nmat = Nmat + Qc;
 
  }
  
  Nmat = Nmat + Qc * arma::inv(ey - QMax);
  
  out1 = Xi * Nmat * o1;
  
  out = out1(0, 0);
    
  return out;
  
}

double CARLSteadyFunc(arma::mat Xi, arma::mat QMax, int t) {
  
  int g = 2 * t + 1;
  arma::mat o1 = ones(g, 1);
  arma::mat ey = eye(g, g);
  arma::mat out1(1, 1);
  double out;

  out1 = Xi * arma::inv(ey - QMax) * o1;
  
  out = out1(0, 0);
  
  return out;
  
}

// [[Rcpp::export]]
double integrandSteady(double U, double V, double L, double lambda, int mm, double nu, double ubc, int tmin, int tmax) {
  
    List tmp = QMaxFunc(U, V, L, lambda, nu, ubc, tmin, tmax, mm);
    
    arma::mat Qmax = tmp["Q"];
    int ttmp = tmp["t"];
    
    arma::mat xi = XiFunc(ttmp, ttmp);
    
    double out = CARLSteadyFunc(xi, Qmax, ttmp);
	
	  return out;
      
}

// [[Rcpp::export]]
int qintegrandSteady(double q, double U, double V, double L, double lambda, int mm, double nu, double ubc, int tmin, int tmax) {
  
  double CARL = integrandSteady(U, V, L, lambda, mm, nu, ubc, tmin, tmax);
  int out = 0;
  
  if (CARL <= q) {
    out = 1;
  }
  
  return out;
  
}

// [[Rcpp::export]]
double integrandZero(double U, double V, double L, double lambda, int mm, double nu, int nn, double ubc, int tmin, int tmax) {
    
    List tmp = QMaxFunc(U, V, L, lambda, nu, ubc, tmin, tmax, mm);
    
    arma::mat Qmax = tmp["Q"];
    int ttmp = tmp["t"];
    
    arma::mat xi = XiFunc(ttmp, ttmp);
    
    double sigma = sigmaFunc(V, nu, ubc);
    double tau = tauFunc(L, sigma, lambda, 0, ttmp); 
    int hmax = MaxTlower(tau, sigma, lambda, nn); 
    
    double out = CARLZeroFunc(U, V, L, lambda, nu, ubc, hmax, ttmp, mm, xi, Qmax);
	
	  return out;
    
}

// [[Rcpp::export]]
int qintegrandZero(double q, double U, double V, double L, double lambda, int mm, double nu, int nn, double ubc, int tmin, int tmax) {

  double CARL = integrandZero(U, V, L, lambda, mm, nu, nn, ubc, tmin, tmax);
  int out = 0;
  
  if (CARL <= q) {
    out = 1;
  }
  
  return out;
  
}
