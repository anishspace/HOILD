// [[Rcpp::depends("RcppArmadillo", "RcppDist")]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends("RcppProgress")]]

#define BOOST_DISABLE_ASSERTS true

// #define ZTOL (DOUBLE_EPS*10.0)
#define ZTOL (2.220446e-16 * 10.0)

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/math/tr1.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "sitmo.hpp"
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

// FCN prototypes
double samplepg(double);
double exprnd(double);
double tinvgauss(double, double);
double truncgamma();
double randinvg(double);
double aterm(int, double, double);

// set seed
// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::export]]
NumericVector rcpp_pgdraw(NumericVector b, NumericVector c)
{
  int m = b.size();
  int n = c.size();
  NumericVector y(n);
  
  // Setup
  int i, j, bi = 1;
  if (m == 1)
  {
    bi = b[0];
  }
  
  // Sample
  for (i = 0; i < n; i++)
  {
    if (m > 1)
    {
      bi = b[i];
    }
    
    // Sample
    y[i] = 0;
    for (j = 0; j < (int)bi; j++)
    {
      y[i] += samplepg(c[i]);
    }
  }
  
  return y;
}


// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/
  // WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg(double z)
{
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
double t = MATH_2_PI;

// Compute p, q and the ratio q / (q + p)
// (derived from scratch; derivation is not in the original paper)
double K = z*z/2.0 + MATH_PI2/8.0;
double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
double logK = (double)std::log(K);
double Kt = K * t;
double w = (double)std::sqrt(MATH_PI_2);

double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
double ratio = 1.0 / (1.0 + p_over_q); 

double u, X;

// Main sampling loop; page 130 of the Windle PhD thesis
while(1) 
{
  // Step 1: Sample X ? g(x|z)
  u = R::runif(0.0,1.0);
  if(u < ratio) {
    // truncated exponential
    X = t + exprnd(1.0)/K;
  }
  else {
    // truncated Inverse Gaussian
    X = tinvgauss(z, t);
  }
  
  // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for 
  // an odd n or U > Sn(X|z) for an even n
  int i = 1;
  double Sn = aterm(0, X, t);
  double U = R::runif(0.0,1.0) * Sn;
  int asgn = -1;
  bool even = false;
  
  while(1) 
  {
    Sn = Sn + asgn * aterm(i, X, t);
    
    // Accept if n is odd
    if(!even && (U <= Sn)) {
      X = X * 0.25;
      return X;
    }
    
    // Return to step 1 if n is even
    if(even && (U > Sn)) {
      break;
    }
    
    even = !even;
    asgn = -asgn;
    i++;
  }
}
return X;
}

// Generate exponential distribution random variates
double exprnd(double mu)
{
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
  // Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm(int n, double x, double t)
{
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- 
                                                         (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }    
  return (double)exp(f);
}

// Generate inverse gaussian random variates
double randinvg(double mu)
{
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {    
    out = mu*mu / out; 
  }    
  return out;
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables 
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma()
{
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;  
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss(double z, double t)
{
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma 
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }  
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }    
  return X;
}








/****************************************************************************************
  * *********************** Generalized Inverse Gaussian Sampling ************************
  ****************************************************************************************/
  
  double _gig_mode(double lambda, double omega) {
    /*---------------------------------------------------------------------------*/
      /* Compute mode of GIG distribution.                                         */
      /*                                                                           */
      /* Parameters:                                                               */
      /*   lambda .. parameter for distribution                                    */
      /*   omega ... parameter for distribution                                    */
      /*                                                                           */
      /* Return:                                                                   */
      /*   mode                                                                    */
      /*---------------------------------------------------------------------------*/
      
      if (lambda >= 1.) {
        /* mode of fgig(x) */
          return (sqrt((lambda-1.)*(lambda-1.) + omega*omega)+(lambda-1.))/omega;
      } else {
        /* 0 <= lambda < 1: use mode of f(1/x) */
          return omega / (sqrt((1.-lambda)*(1.-lambda) + omega*omega)+(1.-lambda));
      }
  } /* end of _gig_mode() */
  
  /*---------------------------------------------------------------------------*/
  
  void _rgig_ROU_noshift (arma::rowvec &res, int n, double lambda, double lambda_old, 
                          double omega, double alpha) {
    /*---------------------------------------------------------------------------*/
      /* Tpye 1:                                                                   */
      /* Ratio-of-uniforms without shift.                                          */
      /*   Dagpunar (1988), Sect.~4.6.2                                            */
      /*   Lehner (1989)                                                           */
      /*---------------------------------------------------------------------------*/
      
      double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
        double ym, um;     /* location of maximum of x*sqrt(f(x)); umax of MBR */
        double s, t;       /* auxiliary variables */
        double U, V, X;    /* random variables */
        
        int i;             /* loop variable (number of generated random variables) */
        int count = 0;     /* counter for total number of iterations */
          
          /* -- Setup -------------------------------------------------------------- */
          
          /* shortcuts */
          t = 0.5 * (lambda-1.);
          s = 0.25 * omega;
          
          /* mode = location of maximum of sqrt(f(x)) */
            xm = _gig_mode(lambda, omega);
            
            /* normalization constant: c = log(sqrt(f(xm))) */
              nc = t*log(xm) - s*(xm + 1./xm);
              
              /* location of maximum of x*sqrt(f(x)):           */
                /* we need the positive root of                   */
                /*    omega/2*y^2 - (lambda+1)*y - omega/2 = 0    */
                  ym = ((lambda+1.) + sqrt((lambda+1.)*(lambda+1.) + omega*omega))/omega;
                  
                  /* boundaries of minmal bounding rectangle:                   */
                    /* we us the "normalized" density f(x) / f(xm). hence         */
                    /* upper boundary: vmax = 1.                                  */
                      /* left hand boundary: umin = 0.                              */
                        /* right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(xm)) */
                          um = exp(0.5*(lambda+1.)*log(ym) - s*(ym + 1./ym) - nc);
                          
                          /* -- Generate sample ---------------------------------------------------- */
                            
                            for (i=0; i<n; i++) {
                              do {
                                ++count;
                                U = um * unif_rand();        /* U(0,umax) */
                                  V = unif_rand();             /* U(0,vmax) */
                                    X = U/V;
                              }                              /* Acceptance/Rejection */
                                while (((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));
                              
                              /* store random point */
                                res(i) = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
                            }
                          
                          /* -- End ---------------------------------------------------------------- */
                            
                            return;
  } /* end of _rgig_ROU_noshift() */
  
  
  /*---------------------------------------------------------------------------*/
  
  void _rgig_newapproach1 (arma::rowvec &res, int n, double lambda, 
                           double lambda_old, double omega, double alpha) {
    /*---------------------------------------------------------------------------*/
      /* Type 4:                                                                   */
      /* New approach, constant hat in log-concave part.                           */
      /* Draw sample from GIG distribution.                                        */
      /*                                                                           */
      /* Case: 0 < lambda < 1, 0 < omega < 1                                       */
      /*                                                                           */
      /* Parameters:                                                               */
      /*   n ....... sample size (positive integer)                                */
      /*   lambda .. parameter for distribution                                    */
      /*   omega ... parameter for distribution                                    */
      /*                                                                           */
      /* Return:                                                                   */
      /*   random sample of size 'n'                                               */
      /*---------------------------------------------------------------------------*/
      
      /* parameters for hat function */
      double A[3], Atot;  /* area below hat */
      double k0;          /* maximum of PDF */
      double k1, k2;      /* multiplicative constant */
      
      double xm;          /* location of mode */
      double x0;          /* splitting point T-concave / T-convex */
      double a;           /* auxiliary variable */
      
      double U, V, X;     /* random numbers */
      double hx;          /* hat at X */
      
      int i;              /* loop variable (number of generated random variables) */
      int count = 0;      /* counter for total number of iterations */
        
        /* -- Check arguments ---------------------------------------------------- */
        if (lambda >= 1. || omega >1.) {
          stop ("error: invalid parameters");
        }
      
      /* -- Setup -------------------------------------------------------------- */
        
        /* mode = location of maximum of sqrt(f(x)) */
          xm = _gig_mode(lambda, omega);
          
          /* splitting point */
            x0 = omega/(1.-lambda);
            
            /* domain [0, x_0] */
              k0 = exp((lambda-1.)*log(xm) - 0.5*omega*(xm + 1./xm));     /* = f(xm) */
                A[0] = k0 * x0;
                
                /* domain [x_0, Infinity] */
                  if (x0 >= 2./omega) {
                    k1 = 0.;
                    A[1] = 0.;
                    k2 = pow(x0, lambda-1.);
                    A[2] = k2 * 2. * exp(-omega*x0/2.)/omega;
                  } else {
                    /* domain [x_0, 2/omega] */
                      k1 = exp(-omega);
                      A[1] = (lambda == 0.) 
                      ? k1 * log(2./(omega*omega))
                      : k1 / lambda * ( pow(2./omega, lambda) - pow(x0, lambda) );
                      
                      /* domain [2/omega, Infinity] */
                        k2 = pow(2/omega, lambda-1.);
                        A[2] = k2 * 2 * exp(-1.)/omega;
                  }
                
                /* total area */
                  Atot = A[0] + A[1] + A[2];
                  
                  /* -- Generate sample ---------------------------------------------------- */
                    
                    for (i=0; i<n; i++) {
                      do {
                        ++count;
                        
                        /* get uniform random number */
                          V = Atot * unif_rand();
                          
                          do {
                            
                            /* domain [0, x_0] */
                              if (V <= A[0]) {
                                X = x0 * V / A[0];
                                hx = k0;
                                break;
                              }
                            
                            /* domain [x_0, 2/omega] */
                              V -= A[0];
                              if (V <= A[1]) {
                                if (lambda == 0.) {
                                  X = omega * exp(exp(omega)*V);
                                  hx = k1 / X;
                                }
                                else {
                                  X = pow(pow(x0, lambda) + (lambda / k1 * V), 1./lambda);
                                  hx = k1 * pow(X, lambda-1.);
                                }
                                break;
                              }
                              
                              /* domain [max(x0,2/omega), Infinity] */
                                V -= A[1];
                                a = (x0 > 2./omega) ? x0 : 2./omega;
                                X = -2./omega * log(exp(-omega/2. * a) - omega/(2.*k2) * V);
                                hx = k2 * exp(-omega/2. * X);
                                break;
                                
                          } while(0);
                          
                          /* accept or reject */
                            U = unif_rand() * hx;
                            
                            if (log(U) <= (lambda-1.) * log(X) - omega/2. * (X+1./X)) {
                              /* store random point */
                                res(i) = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
                                break;
                            }
                      } while(1);
                      
                    }
                  
                  /* -- End ---------------------------------------------------------------- */
                    
                    return;
  } /* end of _rgig_newapproach1() */
  
  /*---------------------------------------------------------------------------*/
  
  void _rgig_ROU_shift_alt (arma::rowvec &res, int n, double lambda, 
                            double lambda_old, double omega, double alpha) {
    /*---------------------------------------------------------------------------*/
      /* Type 8:                                                                   */
      /* Ratio-of-uniforms with shift by 'mode', alternative implementation.       */
      /*   Dagpunar (1989)                                                         */
      /*   Lehner (1989)                                                           */
      /*---------------------------------------------------------------------------*/
      
      double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
        double s, t;       /* auxiliary variables */
        double U, V, X;    /* random variables */
        
        int i;             /* loop variable (number of generated random variables) */
        int count = 0;     /* counter for total number of iterations */
          
          double a, b, c;    /* coefficent of cubic */
          double p, q;       /* coefficents of depressed cubic */
          double fi, fak;    /* auxiliary results for Cardano's rule */
          
          double y1, y2;     /* roots of (1/x)*sqrt(f((1/x)+m)) */
          
          double uplus, uminus;  /* maximum and minimum of x*sqrt(f(x+m)) */
          
          /* -- Setup -------------------------------------------------------------- */
          
          /* shortcuts */
          t = 0.5 * (lambda-1.);
          s = 0.25 * omega;
          
          /* mode = location of maximum of sqrt(f(x)) */
          xm = _gig_mode(lambda, omega);
          
          /* normalization constant: c = log(sqrt(f(xm))) */
          nc = t*log(xm) - s*(xm + 1./xm);
          
          /* location of minimum and maximum of (1/x)*sqrt(f(1/x+m)):  */
          
          /* compute coeffients of cubic equation y^3+a*y^2+b*y+c=0 */
          a = -(2.*(lambda+1.)/omega + xm);       /* < 0 */
          b = (2.*(lambda-1.)*xm/omega - 1.);
          c = xm;
          
          /* we need the roots in (0,xm) and (xm,inf) */
          
          /* substitute y=z-a/3 for depressed cubic equation z^3+p*z+q=0 */
          p = b - a*a/3.;
          q = (2.*a*a*a)/27. - (a*b)/3. + c;
          
          /* use Cardano's rule */
          fi = acos(-q/(2.*sqrt(-(p*p*p)/27.)));
          fak = 2.*sqrt(-p/3.);
          y1 = fak * cos(fi/3.) - a/3.;
          y2 = fak * cos(fi/3. + 4./3.*M_PI) - a/3.;
          
          /* boundaries of minmal bounding rectangle:                  */
            /* we us the "normalized" density f(x) / f(xm). hence        */
            /* upper boundary: vmax = 1.                                 */
              /* left hand boundary: uminus = (y2-xm) * sqrt(f(y2)) / sqrt(f(xm)) */
                /* right hand boundary: uplus = (y1-xm) * sqrt(f(y1)) / sqrt(f(xm)) */
                  uplus  = (y1-xm) * exp(t*log(y1) - s*(y1 + 1./y1) - nc);
                  uminus = (y2-xm) * exp(t*log(y2) - s*(y2 + 1./y2) - nc);
                  
                  /* -- Generate sample ---------------------------------------------------- */
                    
                    for (i=0; i<n; i++) {
                      do {
                        ++count;
                        U = uminus + unif_rand() * (uplus - uminus);    /* U(u-,u+)  */
                          V = unif_rand();                                /* U(0,vmax) */
                            X = U/V + xm;
                      }                                         /* Acceptance/Rejection */
                        while ((X <= 0.) || ((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));
                      
                      /* store random point */
                        res(i) = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
                    }
                  
                  /* -- End ---------------------------------------------------------------- */
                    
                    return;
  } /* end of _rgig_ROU_shift_alt() */
  
  /*---------------------------------------------------------------------------*/
  // [[Rcpp::export]]
double bessel_k_nuasympt (double x, double nu, int islog, int expon_scaled) {
  /*---------------------------------------------------------------------------*/
    /* Asymptotic expansion of Bessel K_nu(x) function                           */
    /* when BOTH  nu and x  are large.                                           */
    /*                                                                           */
    /* parameters:                                                               */
    /*   x            ... argument for K_nu()                                    */
    /*   nu           ... order or Bessel function                               */
    /*   islog        ... return logarithm of result TRUE and result when FALSE  */
    /*   expon_scaled ... return exp(-x)*K_nu(x) when TRUE and K_nu(x) when FALSE*/
    /*                                                                           */
    /*---------------------------------------------------------------------------*/
    /*                                                                           */
    /* references:                                                               */
    /* ##  Abramowitz & Stegun , p.378, __ 9.7.8. __                             */
    /*                                                                           */
    /* ## K_nu(nu * z) ~ sqrt(pi/(2*nu)) * exp(-nu*eta)/(1+z^2)^(1/4)            */
    /* ##                   * {1 - u_1(t)/nu + u_2(t)/nu^2 - ... }               */
    /*                                                                           */
    /* ## where   t = 1 / sqrt(1 + z^2),                                         */
    /* ##       eta = sqrt(1 + z^2) + log(z / (1 + sqrt(1+z^2)))                 */
    /* ##                                                                        */
    /* ## and u_k(t)  from  p.366  __ 9.3.9 __                                   */
    /*                                                                           */
    /* ## u0(t) = 1                                                              */
    /* ## u1(t) = (3*t - 5*t^3)/24                                               */
    /* ## u2(t) = (81*t^2 - 462*t^4 + 385*t^6)/1152                              */
    /* ## ...                                                                    */
    /*                                                                           */
    /* ## with recursion  9.3.10    for  k = 0, 1, .... :                        */
    /* ##                                                                        */
    /* ## u_{k+1}(t) = t^2/2 * (1 - t^2) * u'_k(t) +                             */
    /* ##            1/8  \int_0^t (1 - 5*s^2)* u_k(s) ds                        */
    /*---------------------------------------------------------------------------*/
    /*                                                                           */
    /* Original implementation in R code (R package "Bessel" v. 0.5-3) by        */
    /*   Martin Maechler, Date: 23 Nov 2009, 13:39                               */
    /*                                                                           */
    /* Translated into C code by Kemal Dingic, Oct. 2011.                        */
    /*                                                                           */
    /* Modified by Josef Leydold on Tue Nov  1 13:22:09 CET 2011                 */
    /*                                                                           */
    /*---------------------------------------------------------------------------*/
    
    #define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
    double z;                   /* rescaled argument for K_nu() */
    double sz, t, t2, eta;      /* auxiliary variables */
    double d, u1t,u2t,u3t,u4t;  /* (auxiliary) results for Debye polynomials */
    double res;                 /* value of log(K_nu(x)) [= result] */
    
    /* rescale: we comute K_nu(z * nu) */
    z = x / nu;
    
    /* auxiliary variables */
      sz = hypot(1,z);   /* = sqrt(1+z^2) */
        t = 1. / sz;
        t2 = t*t;
        
        eta = (expon_scaled) ? (1./(z + sz)) : sz;
        eta += log(z) - log1p(sz);                  /* = log(z/(1+sz)) */
          
          /* evaluate Debye polynomials u_j(t) */
          u1t = (t * (3. - 5.*t2))/24.;
          u2t = t2 * (81. + t2*(-462. + t2 * 385.))/1152.;
          u3t = t*t2 * (30375. + t2 * (-369603. + t2 * (765765. - t2 * 425425.)))/414720.;
          u4t = t2*t2 * (4465125. 
                         + t2 * (-94121676.
                                 + t2 * (349922430. 
                                         + t2 * (-446185740. 
                                                 + t2 * 185910725.)))) / 39813120.;
          d = (-u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;
          
          /* log(K_nu(x)) */
            res = log(1.+d) - nu*eta - 0.5*(log(2.*nu*sz) - M_LNPI);
            
            return (islog ? res : exp(res));
} /* end of _unur_bessel_k_nuasympt() */
  
  // [[Rcpp::export]]
arma::rowvec rgig_cpp(int n, double lambda, double chi, double psi) {
  /*---------------------------------------------------------------------------*/
    /* Draw sample from GIG distribution.                                        */
    /* without calling GetRNGstate() ... PutRNGstate()                           */
    /*                                                                           */
    /* Parameters:                                                               */
    /*   n ....... sample size (positive integer)                                */
    /*   lambda .. parameter for distribution                                    */
    /*   chi   ... parameter for distribution                                    */
    /*   psi   ... parameter for distribution                                    */
    /*                                                                           */
    /* Return:                                                                   */
    /*   random sample of size 'n'                                               */
    /*---------------------------------------------------------------------------*/
    
    double omega, alpha;     /* parameters of standard distribution */
    arma::rowvec res(n);
  int i;
  
  /* check sample size */
    if (n<=0) {
      stop("sample size 'n' must be positive integer.");
    }
  
  /* check GIG parameters: */
    if ( !(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
         (chi <  0. || psi < 0)      || 
         (chi == 0. && lambda <= 0.) ||
         (psi == 0. && lambda >= 0.) ) {
      stop("invalid parameters for GIG distribution: lambda=%g, chi=%g, psi=%g",
           lambda, chi, psi);
    }
  
  if (chi < ZTOL) { 
    /* special cases which are basically Gamma and Inverse Gamma distribution */
      if (lambda > 0.0) {
        for (i=0; i<n; i++) res(i) = R::rgamma(lambda, 2.0/psi); 
      }
    else {
      for (i=0; i<n; i++) res(i) = 1.0/R::rgamma(-lambda, 2.0/psi); 
    }    
  }
  
  else if (psi < ZTOL) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
      if (lambda > 0.0) {
        for (i=0; i<n; i++) res(i) = 1.0/R::rgamma(lambda, 2.0/chi); 
      }
    else {
      for (i=0; i<n; i++) res(i) = R::rgamma(-lambda, 2.0/chi); 
    }    
    
  }
  
  else {
    double lambda_old = lambda;
    if (lambda < 0.) lambda = -lambda;
    alpha = sqrt(chi/psi);
    omega = sqrt(psi*chi);
    
    /* run generator */
      do {
        if (lambda > 2. || omega > 3.) {
          /* Ratio-of-uniforms with shift by 'mode', alternative implementation */
            _rgig_ROU_shift_alt(res, n, lambda, lambda_old, omega, alpha);
          break;
        }
        
        if (lambda >= 1.-2.25*omega*omega || omega > 0.2) {
          /* Ratio-of-uniforms without shift */
            _rgig_ROU_noshift(res, n, lambda, lambda_old, omega, alpha);
          break;
        }
        
        if (lambda >= 0. && omega > 0.) {
          /* New approach, constant hat in log-concave part. */
            _rgig_newapproach1(res, n, lambda, lambda_old, omega, alpha);
          break;
        }
        
        /* else */
          stop("parameters must satisfy lambda>=0 and omega>0.");
        
      } while (0);
  }
  
  /* return result */
    return res;
  
}




/****************************************************************************************
  * **************************** Multivariate Normal Sampling ****************************
  ****************************************************************************************/
  
  // [[Rcpp::export]]
arma::mat rmvn_cpp(uint32_t n,  
                   arma::colvec &mu_t,  
                   arma::mat &sigma,
                   int ncores,
                   bool isChol)
{ 
  // uint32_t n = as<uint32_t>(n_);
  // arma::rowvec mu = as<arma::rowvec>(mu_);  
  // arma::mat sigma = as<arma::mat>(sigma_); 
  // unsigned int  ncores = as<unsigned int>(ncores_); 
  // bool isChol = as<bool>(isChol_); 
  
  arma::rowvec mu = mu_t.t();
  unsigned int d = mu.n_elem;
  
  arma::mat A(n,d);
  
  try {
    if(n < 1) stop("n should be a positive integer");
    if(ncores == 0) stop("ncores has to be positive");
    if(d != sigma.n_cols) stop("mu.n_elem != sigma.n_cols");
    if(d != sigma.n_rows) stop("mu.n_elem != sigma.n_rows");
    if(d != A.n_cols) stop("mu.n_elem != A.ncol()");
    if(n != A.n_rows) stop("n != A.nrow()");
    
    // The A matrix that will be filled with firstly with standard normal rvs,
    // and finally with multivariate normal rvs.
    // We A wrap into a arma::mat "tmp" without making a copy.
    arma::mat tmp( A.begin(), A.n_rows, A.n_cols, false );
    
    RNGScope scope; // Declare RNGScope after the output in order to avoid a known Rcpp bug.
    
    // Calculate cholesky dec unless sigma is already a cholesky dec.
    arma::mat cholDec;
    if( isChol == false ) {
      cholDec = trimatu( arma::chol(sigma) );
    }
    else{
      cholDec = trimatu( sigma );
    }
    
    // What I do to seed the sitmo::prng_engine is tricky. I produce "ncores" uniform 
    // numbers between 1 and the largest uint32_t, which are the seeds. I put the first 
    // one in "coreSeed". If there is no support for OpenMP only this seed will be used, 
    // as the computations will be sequential. If there is support for OpenMP, "coreSeed" 
    // will be over-written, so that each core will get its own seed. 
    // NumericVector seeds = runif(ncores, 1.0, std::numeric_limits<uint32_t>::max());
    
    NumericVector seeds(ncores);
    seeds[0] = runif(1, 1.0, std::numeric_limits<uint32_t>::max())[0];
    for (unsigned int j = 0;  j < ncores - 1; j++){
      seeds[j+1] = seeds[j] - 1.0;
      if (seeds[j+1] < 1.0) seeds[j] = std::numeric_limits<uint32_t>::max() - 1.0;
    }
    
    #ifdef _OPENMP
    #pragma omp parallel num_threads(ncores) if(ncores > 1)
    {
      #endif
      
      double acc;
      uint32_t irow, icol, ii;
      arma::rowvec work(d);
      
      uint32_t coreSeed = static_cast<uint32_t>(seeds[0]);
      
      // (Optionally) over-writing the seed here
      #ifdef _OPENMP
      coreSeed = static_cast<uint32_t>( seeds[omp_get_thread_num()] );
      #endif
      
      sitmo::prng_engine engine( coreSeed );
      boost::normal_distribution<> normal(0.0, 1.0);
      
      // Filling "out" with standard normal rvs
      #ifdef _OPENMP
      #pragma omp for schedule(static)
      #endif
      for (irow = 0; irow < n; irow++) 
        for (icol = 0; icol < d; icol++) 
          A(irow, icol) = normal(engine);
      
      // Multiplying "out"" by cholesky decomposition of covariance and adding the
  // mean to obtain the desired multivariate normal data rvs.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for(irow = 0; irow < n; irow++) {
    for(icol = 0; icol < d; icol++) {
      acc = 0.0;
      for(ii = 0; ii <= icol; ii++) acc += tmp.at(irow, ii) * cholDec.at(ii, icol); 
      work.at(icol) = acc; 
    }
    tmp(arma::span(irow), arma::span::all) = work + mu;       
  }
  
#ifdef _OPENMP
}
#endif

//return R_NilValue;

  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  
  return A;
}

static double const log2pi = std::log(2.0*M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvn_cpp(arma::mat const &x,  
                   arma::rowvec const &mean,  
                   arma::mat const &sigma, 
                   bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

// [[Rcpp::export]]
double rgengamma_cpp(double shape, double scale, double p) {
  double sample = R::rgamma(shape, scale);
  return (std::pow(sample, p));
}

// [[Rcpp::export]]
double dgengamma_cpp(double x, double shape, double scale, double p, bool logd) {
  return (R::dgamma(std::pow(x,1./p), shape, scale, logd));
}

// [[Rcpp::export]]
arma::colvec dgengamma_vec_cpp(arma::vec &x, double shape, double scale, 
                               double p, bool logd) {
  NumericVector density_vec;
  density_vec = dgamma(pow(NumericVector(x.begin(), x.end()),1./p),
                       shape, scale, logd);
  arma::colvec density_vec_arma(density_vec.begin(), density_vec.size(), false);
  return(density_vec_arma);
}

double sum(arma::vec &v) {
  return arma::sum(v);
}

// [[Rcpp::export]]
arma::mat h_mat(arma::colvec tt, double rho) {
  int ni = tt.n_elem;
  arma::mat H_mat = arma::zeros(ni, ni);
  for (int i=0; i<ni; i++) {
    for (int j=0; j<ni; j++) {
      H_mat(i,j) = std::pow(rho, std::abs(tt(j)-tt(i)));
    }
  }
  return H_mat;
}

// [[Rcpp::export]]
arma::uvec D2B(int num, int size) {
  std::string err_msg = "error: need more than " + std::to_string(size) + 
    " bits to represent " + std::to_string(num) + "\n";
  if (num >= std::pow(2, size)) stop(err_msg);
  arma::colvec bin_pattern = arma::zeros(size);
  int temp = num;
  int pos = size-1;
  while (temp >= 2) {
    bin_pattern(pos--) = temp % 2;
    temp /= 2;
  }
  bin_pattern(pos) = temp;
  return (arma::conv_to<arma::uvec>::from(bin_pattern));
}



/********************************************************
 ********************************************************/






// [[Rcpp::export]]
double compute_pdf(
    double quad, arma::colvec alpha_arr, arma::uword k_z, double eta_z, 
    double sigma_02, arma::mat &H_inv, arma::uword n_i
) {
  double pdf;
  double zeta_k = 1./std::pow(alpha_arr(k_z), 2) - n_i/2.;
  // Rcout << "zeta_k: " << zeta_k << "\t, 2nd: " <<
  //   std::sqrt(2./sigma_02*quad)/(alpha_arr(k_z)*std::pow(eta_z, k_z)) << "\n";
  pdf =
    (-zeta_k-n_i) * std::log(
        alpha_arr(k_z)*std::pow(eta_z,k_z)*std::sqrt(sigma_02)
    ) +
      std::log(std::sqrt(arma::det(H_inv))) +
      -(n_i/2.+zeta_k/2.-1.)*std::log(2) +
      -(n_i/2.) * std::log(arma::datum::pi) +
      -std::lgamma(1./std::pow(alpha_arr(k_z),2.)) +
      zeta_k/2.*std::log(quad);
  if (zeta_k < 0) {
    pdf += std::log(boost::math::cyl_bessel_k(
      zeta_k,
      std::sqrt(2./sigma_02*quad)/(alpha_arr(k_z)*std::pow(eta_z, k_z))
    ));
  } else {
    pdf += bessel_k_nuasympt(
      std::sqrt(2./sigma_02*quad)/(alpha_arr(k_z)*std::pow(eta_z, k_z)),
      zeta_k,
      true,
      false
    );
  }
  return (pdf);
}



// [[Rcpp::export]]
double compute_pdf_II(
    double quad, arma::colvec alpha_arr, arma::uword k_z, 
    arma::colvec sigma2, arma::mat H_inv, double n_i
) {
  double pdf;
  if (alpha_arr(0) == 0 && k_z == 0) {
    pdf = 
      -(n_i/2.) * std::log(2*arma::datum::pi*sigma2(0)) +
      0.5 * std::log(arma::det(H_inv)) +
      -0.5 * quad/sigma2(0);
  } else {
    double zeta_k = 1./alpha_arr(k_z) - n_i/2.;
    // Rcout << "zeta_k: " << zeta_k << "\t, 2nd: " <<
    //   std::sqrt(2./sigma_02*quad)/(alpha_arr(k_z)*std::pow(eta_z, k_z)) << "\n";
    pdf =
      -(zeta_k+n_i) * std::log(
          std::sqrt(alpha_arr(k_z)*sigma2(k_z))
      ) +
        std::log(std::sqrt(arma::det(H_inv))) +
        -((n_i/2.)+zeta_k/2.-1.)*std::log(2) +
        -(n_i/2.) * std::log(arma::datum::pi) +
        -std::lgamma(1./alpha_arr(k_z)) +
        zeta_k/2.*std::log(quad);
    if (zeta_k < 50) {
      pdf += std::log(boost::math::cyl_bessel_k(
        zeta_k,
        std::sqrt(2./sigma2(k_z)*quad/alpha_arr(k_z))
      ));
    } else {
      pdf += bessel_k_nuasympt(
        std::sqrt(2./sigma2(k_z)*quad/alpha_arr(k_z)),
        zeta_k,
        true,
        false
      );
    }
  }
  return (pdf);
}


// [[Rcpp::export]]
SEXP EVIII_RE_sampler_II(
    List data, Nullable<Rcpp::List> true_params_, 
    arma::uvec fixed_eff_ind, arma::uvec rand_eff_ind,
    arma::uvec fixed_eff_u_ind, arma::uvec fixed_eff_W_ind, 
    arma::uvec fixed_eff_z_ind, 
    int n_sim, const int MAX_SAMPLE, int seed, 
    arma::mat Lambda,
    double eta_u, double eta_w, double eta_z, 
    bool fix_indicators, Nullable<Rcpp::List> indicators_,
    int n_iter1 = 500, int n_iter2 = 2500, int n_iter3 = 4000,
    double w_perc = 3., bool is_outlier_only = false
) {
  set_seed(seed);
  List data_list(data);
  // if (true_params_.isNull()) Rcout << "true_params is NULL";
  // if (indicators_.isNull()) Rcout << "indicators is NULL"; 
  
  List y_mat = as<Rcpp::List>(data_list["data"]);
  List tt = as<Rcpp::List>(data_list["time"]);
  List fixed_eff = as<Rcpp::List>(data_list["fixed.eff"]);
  
  arma::colvec est_u_arr, est_z_arr;
  // Rcout << "Fine1\n";
  List est_w_list = Rcpp::List::create();
  if (fix_indicators & indicators_.isNull())
    stop("Fixed indicators are not available. Aborted.\n");
  if (indicators_.isNotNull()) {
    List indc_list(indicators_);
    est_u_arr = as<arma::colvec>(indc_list["u"]);
    est_w_list = as<Rcpp::List>(indc_list["w"]);
    est_z_arr = as<arma::colvec>(indc_list["z"]);
  }
  // Rcout << "Fine2\n";
  
  // all the following variables should be declared outside and should be
  // defined within if(){}, if we want to access those later, to avoid 
  // compiler error
  Rcpp::List true_params_list(true_params_);
  arma::colvec true_beta, true_gamma_u, true_gamma_w, true_gamma_z;
  // arma::colvec true_gamma_u;
  // arma::colvec true_gamma_w;
  // arma::colvec true_gamma_z;
  Rcpp::List true_w_list;
  arma::colvec true_z_arr, true_u_arr;
  arma::mat true_Lambda, true_r_mat;
  double true_sigma02, true_rho;
  arma::colvec true_sigma_e2;
  
  if (true_params_.isNotNull()) {
    true_beta = as<arma::colvec>(true_params_list["beta"]);
    true_gamma_u = as<arma::colvec>(true_params_list["gamma_u"]);
    true_gamma_w = as<arma::colvec>(true_params_list["gamma_w"]);
    true_gamma_z = as<arma::colvec>(true_params_list["gamma_z"]);
    arma::mat true_Lambda = as<arma::mat>(true_params_list["Lambda"]);
    true_sigma02 = true_params_list["sigma_02"];
    true_rho = true_params_list["rho"];
    true_z_arr = as<arma::colvec>(data_list["z"]);
    true_u_arr = as<arma::colvec>(data_list["u"]);
    true_w_list = as<Rcpp::List>(data_list["w"]);
    true_sigma_e2 = as<arma::colvec>(data_list["sigma_i2"]);
    // Rcpp::List true_r_list = as<Rcpp::List>(data_list["r"]);
    true_r_mat = as<arma::mat>(data_list["r_mat"]);
  }
  
  
  const int n = y_mat.size();
  arma::mat J;
  
  /* ---- Initializations  ---- */
  arma::colvec ni =  arma::colvec(n);
  fixed_eff_ind = fixed_eff_ind - 1;
  int d_beta = fixed_eff_ind.n_elem;
  arma::colvec beta = arma::zeros(d_beta);
  arma::uvec irow(1), jrow(1);
  arma::uvec row_zero = {0};
  for (arma::uword i = 0; i < n; i++) {
    ni(i) = as<arma::vec>(tt[i]).n_elem;
  }
  
  // initializing W randomly
  int ni_max = ni.max();
  // arma::mat W = arma::mat(n,ni_max).fill(NA_REAL);
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < ni(i); j++) {
  //     W(i,j) = R::rbinom(1, 0.5);
  //   }
  // }
  arma::mat W = arma::mat(n, ni_max).fill(0);
  fixed_eff_W_ind = fixed_eff_W_ind - 1;
  int d_gam_W = fixed_eff_W_ind.n_elem;
  // arma::colvec gam_W = R::runif(-2, 2.) * arma::ones(d_gam_W);
  arma::colvec gam_W = arma::join_cols(arma::ones(1), arma::zeros(d_gam_W-1));
  arma::mat pi_W = arma::mat(n, ni_max).fill(NA_REAL);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < ni(i); j++) {
      pi_W(i,j) = R::runif(0,1);
    }
  }
  
  // arma::colvec z = RcppArmadillo::sample(arma::linspace(0,1,2), n, true);
  arma::colvec z = arma::colvec(n).fill(0);
  fixed_eff_z_ind = fixed_eff_z_ind - 1;
  int d_gam_h = fixed_eff_z_ind.n_elem;
  // arma::colvec gam_h = R::runif(-2.,2.)*arma::ones(d_gam_h);
  arma::colvec gam_h = arma::join_cols(arma::ones(1), arma::zeros(d_gam_h-1));
  arma::colvec pi_h = R::runif(0,1)*arma::ones(n);
  arma::colvec sigma_e2 = R::runif(0.1,5.)*arma::ones(n);
  bool flag_alpha; // used as an indicator whether proposed alpha is within limits
  // if (alpha == -1) {  // alpha not provided
  //   alpha = 1.;
  // }
  arma::colvec alpha_arr = {0., 1.};
  arma::colvec alpha_arr_new = arma::colvec(2);
  double mh_sig_alpha = 0.3;
  // double sigma_02 = R::runif(0.1,5.);
  arma::colvec sigma2_arr = {1., std::pow(eta_z,2)*1.};
  arma::colvec sigma2_arr_new = arma::colvec(2);
  arma::colvec mh_sig_sigma2 = {0.2, 0.5};
  // double sigma_12 = 0.4;
  // double eta_w = 4.;
  // double sigma_22 = std::pow(eta_w,2.)*sigma_02;
  // double mh_sig_sigma_22 = 0.4;
  // double sigma_12 = R::runif(0.1,5.);
  // double rho = R::runif(0,1);
  double rho = 0.;
  double mh_sig_rho = 0.08;
  // double sigma2 = 1.;
  
  // additonal things for the random effects
  rand_eff_ind = rand_eff_ind - 1;
  int d_rand = rand_eff_ind.n_elem;
  arma::mat r_mat = arma::zeros(n, d_rand);
  arma::colvec mu_r_i = arma::zeros(d_rand);
  arma::mat Sigma_r_i = arma::diagmat(arma::ones(d_rand));
  // arma::mat Lambda = arma::diagmat(arma::ones(d_rand));
  arma::mat Omega = arma::diagmat(arma::ones(d_rand));
  arma::mat Psi = arma::diagmat(arma::ones(d_rand));
  int nu = 4;
  
  
  arma::colvec ones_ni;
  // double eta_u = 1.;
  // arma::colvec u = RcppArmadillo::sample(arma::linspace(0,1,2), n, true);
  arma::colvec u = arma::colvec(n).fill(0);
  fixed_eff_u_ind = fixed_eff_u_ind - 1;
  int d_gam_u = fixed_eff_u_ind.n_elem;
  // arma::colvec gam_u = R::runif(-2.,2.)*arma::ones(d_gam_u);
  arma::colvec gam_u = arma::join_cols(arma::ones(1), arma::zeros(d_gam_u-1));
  arma::colvec pi_u = R::runif(0,1)*arma::ones(n);
  
  // double Lambda11 = Lambda(0,0);
  // double Lambda22 = Lambda(1,1);
  // double Lambda33 = Lambda(2,2);
  // double Lambda12 = Lambda(0,1) / std::pow(Lambda11*Lambda22, 0.5);
  // double Lambda13 = Lambda(0,2) / std::pow(Lambda11*Lambda33, 0.5);
  // double Lambda23 = Lambda(1,2) / std::pow(Lambda22*Lambda33, 0.5);
  // double Lambda11_new, Lambda22_new, Lambda33_new, Lambda12_new, Lambda13_new, Lambda23_new;
  // arma::mat Lambda_new(d_rand, d_rand);
  bool flag_Lambda;
  double mh_sig_Lambda_cor = 0.2;
  double mh_sig_Lambda_var = 0.2;
  
  /* ---- Storing Parameters ---- */
  arma::mat save_beta(n_sim, d_beta);
  arma::mat save_gam_W(n_sim, d_gam_W);
  arma::mat save_gam_h(n_sim, d_gam_h);
  arma::mat save_gam_u(n_sim, d_gam_u);
  arma::mat save_W(n_sim, n*ni_max);
  arma::mat save_pi_W(n_sim, n*ni_max);
  arma::mat save_z(n_sim, n);
  arma::mat save_pi_h(n_sim, n);
  arma::mat save_u(n_sim, n);
  arma::mat save_pi_u(n_sim, n);
  arma::mat save_sigma_e2(n_sim, n);
  arma::mat save_alpha(n_sim, 2);
  arma::mat save_sigma2(n_sim, 2);
  arma::colvec save_sigma_12(n_sim);
  arma::colvec save_sigma_22(n_sim);
  arma::colvec save_rho(n_sim);
  arma::mat save_r_mat(n_sim, n*d_rand);
  arma::mat save_Lambda(n_sim, d_rand*d_rand);
  // arma::colvec save_sigma2(n_sim);
  int acc_beta = 0;
  int acc_rho = 0;
  arma::colvec acc_alpha = {0., 0.};
  int acc_sigma_02 = 0;
  int acc_sigma_22 = 0;
  
  /* ---- Temporary variables---- */
  double zeta_k, zeta_i, zeta_i_new, xi, xi_new;
  arma::colvec beta_new(d_beta);
  arma::colvec mu_beta(d_beta);
  arma::mat inv_Sigma_beta(d_beta,d_beta), Sigma_beta(d_beta,d_beta);
  arma::mat fixed_eff_i, rand_eff_i;
  
  arma::rowvec fixed_eff_i_row;
  arma::colvec mu_gam_W(d_gam_W);
  arma::mat inv_Sigma_gam_W(d_gam_W,d_gam_W), Sigma_gam_W(d_gam_W,d_gam_W);
  arma::mat W_i;
  arma::colvec log_p_vec_0, log_p_vec_1, log_p_vec, p_vec;
  arma::colvec mu_gam_h(d_gam_h);
  arma::mat inv_Sigma_gam_h(d_gam_h,d_gam_h), Sigma_gam_h(d_gam_h,d_gam_h);
  arma::colvec log_p = arma::zeros(8);
  arma::colvec probs = arma::zeros(8);
  arma::colvec W_vec;
  int m_i = 0;
  arma::uvec m(n);
  arma::mat I_W_i;
  arma::uvec pos_bin, pos_vec, pos_vec_0, pos_vec_1;
  arma::uword z_i;
  
  double omega_i, omega_ij;
  arma::colvec zero_ni;
  double log_denom1, log_denom2, quad, quad_new, p1;
  arma::mat H_inv, H_inv_new;
  arma::colvec p;
  double sigma_12_new;
  arma::colvec ones_n = arma::ones(n);
  double mh_sig_sigma_12 = 0.25;
  double alpha_new;
  double sigma_02_new;
  double sigma2_new;
  double mh_log_lr = 0.;
  double rate = 0.;
  double rho_new;
  bool flag_rho;
  double flag_sampling;
  bool flag_sigma;
  
  arma::colvec mu_gam_u(d_gam_u);
  arma::mat inv_Sigma_gam_u(d_gam_u,d_gam_u), Sigma_gam_u(d_gam_u,d_gam_u);
  
  double p_mid = 0.05;
  // double p_low = 0.05;
  double logit_p_mid = std::log(p_mid/(1-p_mid));
  // double logit_p_low = std::log(p_low/(1-p_low));
  // double gam_intercept_sd = 
  //   (logit_p_mid - logit_p_low) / R::qnorm(0.975, 0., 1., true, false);
  double gam_intercept_sd = 1./std::pow(n/10*2*p_mid*(1-p_mid), 0.5);
  // Rcout << "gam_intercept_sd: " << gam_intercept_sd << "\n";
  double p_W_mid = w_perc/100;
  // double p_W_low = 0.01;
  double logit_p_W_mid = std::log(p_W_mid/(1-p_W_mid));
  // double logit_p_W_low = std::log(p_W_low/(1-p_W_low));
  // double gam_W_intercept_sd = 
  //   (logit_p_W_mid - logit_p_W_low) / R::qnorm(0.975, 0., 1., true, false);
  // Rcout << "gam_intercept_sd: " << gam_intercept_sd << "\n";
  double gam_W_intercept_sd = 1./std::pow(arma::accu(ni)/10*2*p_W_mid*(1-p_W_mid), 0.5);
  // Rcout << "gam_W_intercept_sd: " << gam_W_intercept_sd << "\n";
  arma::colvec prior_probs_arr = arma::colvec(n*ni_max*8).fill(NA_REAL);
  arma::colvec post_probs_arr = arma::colvec(n*ni_max*8).fill(NA_REAL);
  arma::colvec llhd_probs_arr = arma::colvec(n*ni_max*8).fill(NA_REAL);
  
  arma::mat save_prior_probs_mat = arma::mat(n_sim, n*ni_max*8).fill(NA_REAL);
  arma::mat save_post_probs_mat = arma::mat(n_sim, n*ni_max*8).fill(NA_REAL);
  arma::mat save_llhd_probs_mat = arma::mat(n_sim, n*ni_max*8).fill(NA_REAL);
  
  
  // Rcout << "logit_p_mid: " << logit_p_mid << "sd: " << gam_intercept_sd << "\n";
  // Rcout << "logit_p_W_mid: " << logit_p_W_mid << "sd: " << gam_W_intercept_sd << "\n";
  
  // int n_iter1 = 500;
  // int n_iter2 = 2500;
  // int n_iter3 = 4000;
  
  // Rcout << "Fine3\n";
  
  std::string err_msg;
  Progress pb(n_sim, true);
  
  // Rcout << "Fine3\n";
  int count_rho, count_alpha, count_sigma;
  for (arma::uword iter=0; iter<n_sim; iter++) {
    if (Progress::check_abort()) {
      err_msg = "Aborted.\n";
      stop(err_msg);
    }
    // Rcout << "iter: " << iter << std::endl;
    flag_rho = true;
    flag_alpha = true;
    flag_sigma = true;
    count_rho = count_alpha = count_sigma = 0;
    while(
      flag_rho || flag_alpha || flag_sigma
    // || flag_Lambda
    ) {
      // Rcout << 
      //   "sigma_02: " << sigma2_arr[0] << 
      //     "\t sigma_12: " << sigma2_arr[1] <<
      //       "\t alpha: " << alpha_arr[1] << "\n";
      if (count_rho > 1000) {
        err_msg = "Failed to sample rho. Aborted.\n";
        stop(err_msg);
      }
      if (count_alpha > 1000) {
        err_msg = "Failed to sample alpha. Aborted.\n";
        stop(err_msg);
      }
      if (count_sigma > 1000) {
        Rcout << 
          "sigma_02: " << sigma2_arr[0] << 
            "\t sigma_12: " << sigma2_arr[1] <<
              "\t alpha: " << alpha_arr[1] << "\n";
        err_msg = "Failed to sample sigma. Aborted.\n";
        stop(err_msg);
      }
      /* part of the update for rho, in case rho_new is outside [0,1] we skip
       the rest of the code to save time, otherwise make flag_rho = false so that this
       step will not be repeated and update rho wwth rho_new later 
       */
      rho_new = R::runif(rho-mh_sig_rho, rho+mh_sig_rho);
      if (rho_new < 0 || rho_new > 1) { count_rho++; continue; }
      flag_rho = false;
      
      /* same thing for alpha_arr as well */
      alpha_arr_new(0) = 0;
      // alpha_arr_new(0) = R::rlnorm(std::log(alpha_arr(0)), mh_sig_alpha);
      alpha_arr_new(1) = R::rlnorm(std::log(alpha_arr(1)), mh_sig_alpha);
      if (alpha_arr_new(1) > n*0.01) { count_alpha++; continue; }
      flag_alpha = false;
      
      /* also restrict sigma2_arr */
      sigma2_arr_new(0) = R::rlnorm(std::log(sigma2_arr(0)), mh_sig_sigma2(0));
      sigma2_arr_new(1) = R::rlnorm(std::log(sigma2_arr(1)), mh_sig_sigma2(1));
      if (sigma2_arr_new(0) > sigma2_arr_new(1)) { count_sigma++; continue; }
      flag_sigma = false;
      
      /* -- Updating gamma_W -- */
      if (iter > n_iter2) {
        // mu_gam = arma::zeros(d_gam); // gamma is noninformative
        // Rcout << "d_gam_W: " << d_gam_W << "\n";
        // inv_Sigma_gam_W = arma::diagmat(1./4.*arma::ones(d_gam_W));
        inv_Sigma_gam_W = arma::diagmat(arma::join_cols(
          1./std::pow(gam_W_intercept_sd,2.)*arma::ones(1), 1./0.01*arma::ones(d_gam_W-1)
        ));
        // mu_gam_W = inv_Sigma_gam_W * arma::zeros(d_gam_W);
        mu_gam_W = inv_Sigma_gam_W * arma::join_cols(
          logit_p_W_mid*arma::ones(1), arma::zeros(d_gam_W-1)
        );
        // Rcout << "d_gam_W: " << d_gam_W << "\n";
        for (arma::uword i=0; i<n; i++) {
          fixed_eff_i = as<arma::mat>(fixed_eff[i]);
          // Rcout << "fixed eff: " << fixed_eff_i_row << "\n";
          for (arma::uword j = 0; j < ni(i); j++) {
            jrow = {j};
            fixed_eff_i_row = fixed_eff_i.submat(jrow, fixed_eff_W_ind);
            omega_ij = 
              rcpp_pgdraw({1.}, {(fixed_eff_i_row*gam_W).eval()(0)})[0];
            inv_Sigma_gam_W += omega_ij * (fixed_eff_i_row.t()*fixed_eff_i_row);
            mu_gam_W += (W(i,j)-0.5) * fixed_eff_i_row.t();
          }
        }
        // Rcout << "d_gam_W: " << d_gam_W << "\n";
        Sigma_gam_W = arma::inv(inv_Sigma_gam_W);
        mu_gam_W = Sigma_gam_W * mu_gam_W;
        // Rcout << "Sigma_gam: " << Sigma_gam << "\n mu_gam: " << mu_gam << "\n";
        gam_W = rmvn_cpp(1, mu_gam_W, Sigma_gam_W, 1, false).row(0).t();
      }
      if (iter < 1200) {
        gam_W = arma::join_cols(
          logit_p_W_mid*arma::ones(1), arma::zeros(d_gam_W-1)
        );
      }
      // Rcout << "gamma_W: " << gam_W.t() << "\n";
      
      
      /* -- Updating gamma_h -- */
      if (iter > n_iter3) {
        // mu_gam = arma::zeros(d_gam); // gamma is noninformative
        // inv_Sigma_gam_h = arma::diagmat(1./4.*arma::ones(d_gam_h));
        inv_Sigma_gam_h = arma::diagmat(arma::join_cols(
          1./std::pow(gam_intercept_sd,2.)*arma::ones(1), 1./0.01*arma::ones(d_gam_h-1)
        ));
        // mu_gam_h = inv_Sigma_gam_h * arma::zeros(d_gam_h);
        mu_gam_h = inv_Sigma_gam_h * arma::join_cols(
          logit_p_mid*arma::ones(1), arma::zeros(d_gam_h-1)
        );
        // Rcout << "d_gam_h: " << d_gam_h << "\n";
        for (arma::uword i=0; i<n; i++) {
          fixed_eff_i = as<arma::mat>(fixed_eff[i]);
          fixed_eff_i_row = fixed_eff_i.submat(row_zero, fixed_eff_z_ind);
          // Rcout << "fixed eff: " << fixed_eff_i_row << "\n";
          omega_i = 
            rcpp_pgdraw({1.}, {(fixed_eff_i_row*gam_h).eval()(0)})[0];
          inv_Sigma_gam_h += omega_i * (fixed_eff_i_row.t()*fixed_eff_i_row);
          mu_gam_h += (z(i)-0.5)*fixed_eff_i_row.t();
        }
        Sigma_gam_h = arma::inv(inv_Sigma_gam_h);
        mu_gam_h = Sigma_gam_h * mu_gam_h;
        // Rcout << "Sigma_gam: " << Sigma_gam << "\n mu_gam: " << mu_gam << "\n";
        gam_h = rmvn_cpp(1, mu_gam_h, Sigma_gam_h, 1, false).row(0).t();
      }
      
      if (iter < 1200) {
        gam_h = arma::join_cols(
          logit_p_mid*arma::ones(1), arma::zeros(d_gam_h-1)
        );
      }
      // Rcout << "gamma_h: " << gam_h.t() << "\n";
      
      /* -- Updating gamma_u -- */
      if (iter > n_iter1) {
        // mu_gam = arma::zeros(d_gam); // gamma is noninformative
        // inv_Sigma_gam_u = arma::diagmat(1./100.*arma::ones(d_gam_u));
        inv_Sigma_gam_u = arma::diagmat(arma::join_cols(
          1./std::pow(gam_intercept_sd,2.)*arma::ones(1), 1./0.01*arma::ones(d_gam_u-1)
        ));
        // mu_gam_u = inv_Sigma_gam_u * arma::zeros(d_gam_u);
        mu_gam_u = inv_Sigma_gam_u * arma::join_cols(
          logit_p_mid*arma::ones(1), arma::zeros(d_gam_u-1)
        );
        // Rcout << "d_gam_u: " << d_gam_u << "\n";
        for (arma::uword i=0; i<n; i++) {
          fixed_eff_i = as<arma::mat>(fixed_eff[i]);
          fixed_eff_i_row = fixed_eff_i.submat(row_zero, fixed_eff_u_ind);
          // Rcout << "fixed eff: " << fixed_eff_i_row << "\n";
          omega_i = 
            rcpp_pgdraw({1.}, {(fixed_eff_i_row*gam_u).eval()(0)})[0];
          inv_Sigma_gam_u += omega_i * (fixed_eff_i_row.t()*fixed_eff_i_row);
          mu_gam_u += (u(i)-0.5)*fixed_eff_i_row.t();
        }
        Sigma_gam_u = arma::inv(inv_Sigma_gam_u);
        mu_gam_u = Sigma_gam_u * mu_gam_u;
        // Rcout << "Sigma_gam: " << Sigma_gam << "\n mu_gam: " << mu_gam << "\n";
        gam_u = rmvn_cpp(1, mu_gam_u, Sigma_gam_u, 1, false).row(0).t();
      }
      // gam_u = {-2.5, 0.25, 0.10};
      if (iter < 1200) {
        gam_u = arma::join_cols(
          logit_p_mid*arma::ones(1), arma::zeros(d_gam_u-1)
        );
      }
      // Rcout << "gamma_u: " << gam_u.t() << "\n";
      
      
      /* -- Updating u,z,w conditional on sigma_e2 -- */
      if (fix_indicators) {
        for (arma::uword i = 0; i < n; i++) {
          u(i) = est_u_arr(i);
          z(i) = est_z_arr(i);
          W_vec = as<arma::colvec>(est_w_list[i]);
          for (arma::uword j = 0; j < ni(i); j++)  { W(i,j) = W_vec(j); }
        }
      }
      if (iter > n_iter1 && !fix_indicators) {
        for (arma::uword i = 0; i < n; i++) {
          J = arma::ones(ni(i), ni(i));
          irow = {i};
          fixed_eff_i = as<arma::mat>(fixed_eff[i]);
          rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
          fixed_eff_i_row = fixed_eff_i.submat(row_zero, fixed_eff_z_ind);
          pi_h(i) = std::exp((fixed_eff_i_row*gam_h).eval()(0)) /
            (1. + std::exp((fixed_eff_i_row*gam_h).eval()(0)));
          fixed_eff_i_row = fixed_eff_i.submat(row_zero, fixed_eff_u_ind);
          pi_u(i) = std::exp((fixed_eff_i_row*gam_u).eval()(0)) /
            (1. + std::exp((fixed_eff_i_row*gam_u).eval()(0)));
          // if (i == 140) Rcout << "rand_eff_i: " << rand_eff_i << "\n";
          // if (i == 140) Rcout << "W(i,): " << W.row(i) << "\n";
          // if (i == 140) Rcout << "y[[141]]: " << as<arma::rowvec>(y_mat[140]) << std::endl;
          // if (i == 140) Rcout << "X_141 * beta: " << (fixed_eff_i * beta).t() << std::endl;
          // if (i == 140) Rcout << "Lambda: " << Lambda << std::endl;
          // if (i == 140) Rcout << "rho: " << rho << std::endl;
          // if (i == 140) Rcout << "Z*Lambda*Z': " << (rand_eff_i*Lambda*rand_eff_i.t()) << std::endl;
          // if (i == 140) Rcout << "Omega(w_i, rho): " << (eta_w*W_i + I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i) << std::endl;
          for (arma::uword j = 0; j < ni(i); j++) {
            jrow = {j};
            fixed_eff_i_row = fixed_eff_i.submat(jrow, fixed_eff_W_ind);
            pi_W(i,j) = std::exp((fixed_eff_i_row*gam_W).eval()(0)) /
              (1. + std::exp((fixed_eff_i_row*gam_W).eval()(0)));
            // Rcout << "pi_W(i,j): " << pi_W(i,j) << "\n";
            // Rcout << "W(i,): " << W.row(i) << "\n";
            for (arma::uword k_W = 0; k_W < 2; k_W++) {
              for (arma::uword k_z = 0; k_z < 2; k_z++) {
                for (arma::uword k_u = 0; k_u < 2; k_u++) {
                  // Rcout << "log_p[index]: " << log_p(2*k_W + k_z) << "\n";
                  // Rcout << "pi_h(i): " << pi_h(i) << ", pi_W(i,j): " << pi_W(i,j) << "\n";
                  // Rcout << "Test: " << k_W*std::log(pi_W(i,j)) + (1-k_W)*std::log(1-pi_W(i,j)) +
                  //   k_z*std::log(pi_h(i)) + (1-k_z)*std::log(1-pi_h(i)) << "\n";
                  log_p(4*k_W + 2*k_z + k_u) =
                    k_W*std::log(pi_W(i,j)) + (1-k_W)*std::log(1-pi_W(i,j)) +
                    k_z*std::log(pi_h(i)) + (1-k_z)*std::log(1-pi_h(i)) +
                    k_u*std::log(pi_u(i)) + (1-k_u)*std::log(1-pi_u(i));
                  prior_probs_arr(8*ni_max*i + 8*j + 4*k_W+2*k_z+k_u) = 
                    log_p(4*k_W + 2*k_z + k_u);
                  
                  W_vec = W.row(i).subvec(0,ni(i)-1).t();
                  // Rcout << "W_vec: " << W_vec << "\n";
                  W_vec(j) = k_W;
                  m_i = arma::accu(W_vec);
                  W_i = arma::diagmat(W_vec);
                  I_W_i = arma::diagmat(arma::ones(ni(i))) - W_i;
                  H_inv = arma::inv(
                    std::pow(std::pow(eta_u, 2.), k_u) * 
                      (rand_eff_i*Lambda*rand_eff_i.t()) +
                      (eta_w*W_i + I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i)
                  );
                  // Rcout << "H_inv: " << H_inv << "\n";
                  // Rcout << "beta: " << beta << "\n";
                  quad = (
                    (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()) *
                      H_inv *
                      (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()).t()
                  ).eval()(0);
                  // Rcout << "alpha: " << alpha << ", sigma_12: " << sigma_12 <<
                  //   ", quad: " << quad << "\n";
                  zeta_k = 1./alpha_arr(k_z) - ni(i)/2.;
                  llhd_probs_arr(i*8*ni_max+8*j + 4*k_W+2*k_z+k_u) = compute_pdf_II (
                    quad, alpha_arr, k_z, sigma2_arr, H_inv, ni(i)
                  ); 
                  log_p(4*k_W + 2*k_z + k_u) += 
                    llhd_probs_arr(i*8*ni_max+8*j + 4*k_W+2*k_z+k_u);
                  post_probs_arr(8*ni_max*i + 8*j + 4*k_W+2*k_z+k_u) = 
                    log_p(4*k_W+2*k_z+k_u);
                  // if (i == 140 & j == 9) Rcout << "W(i,): " << W.row(i) << std::endl;
                  // if (i == 140 & j == 9) Rcout << "Z_i: " << W.row(i) << std::endl;
                  // if (i == 140 & j == 9) Rcout << "W(i,): " << W.row(i) << std::endl;
                  // if (i == 140 & j == 9) Rcout << "H: " << std::pow(std::pow(eta_u, 2.), k_u) * 
                  // (rand_eff_i*Lambda*rand_eff_i.t()) +
                  // (eta_w*W_i + I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i) << std::endl;
                  // if (i == 140 & j == 9) Rcout << "H_inv: " << H_inv << std::endl;
                  // if (i == 140 & j == 9) Rcout << "quad: " << quad << std::endl;
                }
              }
            }
            if (log_p.has_nan() || arma::accu(arma::exp(log_p)) == 0.) {
              Rcout << "i: " << i << "\tlog_p: " << log_p.t() << "\n";
              Rcout << "W_vec: " << W.row(i).subvec(0,ni(i)-1).t() << "\n";
              Rcout << "m_i: " << m_i << "\n";
              Rcout << "beta: " << beta.t() << "\n";
              Rcout << "gamma_W: " << gam_W.t() << "\n";
              Rcout << "gamma_h: " << gam_h.t() << "\n";
              Rcout << "alpha: " << alpha_arr.t() << "\n";
              Rcout << "sigma2_arr: " << sigma2_arr.t() << "\n";
              Rcout << "rho: " << rho << "\n";
              continue;
            }
            
            // log_p = log_p - arma::max(log_p);
            // probs = arma::exp(log_p);
            prior_probs_arr(i*8*ni_max+8*j + arma::linspace<arma::uvec>(0,7,8)) = 
              arma::exp(
                prior_probs_arr(
                  i*8*ni_max+8*j + arma::linspace<arma::uvec>(0,7,8)
                ) -
                  arma::max(prior_probs_arr(
                      i*8*ni_max+8*j + arma::linspace<arma::uvec>(0,7,8)
                  ))
              );
            
            
            if (iter <= n_iter2) {
              // fitting hom/het intercept and hom variance and no outlier model
              log_p(arma::linspace<arma::uvec>(2,7,6)).fill(-arma::datum::inf);
              // probs(2) = 0; probs(3) = 0; 
              // probs(4) = 0; probs(5) = 0; probs(6) = 0; probs(7) = 0;
              prior_probs_arr(i*8*ni_max+8*j + arma::linspace<arma::uvec>(2,7,6)).fill(0);
            } else if (iter <= n_iter3) {
              // fitting hom/het intercept and hom variance with outlier model
              log_p(2) = -arma::datum::inf;
              log_p(3) = -arma::datum::inf;
              log_p(6) = -arma::datum::inf;
              log_p(7) = -arma::datum::inf;
              // probs(2) = 0; probs(3) = 0;
              // probs(6) = 0; probs(7) = 0;
              prior_probs_arr(i*8*ni_max+8*j + 2) = 0;
              prior_probs_arr(i*8*ni_max+8*j + 3) = 0;
              prior_probs_arr(i*8*ni_max+8*j + 6) = 0;
              prior_probs_arr(i*8*ni_max+8*j + 7) = 0;
            }
            if (is_outlier_only) {
              log_p(arma::linspace<arma::uvec>(1,3,3)).fill(-arma::datum::inf);
              log_p(arma::linspace<arma::uvec>(5,7,3)).fill(-arma::datum::inf);
              prior_probs_arr(i*8*ni_max+8*j + 1) = 0;
              prior_probs_arr(i*8*ni_max+8*j + 2) = 0;
              prior_probs_arr(i*8*ni_max+8*j + 3) = 0;
              prior_probs_arr(i*8*ni_max+8*j + 5) = 0;
              prior_probs_arr(i*8*ni_max+8*j + 6) = 0;
              prior_probs_arr(i*8*ni_max+8*j + 7) = 0;
            }
            
            log_p = log_p - arma::max(log_p);
            probs = arma::exp(log_p);
            post_probs_arr(i*8*ni_max+8*j+ arma::linspace<arma::uvec>(0,7,8)) = probs;
            int z_W = RcppArmadillo::sample(arma::linspace(0,7,8), 1, false, probs)(0);
            
            if (z_W == 0)       { W(i,j) = 0; z(i) = 0; u(i) = 0; }
            else if (z_W == 1)  { W(i,j) = 0; z(i) = 0; u(i) = 1; }
            else if (z_W == 2)  { W(i,j) = 0; z(i) = 1; u(i) = 0; }
            else if (z_W == 3)  { W(i,j) = 0; z(i) = 1; u(i) = 1; }
            else if (z_W == 4)  { W(i,j) = 1; z(i) = 0; u(i) = 0; }
            else if (z_W == 5)  { W(i,j) = 1; z(i) = 0; u(i) = 1; }
            else if (z_W == 6)  { W(i,j) = 1; z(i) = 1; u(i) = 0; }
            else if (z_W == 7)  { W(i,j) = 1; z(i) = 1; u(i) = 1; }
          }
          m(i) = arma::accu(arma::conv_to<arma::uvec>::from(W.row(i)).subvec(0, ni(i)-1));
        }
      }
      // Rcout << "W(i,): " << W.row(0) << "\n";
      // Rcout << "m(i): " << m(0) << "\n";
      
      
      /* -- Updating sigma_02-- */
      for (arma::uword ia=0; ia<2; ia++) {
        // sigma2_arr_new = sigma2_arr;
        // sigma2_arr_new(ia) = R::rlnorm(std::log(sigma2_arr(ia)), mh_sig_sigma2(ia));
        // proposal
        mh_log_lr =
          R::dlnorm(sigma2_arr(ia), std::log(sigma2_arr_new(ia)), mh_sig_sigma2(ia), true) -
          R::dlnorm(sigma2_arr_new(ia), std::log(sigma2_arr(ia)), mh_sig_sigma2(ia), true);
        if (ia == 0) {
          // prior
          mh_log_lr +=
            (R::dgamma(1./sigma2_arr_new(0), 0.1, 10., true) - 2*std::log(sigma2_arr_new(0))) -
            (R::dgamma(1./sigma2_arr(0), 0.1, 10., true) - 2*std::log(sigma2_arr(0)));
          // contribution from sigma_1^2
          if (iter > n_iter3) {
            xi = 0.01 * n / alpha_arr(1);
            mh_log_lr += (
              R::dgamma(1./sigma2_arr(1), xi, 1./((xi-1.)*std::pow(eta_z,2)*sigma2_arr_new(0)), true) -
                2*std::log(sigma2_arr(1))
            ) - (
                R::dgamma(1./sigma2_arr(1), xi, 1./((xi-1.)*std::pow(eta_z,2)*sigma2_arr(0)), true) -
                  2*std::log(sigma2_arr(1))
            );
          }
        } else if (ia == 1) {
          // prior
          xi = 0.01 * n / alpha_arr(ia);
          mh_log_lr += (
            R::dgamma(1./sigma2_arr_new(ia), xi, 1./((xi-1.)*std::pow(eta_z,2.)*sigma2_arr(0)), true) -
              2*std::log(sigma2_arr_new(ia))
          ) - (
              R::dgamma(1./sigma2_arr(ia), xi, 1./((xi-1.)*std::pow(eta_z,2.)*sigma2_arr(0)), true) -
                2*std::log(sigma2_arr(ia))
          );
        }
        for (arma::uword i = 0; i < n; i++) {
          // only considering likelihoods from relevant i's such that z(i) = ia
          if (z(i) == ia) {
            fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
            rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
            W_i = arma::diagmat(W.row(i).subvec(0,ni(i)-1));
            I_W_i = arma::diagmat(arma::ones(ni(i))) - W_i;
            J = arma::ones(ni(i), ni(i));
            H_inv = arma::inv(
              std::pow(eta_u,(2.*u(i))) * (rand_eff_i*Lambda*rand_eff_i.t()) +
                (eta_w*W_i+I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i)
            );
            // Rcout << "eta_u: " << eta_u << "\t u(i):" << u(i) << "\n";
            // Rcout << "Lambda: " << Lambda << "\n";
            // Rcout << "W_i: " << W_i << "\n";
            // Rcout << "I_W_i: " << I_W_i << "\n";
            // Rcout << "rho: " << rho << "\n";
            // Rcout << "H_inv: " << H_inv << "\n";
            // Rcout << "beta: " << beta << "\n";
            quad = (
              (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()) *
                H_inv *
                (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()).t()
            ).eval()(0);
            zeta_i = 1./alpha_arr(z(i)) - ni(i)/2.;
            // Rcout << "zeta_i: " << zeta_i << "\tquad: " << quad << "\n";
            // likelihood
            mh_log_lr +=
              compute_pdf_II(quad, alpha_arr, z(i), sigma2_arr_new, H_inv, ni(i)) -
              compute_pdf_II(quad, alpha_arr, z(i), sigma2_arr, H_inv, ni(i));
          }
        }
        if (R::runif(0,1) < std::min(1., std::exp(mh_log_lr))) {
          acc_sigma_02++;
          sigma2_arr(ia) = sigma2_arr_new(ia);
        }
        // sigma2_arr(0) = 0.4;
      }
      // sigma2_arr(0) = 0.4;
      // sigma2_arr(1) = 3.6;
      // Rcout << "sigma2_arr: " << sigma2_arr.t() << "\n";
      
      
      /* -- Updating rho -- */
      // we already have rho_new at the very beginning
      mh_log_lr = 0;
      for (arma::uword i=0; i<n; i++) {
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        W_i = arma::diagmat(W.row(i).subvec(0,ni(i)-1));
        I_W_i = arma::diagmat(arma::ones(ni(i))) - W_i;
        J = arma::ones(ni(i), ni(i));
        
        H_inv = arma::inv(
          std::pow(eta_u, 2.*u(i)) * (rand_eff_i*Lambda*rand_eff_i.t()) +
            (eta_w*W_i + I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i)
        );
        H_inv_new = arma::inv(
          std::pow(eta_u, 2.*u(i)) * (rand_eff_i*Lambda*rand_eff_i.t()) +
            (eta_w*W_i + I_W_i) * h_mat(tt[i],rho_new) * (eta_w*W_i + I_W_i)
        );
        quad =
          ((as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()) *
          H_inv *
          (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()).t()).eval()(0);
        quad_new =
          ((as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()) *
          H_inv_new *
          (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()).t()).eval()(0);
        zeta_i = 1./alpha_arr(z(i)) - ni(i)/2.;
        mh_log_lr +=
          compute_pdf_II(quad_new, alpha_arr, z(i), sigma2_arr, H_inv_new, ni(i)) -
          compute_pdf_II(quad, alpha_arr, z(i), sigma2_arr, H_inv, ni(i));
      }
      if (R::runif(0,1) < std::min(1., std::exp(mh_log_lr))) {
        acc_rho++;
        rho = rho_new;
      }
      // rho = 0.2;
      // if (iter <= 4000) rho = 0.1;
      // Rcout << "rho: " << rho << "\n";
      
      /* -- Updating alpha -- */
      for (arma::uword ia=0; ia<2; ia++) {
        if (ia == 0) { alpha_arr(ia) = 0; continue; }
        alpha_new = alpha_arr_new(ia);
        mh_log_lr =
          (R::dlnorm(alpha_arr(ia), std::log(alpha_arr_new(ia)), mh_sig_alpha, true) +
          R::dexp(alpha_arr_new(ia), 1., true)) -
          (R::dlnorm(alpha_arr_new(ia), std::log(alpha_arr(ia)), mh_sig_alpha, true) +
          R::dexp(alpha_arr(ia), 1., true));
        if (ia == 1) {
          xi_new = 0.01 * n / alpha_arr_new(ia);
          xi = 0.01 * n / alpha_arr(ia);
          mh_log_lr += (
            R::dgamma(
              1./sigma2_arr(ia), xi_new, 1./((xi_new-1)*std::pow(eta_z,2)*sigma2_arr(0)), true
            ) - 2*std::log(sigma2_arr(ia))
          ) - (
              R::dgamma(
                1./sigma2_arr(ia), xi, 1./((xi-1)*std::pow(eta_z,2)*sigma2_arr(0)), true
              ) - 2*std::log(sigma2_arr(ia))
          );
        }
        for (arma::uword i=0; i<n; i++) {
          if (z(i) == ia) {
            fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
            rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
            W_i = arma::diagmat(W.row(i).subvec(0,ni(i)-1));
            I_W_i = arma::diagmat(arma::ones(ni(i))) - W_i;
            J = arma::ones(ni(i), ni(i));
            H_inv = arma::inv(
              std::pow(eta_u, 2.*u(i)) * (rand_eff_i*Lambda*rand_eff_i.t()) +
                (eta_w*W_i + I_W_i)*
                h_mat(tt[i],rho)*
                (eta_w*W_i + I_W_i)
            );
            quad =
              ((as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()) *
              H_inv *
              (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()).t()).eval()(0);
            zeta_i = 1./alpha_arr(ia) - ni(i)/2.;
            zeta_i_new = 1./alpha_arr_new(ia) - ni(i)/2.;
            mh_log_lr += 
              compute_pdf_II(quad, alpha_arr_new, ia, sigma2_arr, H_inv, ni(i)) - 
              compute_pdf_II(quad, alpha_arr, ia, sigma2_arr, H_inv, ni(i));
            
          }
        }
        if (R::runif(0,1) < std::min(1., std::exp(mh_log_lr))) {
          acc_alpha(ia)++;
          alpha_arr(ia) = alpha_new;
        }
      }
      
      
      /* -- Updating Lambda -- */
      Omega = Psi;
      for (arma::uword i=0; i<n; i++) {
        Omega += r_mat.row(i).t() * r_mat.row(i) /
          (sigma_e2(i) * std::pow(eta_u, 2.*u(i)));
      }
      Lambda = riwish(nu+n, Omega);
      
      
      /* -- Updating r_mat -- */
      for (arma::uword i=0; i<n; i++) {
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        W_i = arma::diagmat(W.row(i).subvec(0,ni(i)-1));
        I_W_i = arma::diagmat(arma::ones(ni(i))) - W_i;
        Sigma_r_i = sigma_e2(i) * arma::inv(
          std::pow(1./eta_u,2.*u(i)) * arma::inv(Lambda) +
            rand_eff_i.t() * arma::inv(
                (eta_w*W_i + I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i)
            ) * rand_eff_i
        );
        mu_r_i = Sigma_r_i / sigma_e2(i) * (rand_eff_i.t()) * arma::inv(
          (eta_w*W_i + I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i)
        ) * (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta);
        r_mat.row(i) = rmvn_cpp(1, mu_r_i, Sigma_r_i, 1, false).row(0);
        // r_mat.row(i) = as<arma::rowvec>(true_r_list[i]);
        // r_mat.row(i) = true_r_mat.row(i);
      }
      
      /* -- Updating sigma_e2 -- */
      for (arma::uword i=0; i<n; i++) {
        if (z(i) == 0) { sigma_e2(i) = sigma2_arr(0); continue; }
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        W_i = arma::diagmat(W.row(i).subvec(0,ni(i)-1));
        I_W_i = arma::diagmat(arma::ones(ni(i))) - W_i;
        J = arma::ones(ni(i), ni(i));
        // z_i = arma::conv_to<arma::uword>::from(z(i));
        z_i = z(i);
        sigma_e2(i) = rgig_cpp(
          1,
          1./alpha_arr(z_i) - ni(i)/2.,
          ((as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()) * 
            arma::inv(
              std::pow(eta_u, 2.*u(i)) * (rand_eff_i*Lambda*rand_eff_i.t()) +
                (eta_w*W_i + I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i)
            ) * (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta).t()).t()).eval()(0),
            2./(alpha_arr(z_i)*sigma2_arr(z_i))
        )(0);
      }
      // Rcout << "sum(sigma_e2): " << arma::accu(sigma_e2) << "\n";
      
      
      /* -- Updating beta -- */
      mu_beta = arma::zeros(d_beta);
      inv_Sigma_beta = arma::diagmat(arma::ones(d_beta));
      for (arma::uword i=0; i<n; i++) {
        irow = {i};
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        W_i = arma::diagmat(W.row(i).subvec(0,ni(i)-1));
        I_W_i = arma::diagmat(arma::ones(ni(i))) - W_i;
        J = arma::ones(ni(i), ni(i));
        inv_Sigma_beta += fixed_eff_i.t() * arma::inv(sigma_e2(i) * (
          std::pow(eta_u, 2.*u(i)) * (rand_eff_i*Lambda*rand_eff_i.t()) +
            (eta_w*W_i + I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i)
        )) * fixed_eff_i;
        // mu_beta +=
        //   fixed_eff_i.t() * arma::inv(sigma_e2(i)*h_mat(tt[i], rho)) *
        //   y_mat.submat(irow, as<arma::uvec>(tt[i])-1).t();
        mu_beta += fixed_eff_i.t() * arma::inv(sigma_e2(i) * (
          std::pow(eta_u, 2.*u(i)) * (rand_eff_i*Lambda*rand_eff_i.t()) +
            (eta_w*W_i + I_W_i) * h_mat(tt[i],rho) * (eta_w*W_i + I_W_i)
        )) * (as<arma::colvec>(y_mat[i]));
      }
      Sigma_beta = arma::inv(inv_Sigma_beta);
      mu_beta = Sigma_beta * mu_beta;
      beta = rmvn_cpp(1, mu_beta, Sigma_beta, 1, false).row(0).t();
      // beta = {5., 2., -1., 0.8, -0.04};
      // beta = true_beta;
      // Rcout << "beta: " << beta << "\n";
      
    }
    
    /* ---- Saving Parameters ---- */
    save_beta.row(iter) = beta.t();
    save_gam_h.row(iter) = gam_h.t();
    save_gam_u.row(iter) = gam_u.t();
    save_gam_W.row(iter) = gam_W.t();
    save_sigma_e2.row(iter) = sigma_e2.t();
    save_W.row(iter) = W.as_row();
    save_z.row(iter) = z.t();
    save_u.row(iter) = u.t();
    save_pi_W.row(iter) = pi_W.as_row();
    save_pi_h.row(iter) = pi_h.t();
    save_pi_u.row(iter) = pi_u.t();
    save_prior_probs_mat.row(iter) = prior_probs_arr.t();
    save_post_probs_mat.row(iter) = post_probs_arr.t();
    save_llhd_probs_mat.row(iter) = llhd_probs_arr.t();
    save_alpha.row(iter) = alpha_arr.t();
    save_sigma2.row(iter) = sigma2_arr.t();
    save_rho(iter) = rho;
    save_r_mat.row(iter) = r_mat.as_col().t();
    save_Lambda.row(iter) = Lambda.as_col().t();
    
    // Rcout << "Fine: " << iter << "\n";
    pb.increment();
  }
  
  return (
      List::create(
        Named("acc") = 
          List::create(Named("beta") = (double)acc_beta/n_sim*100,
                       Named("alpha") = acc_alpha/n_sim*100,
                       Named("rho") = (double)acc_rho/n_sim*100,
                       Named("sigma_02") = (double)acc_sigma_02/n_sim*100),
                       Named("mcmc.df") = 
                         List::create(Named("beta") = save_beta,
                                      Named("gamma.z") = save_gam_h,
                                      Named("gamma.u") = save_gam_u,
                                      Named("gamma.w") = save_gam_W,
                                      Named("sigma.e2") = save_sigma_e2,
                                      Named("w") = save_W,
                                      Named("z") = save_z,
                                      Named("u") = save_u,
                                      Named("alpha") = save_alpha,
                                      Named("sigma2") = save_sigma2,
                                      Named("rho") = save_rho,
                                      Named("r_mat") = save_r_mat,
                                      Named("Lambda") = save_Lambda)));
}



// [[Rcpp::export]]
SEXP EVIII_RE_sampler_t1 (
  List data, Nullable<Rcpp::List> true_params_, 
  arma::uvec fixed_eff_ind, arma::uvec rand_eff_ind,
  int n_sim, const int MAX_SAMPLE, int seed, 
  arma::mat Lambda
) {
  // set_seed(seed);
  List data_list(data);
  // if (true_params_.isNull()) Rcout << "true_params is NULL";
  // if (indicators_.isNull()) Rcout << "indicators is NULL"; 
  
  List y_mat = as<Rcpp::List>(data_list["data"]);
  List tt = as<Rcpp::List>(data_list["time"]);
  List fixed_eff = as<Rcpp::List>(data_list["fixed.eff"]);
  
  // all the following variables should be declared outside and should be
  // defined within if(){}, if we want to access those later, to avoid 
  // compiler error
  Rcpp::List true_params_list(true_params_);
  arma::colvec true_beta;
  arma::mat true_Lambda, true_r_mat;
  double true_sigma02, true_rho;
  arma::colvec true_sigma_e2;
  
  if (true_params_.isNotNull()) {
    true_beta = as<arma::colvec>(true_params_list["beta"]);
    arma::mat true_Lambda = as<arma::mat>(true_params_list["Lambda"]);
    true_sigma02 = true_params_list["sigma_02"];
    true_rho = true_params_list["rho"];
    true_sigma_e2 = as<arma::colvec>(data_list["sigma_i2"]);
    // Rcpp::List true_r_list = as<Rcpp::List>(data_list["r"]);
    true_r_mat = as<arma::mat>(data_list["r_mat"]);
  }
  
  
  const int n = y_mat.size();
  arma::mat J;
  
  /* ---- Initializations  ---- */
    arma::colvec ni =  arma::colvec(n);
  fixed_eff_ind = fixed_eff_ind - 1;
  int d_beta = fixed_eff_ind.n_elem;
  arma::colvec beta = arma::zeros(d_beta);
  arma::uvec irow(1), jrow(1);
  arma::uvec row_zero = {0};
  for (arma::uword i = 0; i < n; i++) {
    ni(i) = as<arma::vec>(tt[i]).n_elem;
  }
  
  // initializing W randomly
  int ni_max = ni.max();
  arma::colvec sigma_e2 = R::runif(0.1,5.)*arma::ones(n);
  double sigma_02 = R::runif(0.1,5.);
  arma::colvec tau2 = arma::ones(n);
  double mh_sig_sigma_02 = 0.2;
  // double sigma_12 = 0.4;
  // double eta_w = 4.;
  // double sigma_22 = std::pow(eta_w,2.)*sigma_02;
  // double mh_sig_sigma_22 = 0.4;
  // double sigma_12 = R::runif(0.1,5.);
  // double rho = R::runif(0,1);
  double rho = 0.;
  double mh_sig_rho = 0.08;
  // double sigma2 = 1.;
  
  // additonal things for the random effects
  rand_eff_ind = rand_eff_ind - 1;
  int d_rand = rand_eff_ind.n_elem;
  arma::mat r_mat = arma::zeros(n, d_rand);
  arma::colvec mu_r_i = arma::zeros(d_rand);
  arma::mat Sigma_r_i = arma::diagmat(arma::ones(d_rand));
  arma::mat Omega = arma::diagmat(arma::ones(d_rand));
  arma::mat Psi = arma::diagmat(arma::ones(d_rand));
  
  // arma::colvec ones_ni;
  
  /* ---- Storing Parameters ---- */
    arma::mat save_beta(n_sim, d_beta);
  arma::colvec save_sigma_02(n_sim);
  arma::mat save_tau2(n_sim, n);
  arma::colvec save_rho(n_sim);
  arma::mat save_r_mat(n_sim, n*d_rand);
  arma::mat save_Lambda(n_sim, d_rand*d_rand);
  // arma::colvec save_sigma2(n_sim);
  int acc_beta = 0;
  int acc_rho = 0;
  int acc_sigma_02 = 0;
  
  /* ---- Temporary variables---- */
    double zeta_k, zeta_i, zeta_i_new, xi, xi_new;
  arma::colvec beta_new(d_beta);
  arma::colvec mu_beta(d_beta);
  arma::mat inv_Sigma_beta(d_beta,d_beta), Sigma_beta(d_beta,d_beta);
  arma::mat fixed_eff_i, rand_eff_i;
  
  arma::rowvec fixed_eff_i_row;
  arma::colvec zero_ni;
  double mh_log_lr = 0.;
  double rate = 0.;
  double rho_new;
  bool flag_rho;
  double flag_sampling;
  
  double h1 = 0.01;
  double h2 = 0.01;
  double h1_new, h2_new;
  double nu = 5;
  double t1_new, t2_new;
  
  
  // Rcout << "logit_p_mid: " << logit_p_mid << "sd: " << gam_intercept_sd << "\n";
  // Rcout << "logit_p_W_mid: " << logit_p_W_mid << "sd: " << gam_W_intercept_sd << "\n";
  
  // int n_iter1 = 500;
  // int n_iter2 = 2500;
  // int n_iter3 = 4000;
  
  // Rcout << "Fine3\n";
  
  std::string err_msg;
  Progress pb(n_sim, true);
  
  // Rcout << "Fine3\n";
  int count_rho, count_alpha, count_sigma;
  for (arma::uword iter=0; iter<n_sim; iter++) {
    if (Progress::check_abort()) {
      err_msg = "Aborted.\n";
      stop(err_msg);
    }
    // Rcout << "iter: " << iter << std::endl;
    flag_rho = true;
    count_rho = 0;
    while(flag_rho) {
      if (count_rho > 1000) {
        err_msg = "Failed to sample rho. Aborted.\n";
        stop(err_msg);
      }
      /* part of the update for rho, in case rho_new is outside [0,1] we skip
      the rest of the code to save time, otherwise make flag_rho = false so that this
      step will not be repeated and update rho wwth rho_new later 
      */
        rho_new = R::runif(rho-mh_sig_rho, rho+mh_sig_rho);
      if (rho_new < 0 || rho_new > 1) { count_rho++; continue; }
      flag_rho = false;
      
      /* -- Updating sigma_02-- */
        h1_new = h1;
      h2_new = h2;
      for (arma::uword i=0; i<n; i++) {
        h1_new += ni(i)/2.;
        irow = {i};
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        h2_new += 0.5 * (
          (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta).t() * 
            inv(
              tau2(i)*rand_eff_i*Lambda*rand_eff_i.t() + h_mat(tt[i],rho)
            ) * 
            (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta)
        ).eval()(0);
      }
      sigma_02 = 1./R::rgamma(h1_new, 1./h2_new);
      
      
      /* -- Updating rho -- */
        // we already have rho_new at the very beginning
      mh_log_lr = 0;
      for (arma::uword i=0; i<n; i++) {
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        
        mh_log_lr += dmvn_cpp(
          as<arma::rowvec>(y_mat[i]),
          (fixed_eff_i * beta).t(),
          sigma_02 * (tau2(i)*rand_eff_i*Lambda*rand_eff_i.t() + h_mat(tt[i],rho)),
          true
        )(0) - dmvn_cpp(
          as<arma::rowvec>(y_mat[i]),
          (fixed_eff_i * beta).t(),
          sigma_02 * (tau2(i)*rand_eff_i*Lambda*rand_eff_i.t() + h_mat(tt[i],rho_new)),
          true
        )(0);
      }
      if (R::runif(0,1) < std::min(1., std::exp(mh_log_lr))) {
        acc_rho++;
        rho = rho_new;
      }
      
      
      /* -- Updating Lambda -- */
        Omega = Psi;
      for (arma::uword i=0; i<n; i++) {
        Omega += r_mat.row(i).t() * r_mat.row(i) / (sigma_02 * tau2(i));
      }
      Lambda = riwish(nu+n, Omega);
      // Lambda.row(0) = {1, 0.05};
      // Lambda.row(1) = {0.05, 0.1};
      // Rcout << "Lambda: " << Lambda << "\n";
      
      
      /* -- Updating r_mat -- */
        for (arma::uword i=0; i<n; i++) {
          fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
          rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
          Sigma_r_i = sigma_02 * arma::inv(
            1./tau2(i) * arma::inv(Lambda) +
              rand_eff_i.t() * arma::inv(h_mat(tt[i],rho)) * rand_eff_i
          );
          mu_r_i = Sigma_r_i / (sigma_02) * 
            rand_eff_i.t() * arma::inv(h_mat(tt[i],rho)) * 
            (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta);
          r_mat.row(i) = rmvn_cpp(1, mu_r_i, Sigma_r_i, 1, false).row(0);
          // r_mat.row(i) = as<arma::rowvec>(true_r_list[i]);
          // r_mat.row(i) = true_r_mat.row(i);
        }
      // Rcout << "r_mat[1,]: " << r_mat.row(0) << "\n";
      
      /* -- Updating tau2(i) -- */
        for (arma::uword i=0; i<n; i++) {
          t1_new = nu/2. + ni(i)/2.;
          irow = {i};
          fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
          rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
          t2_new = nu/2. + 0.5 * (
            r_mat.row(i) * inv(sigma_02*Lambda) * r_mat.row(i).t()
          ).eval()(0);
          tau2(i) = 1./R::rgamma(t1_new, 1./t2_new);
        }
      
      
      /* -- Updating beta -- */
        mu_beta = arma::zeros(d_beta);
      inv_Sigma_beta = arma::diagmat(arma::ones(d_beta));
      for (arma::uword i=0; i<n; i++) {
        irow = {i};
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        inv_Sigma_beta += fixed_eff_i.t() * arma::inv(sigma_02 * (
          tau2(i)*rand_eff_i*Lambda*rand_eff_i.t() + h_mat(tt[i],rho)
        )) * fixed_eff_i;
        mu_beta += fixed_eff_i.t() * arma::inv(sigma_02 * (
          tau2(i)*rand_eff_i*Lambda*rand_eff_i.t() + h_mat(tt[i],rho)
        )) * (as<arma::colvec>(y_mat[i]));
      }
      Sigma_beta = arma::inv(inv_Sigma_beta);
      mu_beta = Sigma_beta * mu_beta;
      beta = rmvn_cpp(1, mu_beta, Sigma_beta, 1, false).row(0).t();
      // beta = {5., 2., -1., 0.8, -0.04};
      // beta = true_beta;
      // Rcout << "beta: " << beta << "\n";
      
    }
    
    /* ---- Saving Parameters ---- */
      save_beta.row(iter) = beta.t();
      save_tau2.row(iter) = tau2.t();
      save_sigma_02(iter) = sigma_02;
      save_rho(iter) = rho;
      save_r_mat.row(iter) = r_mat.as_col().t();
      save_Lambda.row(iter) = Lambda.as_col().t();
      
      // Rcout << "Fine: " << iter << "\n";
      pb.increment();
  }
  
  return (
    List::create(
      Named("acc") = 
        List::create(Named("beta") = (double)acc_beta/n_sim*100,
                     Named("rho") = (double)acc_rho/n_sim*100,
                     Named("sigma_02") = (double)acc_sigma_02/n_sim*100),
      Named("mcmc.df") = 
        List::create(Named("beta") = save_beta,
                     Named("tau2") = save_tau2,
                     Named("sigma02") = save_sigma_02,
                     Named("rho") = save_rho,
                     Named("r_mat") = save_r_mat,
                     Named("Lambda") = save_Lambda)));
}



// [[Rcpp::export]]
SEXP EVIII_RE_sampler_t2 (
  List data, Nullable<Rcpp::List> true_params_, 
  arma::uvec fixed_eff_ind, arma::uvec rand_eff_ind,
  int n_sim, const int MAX_SAMPLE, int seed, 
  arma::mat Lambda
) {
  // set_seed(seed);
  List data_list(data);
  // if (true_params_.isNull()) Rcout << "true_params is NULL";
  // if (indicators_.isNull()) Rcout << "indicators is NULL"; 
  
  List y_mat = as<Rcpp::List>(data_list["data"]);
  List tt = as<Rcpp::List>(data_list["time"]);
  List fixed_eff = as<Rcpp::List>(data_list["fixed.eff"]);
  
  // all the following variables should be declared outside and should be
  // defined within if(){}, if we want to access those later, to avoid 
  // compiler error
  Rcpp::List true_params_list(true_params_);
  arma::colvec true_beta;
  arma::mat true_Lambda, true_r_mat;
  double true_sigma02, true_rho;
  arma::colvec true_sigma_e2;
  
  if (true_params_.isNotNull()) {
    true_beta = as<arma::colvec>(true_params_list["beta"]);
    arma::mat true_Lambda = as<arma::mat>(true_params_list["Lambda"]);
    true_sigma02 = true_params_list["sigma_02"];
    true_rho = true_params_list["rho"];
    true_sigma_e2 = as<arma::colvec>(data_list["sigma_i2"]);
    // Rcpp::List true_r_list = as<Rcpp::List>(data_list["r"]);
    true_r_mat = as<arma::mat>(data_list["r_mat"]);
  }
  
  
  const int n = y_mat.size();
  arma::mat J;
  
  /* ---- Initializations  ---- */
    arma::colvec ni =  arma::colvec(n);
  fixed_eff_ind = fixed_eff_ind - 1;
  int d_beta = fixed_eff_ind.n_elem;
  arma::colvec beta = arma::zeros(d_beta);
  arma::uvec irow(1), jrow(1);
  arma::uvec row_zero = {0};
  for (arma::uword i = 0; i < n; i++) {
    ni(i) = as<arma::vec>(tt[i]).n_elem;
  }
  
  // initializing W randomly
  int ni_max = ni.max();
  arma::colvec sigma_e2 = R::runif(0.1,5.)*arma::ones(n);
  double sigma_02 = R::runif(0.1,5.);
  arma::colvec tau2 = arma::ones(n);
  double mh_sig_sigma_02 = 0.2;
  // double sigma_12 = 0.4;
  // double eta_w = 4.;
  // double sigma_22 = std::pow(eta_w,2.)*sigma_02;
  // double mh_sig_sigma_22 = 0.4;
  // double sigma_12 = R::runif(0.1,5.);
  // double rho = R::runif(0,1);
  double rho = 0.;
  double mh_sig_rho = 0.08;
  // double sigma2 = 1.;
  
  // additonal things for the random effects
  rand_eff_ind = rand_eff_ind - 1;
  int d_rand = rand_eff_ind.n_elem;
  arma::mat r_mat = arma::zeros(n, d_rand);
  arma::colvec mu_r_i = arma::zeros(d_rand);
  arma::mat Sigma_r_i = arma::diagmat(arma::ones(d_rand));
  arma::mat Omega = arma::diagmat(arma::ones(d_rand));
  arma::mat Psi = arma::diagmat(arma::ones(d_rand));
  
  // arma::colvec ones_ni;
  
  /* ---- Storing Parameters ---- */
    arma::mat save_beta(n_sim, d_beta);
  arma::colvec save_sigma_02(n_sim);
  arma::mat save_tau2(n_sim, n);
  arma::colvec save_rho(n_sim);
  arma::mat save_r_mat(n_sim, n*d_rand);
  arma::mat save_Lambda(n_sim, d_rand*d_rand);
  // arma::colvec save_sigma2(n_sim);
  int acc_beta = 0;
  int acc_rho = 0;
  int acc_sigma_02 = 0;
  
  /* ---- Temporary variables---- */
    double zeta_k, zeta_i, zeta_i_new, xi, xi_new;
  arma::colvec beta_new(d_beta);
  arma::colvec mu_beta(d_beta);
  arma::mat inv_Sigma_beta(d_beta,d_beta), Sigma_beta(d_beta,d_beta);
  arma::mat fixed_eff_i, rand_eff_i;
  
  arma::rowvec fixed_eff_i_row;
  arma::colvec zero_ni;
  double mh_log_lr = 0.;
  double rate = 0.;
  double rho_new;
  bool flag_rho;
  double flag_sampling;
  
  double h1 = 0.01;
  double h2 = 0.01;
  double h1_new, h2_new;
  double nu = 5;
  double t1_new, t2_new;
  
  
  // Rcout << "logit_p_mid: " << logit_p_mid << "sd: " << gam_intercept_sd << "\n";
  // Rcout << "logit_p_W_mid: " << logit_p_W_mid << "sd: " << gam_W_intercept_sd << "\n";
  
  // int n_iter1 = 500;
  // int n_iter2 = 2500;
  // int n_iter3 = 4000;
  
  // Rcout << "Fine3\n";
  
  std::string err_msg;
  Progress pb(n_sim, true);
  
  // Rcout << "Fine3\n";
  int count_rho, count_alpha, count_sigma;
  for (arma::uword iter=0; iter<n_sim; iter++) {
    if (Progress::check_abort()) {
      err_msg = "Aborted.\n";
      stop(err_msg);
    }
    // Rcout << "iter: " << iter << std::endl;
    flag_rho = true;
    count_rho = 0;
    while(flag_rho) {
      if (count_rho > 1000) {
        err_msg = "Failed to sample rho. Aborted.\n";
        stop(err_msg);
      }
      /* part of the update for rho, in case rho_new is outside [0,1] we skip
      the rest of the code to save time, otherwise make flag_rho = false so that this
      step will not be repeated and update rho wwth rho_new later 
      */
        rho_new = R::runif(rho-mh_sig_rho, rho+mh_sig_rho);
      if (rho_new < 0 || rho_new > 1) { count_rho++; continue; }
      flag_rho = false;
      
      /* -- Updating sigma_02-- */
        h1_new = h1/2.;
      h2_new = h2/2.;
      for (arma::uword i=0; i<n; i++) {
        h1_new += ni(i)/2.;
        irow = {i};
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        h2_new += 0.5 * (
          (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta).t() * 
            inv(
              rand_eff_i*Lambda*rand_eff_i.t() + tau2(i)*h_mat(tt[i],rho)
            ) * 
            (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta)
        ).eval()(0);
      }
      sigma_02 = 1./R::rgamma(h1_new, 1./h2_new);
      
      
      /* -- Updating rho -- */
        // we already have rho_new at the very beginning
      mh_log_lr = 0;
      for (arma::uword i=0; i<n; i++) {
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        
        mh_log_lr += dmvn_cpp(
          as<arma::rowvec>(y_mat[i]),
          (fixed_eff_i * beta).t(),
          sigma_02* (rand_eff_i*Lambda*rand_eff_i.t() + tau2(i)*h_mat(tt[i],rho)),
          true
        )(0) - dmvn_cpp(
          as<arma::rowvec>(y_mat[i]),
          (fixed_eff_i * beta).t(),
          sigma_02* (rand_eff_i*Lambda*rand_eff_i.t() + tau2(i)*h_mat(tt[i],rho_new)),
          true
        )(0);
      }
      if (R::runif(0,1) < std::min(1., std::exp(mh_log_lr))) {
        acc_rho++;
        rho = rho_new;
      }
      
      
      /* -- Updating Lambda -- */
        Omega = Psi;
      for (arma::uword i=0; i<n; i++) {
        Omega += r_mat.row(i).t() * r_mat.row(i) / sigma_02;
      }
      Lambda = riwish(nu+n, Omega);
      // Lambda.row(0) = {1, 0.05};
      // Lambda.row(1) = {0.05, 0.1};
      // Rcout << "Lambda: " << Lambda << "\n";
      
      
      /* -- Updating r_mat -- */
        for (arma::uword i=0; i<n; i++) {
          fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
          rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
          Sigma_r_i = sigma_02 * arma::inv(
            arma::inv(Lambda) +
              rand_eff_i.t() * arma::inv(tau2(i)*h_mat(tt[i],rho)) * rand_eff_i
          );
          mu_r_i = Sigma_r_i / sigma_02 * 
            rand_eff_i.t() * arma::inv(tau2(i)*h_mat(tt[i],rho)) * 
            (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta);
          r_mat.row(i) = rmvn_cpp(1, mu_r_i, Sigma_r_i, 1, false).row(0);
          // r_mat.row(i) = as<arma::rowvec>(true_r_list[i]);
          // r_mat.row(i) = true_r_mat.row(i);
        }
      // Rcout << "r_mat[1,]: " << r_mat.row(0) << "\n";
      
      /* -- Updating tau2(i) -- */
        for (arma::uword i=0; i<n; i++) {
          t1_new = nu/2. + ni(i)/2.;
          irow = {i};
          fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
          rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
          t2_new = nu/2. + 0.5 * (
            (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta - rand_eff_i*r_mat.row(i).t()).t() * 
              inv(sigma_02 * h_mat(tt[i],rho)) * 
              (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta - rand_eff_i*r_mat.row(i).t())
          ).eval()(0);
          tau2(i) = 1./R::rgamma(t1_new, 1./t2_new);
        }
      
      
      /* -- Updating beta -- */
        mu_beta = arma::zeros(d_beta);
      inv_Sigma_beta = arma::diagmat(arma::ones(d_beta));
      for (arma::uword i=0; i<n; i++) {
        irow = {i};
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        inv_Sigma_beta += fixed_eff_i.t() * arma::inv(sigma_02 * (
          rand_eff_i*Lambda*rand_eff_i.t() + tau2(i)*h_mat(tt[i],rho)
        )) * fixed_eff_i;
        mu_beta += fixed_eff_i.t() * arma::inv(sigma_02 * (
          rand_eff_i*Lambda*rand_eff_i.t() + tau2(i)*h_mat(tt[i],rho)
        )) * (as<arma::colvec>(y_mat[i]));
      }
      Sigma_beta = arma::inv(inv_Sigma_beta);
      mu_beta = Sigma_beta * mu_beta;
      beta = rmvn_cpp(1, mu_beta, Sigma_beta, 1, false).row(0).t();
      // beta = {5., 2., -1., 0.8, -0.04};
      // beta = true_beta;
      // Rcout << "beta: " << beta << "\n";
      
    }
    
    /* ---- Saving Parameters ---- */
      save_beta.row(iter) = beta.t();
      save_tau2.row(iter) = tau2.t();
      save_sigma_02(iter) = sigma_02;
      save_rho(iter) = rho;
      save_r_mat.row(iter) = r_mat.as_col().t();
      save_Lambda.row(iter) = Lambda.as_col().t();
      
      // Rcout << "Fine: " << iter << "\n";
      pb.increment();
  }
  
  return (
    List::create(
      Named("acc") = 
        List::create(Named("beta") = (double)acc_beta/n_sim*100,
                     Named("rho") = (double)acc_rho/n_sim*100,
                     Named("sigma_02") = (double)acc_sigma_02/n_sim*100),
      Named("mcmc.df") = 
        List::create(Named("beta") = save_beta,
                     Named("tau2") = save_tau2,
                     Named("sigma02") = save_sigma_02,
                     Named("rho") = save_rho,
                     Named("r_mat") = save_r_mat,
                     Named("Lambda") = save_Lambda)));
}


// [[Rcpp::export]]
SEXP EVIII_RE_sampler_t3 (
  List data, Nullable<Rcpp::List> true_params_, 
  arma::uvec fixed_eff_ind, arma::uvec rand_eff_ind,
  int n_sim, const int MAX_SAMPLE, int seed, 
  arma::mat Lambda
) {
  // set_seed(seed);
  List data_list(data);
  // if (true_params_.isNull()) Rcout << "true_params is NULL";
  // if (indicators_.isNull()) Rcout << "indicators is NULL"; 
  
  List y_mat = as<Rcpp::List>(data_list["data"]);
  List tt = as<Rcpp::List>(data_list["time"]);
  List fixed_eff = as<Rcpp::List>(data_list["fixed.eff"]);
  
  // all the following variables should be declared outside and should be
  // defined within if(){}, if we want to access those later, to avoid 
  // compiler error
  Rcpp::List true_params_list(true_params_);
  arma::colvec true_beta;
  arma::mat true_Lambda, true_r_mat;
  double true_sigma02, true_rho;
  arma::colvec true_sigma_e2;
  
  if (true_params_.isNotNull()) {
    true_beta = as<arma::colvec>(true_params_list["beta"]);
    arma::mat true_Lambda = as<arma::mat>(true_params_list["Lambda"]);
    true_sigma02 = true_params_list["sigma_02"];
    true_rho = true_params_list["rho"];
    true_sigma_e2 = as<arma::colvec>(data_list["sigma_i2"]);
    // Rcpp::List true_r_list = as<Rcpp::List>(data_list["r"]);
    true_r_mat = as<arma::mat>(data_list["r_mat"]);
  }
  
  
  const int n = y_mat.size();
  arma::mat J;
  
  /* ---- Initializations  ---- */
    arma::colvec ni =  arma::colvec(n);
  fixed_eff_ind = fixed_eff_ind - 1;
  int d_beta = fixed_eff_ind.n_elem;
  arma::colvec beta = arma::zeros(d_beta);
  arma::uvec irow(1), jrow(1);
  arma::uvec row_zero = {0};
  for (arma::uword i = 0; i < n; i++) {
    ni(i) = as<arma::vec>(tt[i]).n_elem;
  }
  
  // initializing W randomly
  int ni_max = ni.max();
  arma::colvec sigma_e2 = R::runif(0.1,5.)*arma::ones(n);
  double sigma_02 = R::runif(0.1,5.);
  arma::colvec tau12 = arma::ones(n);
  arma::colvec tau22 = arma::ones(n);
  double mh_sig_sigma_02 = 0.2;
  // double sigma_12 = 0.4;
  // double eta_w = 4.;
  // double sigma_22 = std::pow(eta_w,2.)*sigma_02;
  // double mh_sig_sigma_22 = 0.4;
  // double sigma_12 = R::runif(0.1,5.);
  // double rho = R::runif(0,1);
  double rho = 0.;
  double mh_sig_rho = 0.08;
  // double sigma2 = 1.;
  
  // additonal things for the random effects
  rand_eff_ind = rand_eff_ind - 1;
  int d_rand = rand_eff_ind.n_elem;
  arma::mat r_mat = arma::zeros(n, d_rand);
  arma::colvec mu_r_i = arma::zeros(d_rand);
  arma::mat Sigma_r_i = arma::diagmat(arma::ones(d_rand));
  arma::mat Omega = arma::diagmat(arma::ones(d_rand));
  arma::mat Psi = arma::diagmat(arma::ones(d_rand));
  
  // arma::colvec ones_ni;
  
  /* ---- Storing Parameters ---- */
    arma::mat save_beta(n_sim, d_beta);
  arma::colvec save_sigma_02(n_sim);
  arma::mat save_tau12(n_sim, n);
  arma::mat save_tau22(n_sim, n);
  arma::colvec save_rho(n_sim);
  arma::mat save_r_mat(n_sim, n*d_rand);
  arma::mat save_Lambda(n_sim, d_rand*d_rand);
  // arma::colvec save_sigma2(n_sim);
  int acc_beta = 0;
  int acc_rho = 0;
  int acc_sigma_02 = 0;
  
  /* ---- Temporary variables---- */
    double zeta_k, zeta_i, zeta_i_new, xi, xi_new;
  arma::colvec beta_new(d_beta);
  arma::colvec mu_beta(d_beta);
  arma::mat inv_Sigma_beta(d_beta,d_beta), Sigma_beta(d_beta,d_beta);
  arma::mat fixed_eff_i, rand_eff_i;
  
  arma::rowvec fixed_eff_i_row;
  arma::colvec zero_ni;
  double mh_log_lr = 0.;
  double rate = 0.;
  double rho_new;
  bool flag_rho;
  double flag_sampling;
  
  double h1 = 0.01;
  double h2 = 0.01;
  double h1_new, h2_new;
  double nu = 5;
  double t1_new, t2_new;
  
  
  // Rcout << "logit_p_mid: " << logit_p_mid << "sd: " << gam_intercept_sd << "\n";
  // Rcout << "logit_p_W_mid: " << logit_p_W_mid << "sd: " << gam_W_intercept_sd << "\n";
  
  // int n_iter1 = 500;
  // int n_iter2 = 2500;
  // int n_iter3 = 4000;
  
  // Rcout << "Fine3\n";
  
  std::string err_msg;
  Progress pb(n_sim, true);
  
  // Rcout << "Fine3\n";
  int count_rho, count_alpha, count_sigma;
  for (arma::uword iter=0; iter<n_sim; iter++) {
    if (Progress::check_abort()) {
      err_msg = "Aborted.\n";
      stop(err_msg);
    }
    // Rcout << "iter: " << iter << std::endl;
    flag_rho = true;
    count_rho = 0;
    while(flag_rho) {
      if (count_rho > 1000) {
        err_msg = "Failed to sample rho. Aborted.\n";
        stop(err_msg);
      }
      /* part of the update for rho, in case rho_new is outside [0,1] we skip
      the rest of the code to save time, otherwise make flag_rho = false so that this
      step will not be repeated and update rho wwth rho_new later 
      */
        rho_new = R::runif(rho-mh_sig_rho, rho+mh_sig_rho);
      if (rho_new < 0 || rho_new > 1) { count_rho++; continue; }
      flag_rho = false;
      
      /* -- Updating sigma_02-- */
        h1_new = h1/2.;
      h2_new = h2/2.;
      for (arma::uword i=0; i<n; i++) {
        h1_new += ni(i)/2.;
        irow = {i};
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        h2_new += 0.5 * (
          (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta).t() * 
            inv(
              tau22(i) * rand_eff_i*Lambda*rand_eff_i.t() + tau12(i)*h_mat(tt[i],rho)
            ) * 
            (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta)
        ).eval()(0);
      }
      sigma_02 = 1./R::rgamma(h1_new, 1./h2_new);
      
      
      /* -- Updating rho -- */
        // we already have rho_new at the very beginning
      mh_log_lr = 0;
      for (arma::uword i=0; i<n; i++) {
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        
        mh_log_lr += dmvn_cpp(
          as<arma::rowvec>(y_mat[i]),
          (fixed_eff_i * beta).t(),
          sigma_02* (tau22(i)*rand_eff_i*Lambda*rand_eff_i.t() + tau12(i)*h_mat(tt[i],rho)),
          true
        )(0) - dmvn_cpp(
          as<arma::rowvec>(y_mat[i]),
          (fixed_eff_i * beta).t(),
          sigma_02* (
            tau22(i)*rand_eff_i*Lambda*rand_eff_i.t() + tau12(i)*h_mat(tt[i],rho_new)
          ),
          true
        )(0);
      }
      if (R::runif(0,1) < std::min(1., std::exp(mh_log_lr))) {
        acc_rho++;
        rho = rho_new;
      }
      
      
      /* -- Updating Lambda -- */
        Omega = Psi;
      for (arma::uword i=0; i<n; i++) {
        Omega += r_mat.row(i).t() * r_mat.row(i) / (sigma_02 * tau22(i));
      }
      Lambda = riwish(nu+n, Omega);
      // Lambda.row(0) = {1, 0.05};
      // Lambda.row(1) = {0.05, 0.1};
      // Rcout << "Lambda: " << Lambda << "\n";
      
      
      /* -- Updating r_mat -- */
        for (arma::uword i=0; i<n; i++) {
          fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
          rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
          Sigma_r_i = sigma_02 * arma::inv(
            1./tau22(i) * arma::inv(Lambda) +
              rand_eff_i.t() * 1./tau12(i) * arma::inv(h_mat(tt[i],rho)) * rand_eff_i
          );
          mu_r_i = Sigma_r_i / (sigma_02 * tau12(i)) * 
            rand_eff_i.t() * arma::inv(h_mat(tt[i],rho)) * 
            (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta);
          r_mat.row(i) = rmvn_cpp(1, mu_r_i, Sigma_r_i, 1, false).row(0);
          // r_mat.row(i) = as<arma::rowvec>(true_r_list[i]);
          // r_mat.row(i) = true_r_mat.row(i);
        }
      // Rcout << "r_mat[1,]: " << r_mat.row(0) << "\n";
      
      /* -- Updating tau22(i) -- */
        for (arma::uword i=0; i<n; i++) {
          t1_new = nu/2. + d_rand/2.;
          irow = {i};
          fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
          rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
          t2_new = nu/2. + 0.5 * (
            r_mat.row(i) * 
              inv(sigma_02 * Lambda) * 
              r_mat.row(i).t()
          ).eval()(0);
          tau22(i) = 1./R::rgamma(t1_new, 1./t2_new);
        }
      /* -- Updating tau12(i) -- */
        for (arma::uword i=0; i<n; i++) {
          t1_new = nu/2. + ni(i)/2.;
          irow = {i};
          fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
          rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
          t2_new = nu/2. + 0.5 * (
            (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta - rand_eff_i*r_mat.row(i).t()).t() * 
              inv(sigma_02 * h_mat(tt[i],rho)) * 
              (as<arma::colvec>(y_mat[i]) - fixed_eff_i*beta - rand_eff_i*r_mat.row(i).t())
          ).eval()(0);
          tau12(i) = 1./R::rgamma(t1_new, 1./t2_new);
        }
      
      
      /* -- Updating beta -- */
        mu_beta = arma::zeros(d_beta);
      inv_Sigma_beta = arma::diagmat(arma::ones(d_beta));
      for (arma::uword i=0; i<n; i++) {
        irow = {i};
        fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind);
        rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind);
        inv_Sigma_beta += fixed_eff_i.t() * arma::inv(sigma_02* (
          tau22(i) * rand_eff_i*Lambda*rand_eff_i.t() + tau12(i) * h_mat(tt[i],rho)
        )) * fixed_eff_i;
        mu_beta += fixed_eff_i.t() * arma::inv(sigma_02* (
          tau22(i) * rand_eff_i*Lambda*rand_eff_i.t() + tau12(i) * h_mat(tt[i],rho)
        )) * (as<arma::colvec>(y_mat[i]));
      }
      Sigma_beta = arma::inv(inv_Sigma_beta);
      mu_beta = Sigma_beta * mu_beta;
      beta = rmvn_cpp(1, mu_beta, Sigma_beta, 1, false).row(0).t();
      // beta = {5., 2., -1., 0.8, -0.04};
      // beta = true_beta;
      // Rcout << "beta: " << beta << "\n";
      
    }
    
    /* ---- Saving Parameters ---- */
      save_beta.row(iter) = beta.t();
      save_tau12.row(iter) = tau12.t();
      save_tau22.row(iter) = tau22.t();
      save_sigma_02(iter) = sigma_02;
      save_rho(iter) = rho;
      save_r_mat.row(iter) = r_mat.as_col().t();
      save_Lambda.row(iter) = Lambda.as_col().t();
      
      // Rcout << "Fine: " << iter << "\n";
      pb.increment();
  }
  
  return (
    List::create(
      Named("acc") = 
        List::create(Named("beta") = (double)acc_beta/n_sim*100,
                     Named("rho") = (double)acc_rho/n_sim*100,
                     Named("sigma_02") = (double)acc_sigma_02/n_sim*100),
      Named("mcmc.df") = 
        List::create(Named("beta") = save_beta,
                     Named("tau12") = save_tau12,
                     Named("tau22") = save_tau22,
                     Named("sigma02") = save_sigma_02,
                     Named("rho") = save_rho,
                     Named("r_mat") = save_r_mat,
                     Named("Lambda") = save_Lambda)));
}




//[[Rcpp::export]]
double sum_of_exp(arma::colvec log_arr) {
  double sum_arr, max_log_arr;
  max_log_arr = arma::max(log_arr);
  // Rcout << "min: " << min_log_arr;
  sum_arr = max_log_arr + std::log(arma::accu(arma::exp(log_arr - max_log_arr)));
  return (sum_arr);
}


// [[Rcpp::export]]
double compute_WAIC_EVIII_RE_marginal_1(
  List dat, List model, List W_ind_list, arma::uvec &ind, 
  arma::uvec &fixed_eff_ind, arma::uvec &rand_eff_ind, 
  arma::uvec fixed_eff_u_ind, arma::uvec &fixed_eff_W_ind,
  arma::uvec &fixed_eff_z_ind, 
  double eta_w, double eta_u, double eta_z, int mod_type, int is_conditional
) {
  set_seed(123456);
  List data_list(dat);
  List y_mat = as<Rcpp::List>(data_list["data"]);
  List fixed_eff = as<Rcpp::List>(data_list["fixed.eff"]);
  List tt = as<Rcpp::List>(data_list["time"]);
  ind = ind - 1;
  
  const int n = y_mat.size();
  
  int n_sim = ind.n_elem;
  arma::colvec ni(n);
  for (arma::uword i = 0; i < n; i++) {
    ni(i) = as<arma::vec>(tt[i]).n_elem;
  }
  int ni_max = ni.max();
  
  List model_list(model);
  arma::mat beta = as<arma::mat>(model_list["beta"]);
  beta = beta.rows(ind);
  // arma::mat r_mat = as<arma::mat>(model_list["r_mat"]);
  // r_mat = r_mat.rows(ind);
  arma::mat sigma2 = as<arma::mat>(model_list["sigma2"]);
  sigma2 = sigma2.rows(ind);
  arma::mat alpha = as<arma::mat>(model_list("alpha"));
  alpha = alpha.rows(ind);
  arma::colvec rho = as<arma::colvec>(model_list("rho"));
  rho = rho(ind);
  arma::mat gamma_h = as<arma::mat>(model_list("gamma.z"));
  gamma_h = gamma_h.rows(ind);
  arma::mat gamma_W = as<arma::mat>(model_list("gamma.w"));
  gamma_W = gamma_W.rows(ind);
  arma::mat gamma_u = as<arma::mat>(model_list("gamma.u"));
  gamma_u = gamma_u.rows(ind);
  arma::mat Lambda = as<arma::mat>(model_list("Lambda"));
  Lambda = Lambda.rows(ind);
  // arma::mat sigma_i2_mat = as<arma::mat>(model_list("sigma.e2"));
  // sigma_i2_mat = sigma_i2_mat.rows(ind);
  arma::mat u_mat = as<arma::mat>(model_list("u"));
  u_mat = u_mat.rows(ind);
  arma::mat W_mat = as<arma::mat>(model_list("w"));
  W_mat = W_mat.rows(ind);
  arma::mat z_mat = as<arma::mat>(model_list("z"));
  z_mat= z_mat.rows(ind);
  arma::colvec W_i_vec = arma::colvec(ni_max);
  double pi_h;
  arma::colvec pi_W;
  double pi_u;
  arma::mat J;
  List W_ind(W_ind_list);
  arma::colvec W_ind_i;
  double prior_i_out_marg = 0.;
  int d_rand = rand_eff_ind.n_elem;
  double zeta;
  arma::colvec sigma2_arr(2);
  
  arma::uvec irow(1);
  arma::mat fixed_eff_i;
  arma::mat rand_eff_i;
  arma::mat fixed_eff_z_i;
  arma::mat fixed_eff_u_i;
  arma::mat fixed_eff_W_i;
  arma::colvec fixed_eff_i_row;
  
  arma::mat Lambda_mat;
  arma::uvec pos_vec, pos_bin;
  arma::mat W_i, I_W_i;
  int m_i;
  arma::mat H_inv;
  double quad;
  double cpo_i;
  arma::colvec cpo_arr = arma::colvec(n);
  arma::colvec waic_arr = arma::colvec(n);
  double LPML = 0.;
  double prior_i_out, l_ij, ll_ij, ll_ij_z0_u0, ll_ij_z0_u1, ll_ij_z1_u0, ll_ij_z1_u1;
  arma::colvec prior_i_out_arr;
  arma::colvec ll_ij_arr;
  arma::colvec ll_i_arr = arma::colvec(n_sim);
  arma::colvec f_arr;
  arma::mat ll_00 = arma::zeros(n_sim,n);
  arma::mat ll_01 = arma::zeros(n_sim,n);
  arma::mat ll_10 = arma::zeros(n_sim,n);
  arma::mat ll_11 = arma::zeros(n_sim,n);
  arma::colvec ll = arma::zeros(n_sim);
  
  Progress pb(n, true);
  for (arma::uword i=0; i<n; i++) {
    irow = {i};
    fixed_eff_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_ind-1);
    rand_eff_i = (as<arma::mat>(fixed_eff[i])).cols(rand_eff_ind-1);
    fixed_eff_z_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_z_ind-1);
    fixed_eff_u_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_u_ind-1);
    fixed_eff_W_i = (as<arma::mat>(fixed_eff[i])).cols(fixed_eff_W_ind-1);
    W_ind_i = (as<arma::colvec>)(W_ind[i]);
    cpo_i = 0;
    
    for (arma::uword iter = 0; iter < n_sim; iter++) {
      l_ij = 0;
      pi_h = 1./(1. + std::exp(-arma::dot(fixed_eff_z_i.row(0), gamma_h.row(iter).t())));
      pi_u = 1./(1. + std::exp(-arma::dot(fixed_eff_u_i.row(0), gamma_u.row(iter).t())));
      pi_W = arma::exp((fixed_eff_W_i*gamma_W.row(iter).t())) /
        (1. + arma::exp((fixed_eff_W_i*gamma_W.row(iter).t())));
      
      int pos_max = std::pow(2, ni(i));
      
      Lambda_mat = arma::reshape(Lambda.row(iter).t(), d_rand, d_rand);
      if (is_conditional != 0) {
        W_i_vec = W_mat.row(iter).subvec(ni_max*i,ni_max*(i+1)-1).t();
        W_i = arma::diagmat(W_i_vec.subvec(0,ni(i)-1));
        I_W_i = arma::diagmat(arma::ones(ni(i))) - W_i;
        prior_i_out = 
          (1-u_mat(iter,i))*log(1-pi_u) + u_mat(iter,i)*log(pi_u) +
          (1-z_mat(iter,i))*log(1-pi_h) + z_mat(iter,i)*log(pi_h);
        for (int j=0; j < ni(i); j++) {
          prior_i_out +=
            (1-W_i_vec(j))*(1-pi_W(j)) + W_i_vec(j)*pi_W(j);
        }
        H_inv = arma::inv(
          std::pow(eta_u, 2.*u_mat(iter,i)) * (rand_eff_i * Lambda_mat * rand_eff_i.t()) +
            (eta_w*W_i + I_W_i) * h_mat(tt[i],rho(iter)) * (eta_w*W_i + I_W_i)
        );
        quad =
          ((as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()) *
             H_inv *
             (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()).t()).eval()(0);
        sigma2_arr = sigma2.row(iter).t();
        if (is_conditional == 1) {
          ll_i_arr(iter) = compute_pdf_II(
            quad, alpha.row(iter).t(), z_mat(iter,i), sigma2_arr, H_inv, ni(i)
          );
        } 
        // else if (is_conditional == 2) {
        //  // this is a workaround
        //  sigma2_arr(0) = sigma_i2_mat(iter,i);
        //  ll_i_arr(iter) = compute_pdf_II(
        //    quad, alpha.row(iter).t(), 0, sigma2_arr, H_inv, ni(i)
        //  );
        // }
        
      } else if (is_conditional == 0){
        ll_ij_arr = arma::colvec(pos_max).fill(-arma::datum::inf);
        prior_i_out_arr = arma::colvec(pos_max).fill(-arma::datum::inf);
        for (int pos = 0; pos < pos_max; pos++) {
          pos_bin = D2B(pos, ni(i));
          W_i = arma::diagmat(arma::conv_to<arma::colvec>::from(pos_bin));
          I_W_i = arma::diagmat(arma::ones(ni(i))) - W_i;
          
          prior_i_out = 1.;
          // when we allow outliers in the model
          if (mod_type == 3 || mod_type == 4 || mod_type == 5) {
            if (arma::accu(arma::find(W_ind_i == pos)) == 0) {
              continue;
            }
            // ll_ij_arr(pos) = 0;
            for (int j = 0; j < ni(i); j++) {
              prior_i_out *= 
                std::pow(pi_W(j), pos_bin(j)) * std::pow(1-pi_W(j), 1-pos_bin(j));
            }
          }
          
          // prior_i_out_marg += prior_i_out;
          prior_i_out_arr(pos) = std::log(prior_i_out);
          ll_ij_arr(pos) = prior_i_out_arr(pos);
          
          // 00
          H_inv = arma::inv(
            (rand_eff_i * Lambda_mat * rand_eff_i.t()) +
              (eta_w*W_i + I_W_i) * h_mat(tt[i],rho(iter)) * (eta_w*W_i + I_W_i)
          );
          quad =
            ((as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()) *
               H_inv *
               (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()).t()).eval()(0);
          
          zeta = 1./std::pow(alpha(iter,0), 2) - ni(i)/2.;
          ll_ij_z0_u0 = compute_pdf_II(
            quad, alpha.row(iter).t(), 0, sigma2.row(iter).t(), H_inv, ni(i)
          );
          
          // 01
          H_inv = arma::inv(
            std::pow(eta_u, 2.) * (rand_eff_i * Lambda_mat * rand_eff_i.t()) +
              (eta_w*W_i + I_W_i) * h_mat(tt[i],rho(iter)) * (eta_w*W_i + I_W_i)
          );
          quad =
            ((as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()) *
               H_inv *
               (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()).t()).eval()(0);
          zeta = 1./std::pow(alpha(iter,0), 2) - ni(i)/2.;
          ll_ij_z0_u1 = compute_pdf_II(
            quad, alpha.row(iter).t(), 0, sigma2.row(iter).t(), H_inv, ni(i)
          );
          
          // 10
          H_inv = arma::inv(
            (rand_eff_i * Lambda_mat * rand_eff_i.t()) +
              (eta_w*W_i + I_W_i) * h_mat(tt[i],rho(iter)) * (eta_w*W_i + I_W_i)
          );
          quad =
            ((as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()) *
               H_inv *
               (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()).t()).eval()(0);
          zeta = 1./std::pow(alpha(iter,1), 2) - ni(i)/2.;
          ll_ij_z1_u0 = compute_pdf_II(
            quad, alpha.row(iter).t(), 1, sigma2.row(iter).t(), H_inv, ni(i)
          );
          
          // 11
          H_inv = arma::inv(
            std::pow(eta_u, 2.) * (rand_eff_i * Lambda_mat * rand_eff_i.t()) +
              (eta_w*W_i + I_W_i) * h_mat(tt[i],rho(iter)) * (eta_w*W_i + I_W_i)
          );
          quad =
            ((as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()) *
               H_inv *
               (as<arma::rowvec>(y_mat[i]) - (fixed_eff_i * beta.row(iter).t()).t()).t()).eval()(0);
          zeta = 1./std::pow(alpha(iter,1), 2) - ni(i)/2.;
          ll_ij_z1_u1 = compute_pdf_II(
            quad, alpha.row(iter).t(), 1, sigma2.row(iter).t(), H_inv, ni(i)
          );
          
          if (mod_type == 1 || mod_type == 5) {
            ll_ij_arr(pos) += ll_ij_z0_u0;
          } else if (mod_type == 2 || mod_type == 3) {
            f_arr = arma::colvec(2);
            f_arr(0) = log(1-pi_u) + ll_ij_z0_u0;
            f_arr(1) = log(pi_u) + ll_ij_z0_u1;
            ll_ij_arr(pos) += sum_of_exp(f_arr);
          } else if (mod_type == 4) {
            f_arr = arma::colvec(4);
            f_arr(0) = log(pi_h) + log(pi_u) + ll_ij_z1_u1;
            f_arr(1) = log(pi_h) + log(1-pi_u) + ll_ij_z1_u0;
            f_arr(2) = log(1-pi_h) + log(pi_u) + ll_ij_z0_u1;
            f_arr(3) = log(1-pi_h) + log(1-pi_u) + ll_ij_z0_u0;
            ll_ij_arr(pos) += sum_of_exp(f_arr);
          }
          // Rcout << "f_arr: " << f_arr << "\n";
          if (mod_type == 1 || mod_type == 2) break;
        }
        ll_i_arr(iter) = sum_of_exp(ll_ij_arr) - sum_of_exp(prior_i_out_arr);
      }
      
      if (u_mat(iter,i) == 0 && z_mat(iter,i) == 0) {
        ll_00(iter,i) = ll_i_arr(iter);
      } else if (u_mat(iter,i) == 0 && z_mat(iter,i) == 1) {
        ll_01(iter,i) = ll_i_arr(iter);
      } else if (u_mat(iter,i) == 1 && z_mat(iter,i) == 0) {
        ll_10(iter,i) = ll_i_arr(iter);
      } else if (u_mat(iter,i) == 1 && z_mat(iter,i) == 1) {
        ll_11(iter,i) = ll_i_arr(iter);
      }
      
      ll(iter) += ll_i_arr(iter);
    }
    
    cpo_arr(i) = - sum_of_exp(-ll_i_arr) + std::log(n_sim);
    waic_arr(i) = - sum_of_exp(ll_i_arr) + std::log(n_sim) + 2*arma::mean(ll_i_arr);
    pb.increment();
  }
  LPML = arma::accu(cpo_arr);
  return (2*arma::accu(waic_arr));
}


