// STL
#include <iostream>
#include <algorithm> // random_shuffle, reverse, sort, ...
#include <cmath>
// SCYTHE
#include "matrix.h" 
#include "distributions.h"
#include "ide.h" 
#include "la.h"
#include "mersenne.h"
#include "rng.h"
#include "stat.h" 
#include "smath.h" 
// R interface
#include <R.h>           //  Rprintf()
#include <R_ext/Utils.h> //  user interrupts
#include <Rdefines.h>
#include <Rinternals.h>
// 
#include "extra.h"
#include "cumres.h"


//*************************************************


using namespace scythe;
using namespace std;

double KolmogorovSmirnov(const scythe::Matrix<double> &W) {
  return(scythe::max(scythe::fabs(W))); // sup |W(x)|

}

double CramerVonMises(const scythe::Matrix<double> &W, const scythe::Matrix<double> &x) {
  unsigned n = x.size();
  scythe::Matrix<double> delta(W.rows(),W.cols());
  for (unsigned i=0; i<(n-1); i++) {
    delta[i] = x[i+1]-x[i];       
  } // OBS: last element in delta is zero
  return(scythe::sum(delta%W%W)); // int |W(x)|^2 dx
}

Matrix<double> Kz(const Matrix<double> &x, const double b) {
  unsigned n = x.size();
  Matrix<double> xo = sort(x);
  Matrix <double, Col> ones(n, n);
  for (unsigned i=0; i<n; i++) {
    for (unsigned j=0; j<n; j++) {
      if (b>0) {
	if (x[j]<= xo[i] && x[j]> xo[i]-b) ones(i,j) = 1;	
      } else {
	if (x[j]<= xo[i]) ones(i,j) = 1;	
      }
    }
  }
  return(ones);
}

Matrix<double> Cpred(const Matrix<double> &cum, 
		     const Matrix<double> &x) {
  double timex,sc1,sc2;
  unsigned n = cum.rows();
  unsigned p = cum.cols();
  unsigned nx = x.rows();
  Matrix<double> pred(nx,p);  
  double smax = cum(n-1,0);
  for (unsigned s=0; s<nx; s++) {
    timex=x[s]; pred(s,0)=timex; double c=n-1;
    sc1=cum(n-1,0); sc2=smax+x[nx-1];    
    while ((!((timex<sc2) && (timex>=sc1))) && (c>0)) {
      sc1=cum[c-1]; sc2=cum[c]; c=c-1; }
    for(unsigned j=1; j<p; j++) pred(s,j) = cum(c,j);
  }
  return(pred);
}


Matrix<double> Wi(const Matrix<double> &r, 
		  const Matrix<double> &x,
		  const Matrix<double> &betaiid,
		   const Matrix<double> &etaraw, const double b,
		  Matrix<double> &Wobs) {
  unsigned n = r.rows();
  
  Matrix<double> ones = Kz(x,b);

  Matrix<double> rs = ones;
   for (unsigned i=0; i<n; i++) {
    rs(i,_) = t(r) % ones(i,_);
  }
//   cerr << "ones=" << ones << endl;
//   cerr << "r=" << r << endl;
  Wobs = ones*r/sqrt(n); // We expect Wobs at call to contain column
//   cerr << "Wobs=" << ones*r << endl;
  //  cerr << "r=" << r << endl;
//   cerr << "n=" << n << endl;
  // matrix of same length as residuals
  Matrix<double> eta = ones*etaraw;
  Matrix<double> W1 = eta*betaiid;
  //  cerr << "rs=" << rs << endl;
  //  cerr << "W1=" << W1 << endl;
  
  Matrix<double> res = (rs-W1);
  //  for (unsigned i=0; i<n; i++) {
  //    res(i,_) = t(res(i,_))%r;
  //  }
  return(res/sqrt(n));
}

    
extern "C" {
  void W(const unsigned *R,
	 const double *b,
	 const unsigned *n, 
	 const unsigned *npar,
	 const double *xdata,
	 const double *rdata, 
	 const double *betaiiddata, 
	 const double *etarawdata,
	 const unsigned *plotnum,
	 const int *seed,
 	 double *KS,
	 double *CvM,
	 double *Wsd,
	 double *cvalues,
	 double *Ws,
	 double *Wobsdata
	 ) {

    Matrix<double> x(*n, 1, xdata);
    Matrix<double, Col> etaraw(*n, *npar, etarawdata);    
    Matrix<double> r(*n, 1, rdata);
    Matrix<double, Col> betaiid(*npar, *n, betaiiddata);   

    Matrix<double> xo = sort(x);
    Matrix<double> Wobs = r;

//     cerr << "betaiid=" << betaiid << endl;


    Matrix<double> W = Wi(r, x, betaiid, etaraw, *b, Wobs);
    //     cerr << "W=" << W << endl;
//    return(0);
    
    Matrix<double> sdW = apply(W,1,ss2);

    double KSobs = KolmogorovSmirnov(Wobs);
    double CvMobs = CramerVonMises(Wobs,xo); //xo

    unsigned KScount=0; 
    unsigned CvMcount=0;
    mersenne myrng; myrng.initialize((unsigned long)*seed);
    Matrix<double> N(*n,1); 
    for (unsigned j=0; j<(*n); j++) N[j]=myrng.rnorm(0, 1);

    //    cerr << "Starting simulation\n";
    //    cerr << "W=" << W <<"\n";
    Matrix<double> Res(min((double)*plotnum,(double)*R),*n);
    for (unsigned i=0; i< *R; i++) {
      for (unsigned j=0; j<(*n); j++) N[j]=myrng.rnorm(0, 1);
      Matrix<double> What = W*(N);
      //      cerr << "What=" << What << endl << endl << endl;
      cvalues[i] = max(fabs(What/sdW));
      double KShat = KolmogorovSmirnov(What);
      double CvMhat = CramerVonMises(What,xo); //xo
      //      double KShat = 0, CvMhat = 0;
      //cerr << "!\n";
      if (KShat>KSobs) KScount++;
      if (CvMhat>CvMobs) CvMcount++;
      if ((unsigned)i< *plotnum) {
	Res(i,_) = What;
      }
    }
    //    cerr << "Done\n";
    
    for (unsigned i=0; i< *n; i++) {
      Wsd[i] = sdW[i];
      Wobsdata[i] = Wobs[i];
    }    

    KS[0] = (double)KScount/(double)(*R);
    CvM[0] = (double)CvMcount/(double)(*R);
    for (unsigned r=0; r< Res.rows(); r++) {
      for (unsigned s=0; s< Res.cols(); s++) {	
	unsigned pos = r*(Res.cols())+s;
	Ws[pos] = Res(r,s);
      }
    }
    
  }
} // extern "C"

