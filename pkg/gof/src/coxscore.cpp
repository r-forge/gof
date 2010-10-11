// STL
#include <iostream>
#include <algorithm> // random_shuffle, reverse, sort, ...
#include <cmath>
// SCYTHE
#include "mersenne.h"
#include "rng.h"
#include "distributions.h"
#include "ide.h" 
#include "la.h"
#include "matrix.h" 
#include "stat.h" 
#include "smath.h" 
// R interface
#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts
#include <Rdefines.h>
#include <Rinternals.h>
#include "cumres.h"
#include "extra.h"

using namespace scythe;
using namespace std;

Matrix<double> Wscorerate_cox(unsigned Var,
			      const Matrix<double> &X, 
			      const Matrix<double> &schoen, 
			      const Matrix<double> &RR, 
			      const Matrix<double> &E,
			      const Matrix<double> &S_0,
			      const Matrix<double> &cumhaz,
			      const Matrix<double> &beta_iid,
			      const Matrix<double> &It,			  
			      const Matrix<unsigned> &index_dtimes,
			      const Matrix<double> &time) {
  Matrix<double> dtimes = chrows(time, index_dtimes);  
  unsigned n=X.rows(); unsigned p=X.cols(); unsigned nd=dtimes.size();

  Matrix<double> schoendN(n,p); // Initialized as zero
  for (unsigned i=0; i<nd; i++) 
    schoendN(index_dtimes[i],_) = schoen(i,_);
  
  Matrix<double> dMscorerate(n,nd); // component 'Var' of score
  for (unsigned i=0; i<nd; i++) {      
    Matrix<double> risk =  (time>=dtimes[i]);

    dMscorerate(_,i) = -risk%RR%(X(_,Var)-E(i,Var))/S_0[i];
    dMscorerate(index_dtimes[i],i) = schoen(i,Var)+dMscorerate(index_dtimes[i],i);
  }
  Matrix<double> Mscorerate = (cumsum(t(dMscorerate)));
  
  Matrix<double> Wscorerate(n,nd); 
  for (unsigned i=0; i<n; i++) {
    for (unsigned j=0; j<nd; j++) {
      Matrix<double> Itd = It(j,_); Itd.resize(p,p);
      Matrix<double> betaiidrow = t(beta_iid(i,_));
      //      Wscorerate(i,j) = Mscorerate(i,j) - (Itd*betaiidrow)[Var];
      Wscorerate(i,j) = Mscorerate(j,i) - (Itd*betaiidrow)[Var];
    }
  }   


  return(Wscorerate);
}

extern "C" {
  void coxscoreW(const int *R, // Number of realizations
		 const int *n, // Number of individuals
		 const int *nd, // Number of death
		 const int *p, // Number of  parameters
		 const double *beta_data,  // nxp, parameter vector
		 const double *time_data, // 
		 const unsigned *index_dtimes_data, // Death times
		 const double *X_data,  // nxp, Design matrix
		 const double *beta_iid_data, // nxp 
		 const double *Mt_data, // Martingale residuals
		 const unsigned *paridx, const int *nparidx,
		 const unsigned *Type, 
		 const unsigned *seed,
		 const unsigned *plotnum,
		 double *KS,
		 double *CvM,
		 double *Wsd,
		 double *cvalues,
		 double *Ws,
		 double *W,
		 double *WWW // Not used;
		 ) {

    Matrix<double, Col> X(*n, *p, X_data);
    Matrix<double, Col> beta_iid(*n, *p, beta_iid_data);
    Matrix<double> beta(*p, 1, beta_data);
    Matrix<double> time(*n, 1, time_data);     
    Matrix<unsigned> index_dtimes(*nd, 1, index_dtimes_data); 
    Matrix<double> dtimes = chrows(time, index_dtimes);
    Matrix<unsigned> index_times = seqa(*n-1,-1,*n); // n:1 - 1  

    Matrix<double> xbeta = X*beta;
    Matrix<double> RR = exp(xbeta);
    Matrix<double> S_0 = chrows(reverse(cumsum(reverse(RR))), index_dtimes);
    //    cerr << "###" << endl;
    Matrix<double> cumhaz = cumsum(1/S_0);
    //    cerr << "###" << endl;
    Matrix<double> cumhaztime = Cpred(cbind(dtimes,cumhaz), time);
    //    cerr << "..." << endl;
    Matrix<double> XRR(*n,*p);
    unsigned p2 = (*p)*(*p);
    Matrix<double> XRRX(*n,p2);
    for (unsigned i=0; i<*n; i++) { 
      XRR(i,_) = RR[i]*X(i,_); 
      // Matrix<double> newr = RR[i]*crossprod(X(i,_)); newr.resize(1,p2);
      XRRX(i,_) = RR[i]*crossprod(X(i,_));
    }
    //    cerr << "..." << endl;
    Matrix<double> XRR_ = chrows(XRR, index_times);
    Matrix<double> XRRX_ = chrows(XRRX, index_times);
    Matrix<double> S_1  = chrows(chrows(cumsum(XRR_),index_times), index_dtimes); // Derivative of score,  n x p ## Martinussen & Scheike p. 182   
    Matrix<double> S_2 = chrows(chrows(cumsum(XRRX_),index_times), index_dtimes); // Second derivative, n x p^2 
    Matrix<double> E = multCol(S_1,1/S_0);
    Matrix<double> intS1oS02 = multCol(E, 1/S_0);
    intS1oS02 = cumsum(intS1oS02); // nd x p
    //    cerr << "..." << endl;
    Matrix<double> intS1oS02time = Cpred(cbind(dtimes,intS1oS02),time);
    Matrix<double> E_2(*nd,p2,false);
    for (unsigned i=0; i<*nd; i++) { 
      E_2(i,_) = crossprod(E(i,_)); // E(t,beta)^{\otimes 2} p. 184
    }
    Matrix<double> It = multCol(S_2, 1/S_0); 
    It = cumsum(It-E_2); // Martinussen & Scheike p. 184
    Matrix<double> Itau = It(*nd-1,_); Itau.resize(2,2); Matrix<double> Iinv = inv(Itau);
    Matrix<double> schoen = chrows(X,index_dtimes) - E; // Schoenfeld residuals

    Matrix<double> WW;
    if (*Type>1) {
      // Martingale residuals not yet implemented
    } else {        
      WW = Wscorerate_cox(paridx[0],X,schoen,RR,E,S_0,cumhaz,beta_iid,It,index_dtimes,time);
    }

    Matrix<double> sdW = apply(WW,2,ss2);    
    Matrix<double> Score = cumsum(schoen)(_,paridx[0]);
    double KSobs = KolmogorovSmirnov(Score);
    double CvMobs = CramerVonMises(Score,dtimes);

    unsigned KScount=0; 
    unsigned CvMcount=0;
    Matrix<double> Res(min((double)*plotnum,(double)*R),*nd);
    mersenne myrng; myrng.initialize((unsigned long) *seed);
    for (unsigned i=0; i<*R; i++) {
      Matrix<double> G(1,*n,false);
      for (unsigned j=0; j<G.size(); j++) G[j]=myrng.rnorm(0, 1);
      Matrix<double> ScoreSim = G*WW; 
      cvalues[i] = max(fabs(ScoreSim/t(sdW)));
      double KShat = KolmogorovSmirnov(ScoreSim);
      double CvMhat = CramerVonMises(ScoreSim,dtimes);
      if (KShat>KSobs) KScount++;
      if (CvMhat>CvMobs) CvMcount++;
      if ((unsigned)i< *plotnum) { Res(i,_) = ScoreSim; }
    }
    KS[0] = (double)KScount/(double)(*R);
    CvM[0] = (double)CvMcount/(double)(*R);

    for (unsigned s=0; s< Res.cols(); s++) {	
      for (unsigned r=0; r<Res.rows(); r++) {
	unsigned pos = r*(Res.cols())+s;
	Ws[pos] = Res(r,s);
      }
    }

    for (unsigned i=0; i< *nd; i++) { 
      Wsd[i] = sdW[i];
    }    
    for (unsigned r=0; r< *nd; r++) 
      W[r] = Score[r];
    

  }
} // extern "C"
