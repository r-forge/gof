#ifndef CUMRES_H
#define CUMRES_H

#include <cmath>
#include "matrix.h" 


inline double ss2(const scythe::Matrix<double> &M) {
  return(std::sqrt(scythe::sum(M%M)));
}

double KolmogorovSmirnov(const scythe::Matrix<double> &W);

double CramerVonMises(const scythe::Matrix<double> &W, const scythe::Matrix<double> &x);

scythe::Matrix<double> Kz(const scythe::Matrix<double> &x, const double b);

scythe::Matrix<double> Wi(const scythe::Matrix<double> &r, 
			  const scythe::Matrix<double> &x,
			  const scythe::Matrix<double> &betaiid,
			  const scythe::Matrix<double> &etaraw, const double b, 
			  scythe::Matrix<double> &Wobs);

scythe::Matrix<double> Cpred(const scythe::Matrix<double> &cum, 
				    const scythe::Matrix<double> &x);


#endif /* CUMRES_H */
