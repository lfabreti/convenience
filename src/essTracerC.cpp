#include <Rcpp.h>
using namespace Rcpp;

//' @export

// [[Rcpp::export]]
double essTracerC( NumericVector x ){
  
  double samples = x.size();
  double MAX_LAG = 2000;
  int maxLag = (samples - 1 < MAX_LAG ? samples - 1 : MAX_LAG);
  
  double m = 0;
  for (int i=0; i<samples; i++)
  {
    m += x.at(i);
  }
  double mean = m/samples;
  
  double* gammaStat = new double[maxLag];
  // setting values to 0
  for (size_t i=0; i<maxLag; i++)
  {
    gammaStat[i] = 0;
  }
  double varStat = 0.0;
  
  for (int lag = 0; lag < maxLag; lag++) {
    for (int j = 0; j < samples - lag; j++) {
      double del1 = x.at(j) - mean;
      double del2 = x.at(j + lag) - mean;
      gammaStat[lag] += (del1 * del2);
    }
    
    gammaStat[lag] /= ((double) (samples - lag));
    
    if (lag == 0) {
      varStat = gammaStat[0];
    } else if (lag % 2 == 0) {
      // fancy stopping criterion :)
      if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
        varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
      }
      // stop
      else
        maxLag = lag;
    }
  }
  
  double act = varStat / gammaStat[0];
  
  // effective sample size
  double ess = samples / act;
  
  return ess;
}
