#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix maxEntropyCpp(NumericMatrix aw, NumericVector rs, NumericVector cs) {
  
  double err = 1;
  int i=0;

  int n = rs.length();
  NumericVector rsaw(n);
  NumericVector csaw(n);
  
  LogicalVector rson(n, true);
  LogicalVector cson(n, true);
  
  for (int i=0; i<n; i++) {
    if(rs[i] > 1e-12) rson[i] = false;
    if(cs[i] > 1e-12) cson[i] = false;
  }
    
  while(err > 1e-18) {
    // Get the row sums
    for(int i=0; i<n; i++) {
      rsaw[i]=0;
      if (rson[i]) {
        for (int j=0;j<n;j++) {
          rsaw[i] += aw(i,j);
        }
      } else {
        rsaw[i]=1;
      }
    }
    
    // Apply row sums
    for (int i=0;i<n;i++) {
      for (int j=0;j<n;j++) {
        aw(i,j) *= rs[i] / rsaw[i];
      }
    }
    
    // Now for columns
    // Get the col sums
    for(int i=0; i<n; i++) {
      csaw[i]=0;
      if (cson[i]) {
        for (int j=0;j<n;j++) {
          csaw[i] += aw(j,i);
        }
      } else {
        csaw[i]=1;
      }
    }
    
    // Apply col sums
    for (int i=0;i<n;i++) {
      for (int j=0;j<n;j++) {
        aw(j,i) *= cs[i] / csaw[i];
      }
    }
    
    
    
  }
  
}

