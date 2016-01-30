#include <Rcpp.h>
using namespace Rcpp;

//’ Calculated max entropy solution C++ code
//’
//’ @param aw NumericMatrix Initial matrix with the zero elements already fixed.
//’ @param rs NumericVector Row sums.
//’ @param rs NumericVector Column sums.
//’ @param minError double minimum error for the iterator.
//’ @return NumericMatrix with the solution or the result after 5000 iterations.
//' @export
// [[Rcpp::export]]
NumericMatrix maxEntropyCpp(NumericMatrix aw, NumericVector rs, NumericVector cs,
                            double minError) {
  
  double err = 1;

  int n = rs.length();
  // These vectors hold the current row/col sums
  NumericVector rsaw(n);
  NumericVector csaw(n);
  // These vectors tell us whether any out/in edges are allowed
  LogicalVector rson = rs > 1e-12;
  LogicalVector cson = cs > 1e-12;
  // We also need to catch missing values
  LogicalVector rsna = is_na(rs);
  LogicalVector csna = is_na(cs);
  // Remember negatives are considered missing values
  for (int i=0; i<n; i++) {
    if(rs[i] < -1e-12) rsna[i] = true;
    if(cs[i] < -1e-12) csna[i] = true;
  }
  
  // Clear zero row/cols
  for(int i=0; i<n; i++) {
    for (int j=0;j<n;j++) {
      if((!rson[i]) || (!cson[j]))
           aw(i,j) = 0;
    }
  }
  
  int iter = 0;
  while(err > minError && iter < 5000) {
    // Get the row sums
    for(int i=0; i<n; i++) {
      rsaw[i]=0;
        for (int j=0;j<n;j++) {
          rsaw[i] += aw(i,j);
        }
    }
    
    // Apply row sums
    for (int i=0;i<n;i++) {
      if(rson[i] && !rsna[i]) {
        for (int j=0;j<n;j++) {
          aw(i,j) *= rs[i] / rsaw[i];
        }
      }
    }
    
    // Now for columns
    // Get the col sums
    for(int i=0; i<n; i++) {
      csaw[i]=0;
        for (int j=0;j<n;j++) {
          csaw[i] += aw(j,i);
        }
    }
    
    // Apply col sums (if appropriate)
    for (int i=0;i<n;i++) {
      if(cson[i] && !csna[i]) {
        for (int j=0;j<n;j++) {
          aw(j,i) *= cs[i] / csaw[i];
        }
      }
    }

    // Calculate the error, ignore missings    
    err=0;
    for (int i=0; i<n; i++) {
      if(!rsna[i]) {
        err += std::pow(rs[i] - rsaw[i],2);
      }
      if(!csna[i]) {
        err += std::pow(cs[i] - csaw[i],2);
      }
    }
        
    iter ++;
  }


  return(aw);  
}

