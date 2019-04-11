// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

//' The IT alrogirhm from Ulrich & Gotelli (2010)
//'
//' This function does blah blah blah.
//'
// [[Rcpp::export]]
IntegerMatrix IT (NumericVector rowsums,
                  NumericVector colsums,
                  int nrow,
                  int ncol,
                  int N)
{
  IntegerVector rowSeq = seq_len(nrow);
  IntegerVector colSeq = seq_len(ncol);
  IntegerMatrix Y(nrow, ncol);
  int size = 1;

  IntegerVector rows = RcppArmadillo::sample(rowSeq, N, TRUE, rowsums);
  IntegerVector cols = RcppArmadillo::sample(colSeq, N, TRUE, colsums);

  for(int i = 0; i < N; ++i)
  {
    int row = rows(i) - 1;
    int col = cols(i) - 1;
    Y(row, col) = Y(row, col) + 1;
  }

  return Y;
}
