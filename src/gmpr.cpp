#include <Rcpp.h>

//[[Rcpp::plugins(cpp11)]]

#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;
using namespace Rcpp;

class GMPR {
private:
  const IntegerMatrix comm;
  const int n, p;
  const int minct;
  const int inter_n;
  void diag(vector<double> &square_matrix, const int &n, const double &i);
public:
  vector<double> factor;
  vector<double> size_factor;
  vector<int> NSS;
  GMPR(const IntegerMatrix &otumatrix, int nrow, int pcol,
       int min_ct, int intersect_no);
  GMPR(const IntegerMatrix &otumatrix, int nrow, int pcol);
  ~GMPR();
  void Size_factor(void);
  void Factor(void);
};

GMPR::~GMPR() {};

GMPR::GMPR(const IntegerMatrix &otumatrix, int nrow, int pcol, int min_ct, int intersect_no) :
  comm(otumatrix), n(nrow), p(pcol), minct(min_ct), inter_n(intersect_no), factor(n*n, 0),
  size_factor(n, 0), NSS(n, 0)
{};


GMPR::GMPR(const IntegerMatrix &otumatrix, int nrow, int pcol) :
  comm(otumatrix), n(nrow), p(pcol), minct(2), inter_n(4), factor(n*n, 0),
  size_factor(n, 0), NSS(n, 0)
{};
/*
void GMPR::Size_factor(void) {
for (int i = 0; i < n; i++) {
for (int j = 0; j < n; j++) {
if (abs(factor[i*n + j]) > 1e-10) {
NSS[i] += 1;
size_factor[i] *= factor[i*n + j];
};
};
NSS[i] -= 1;  //substrat itself
if (NSS[i] >= inter_n)
size_factor[i] = pow(size_factor[i], 1.0 / NSS[i]);
};
};
*/
void GMPR::Size_factor(void) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (abs(factor[i*n + j]) > 1e-10) {
        NSS[i] += 1;
        size_factor[i] += log(factor[i*n + j]);
      };
    };
//    NSS[i] -= 1;  //substrat itself
    size_factor[i] = exp(size_factor[i] / NSS[i]);
  };
};


void GMPR::Factor(void) {

  vector<bool> index(p*n, false);
  void diag(vector<double> &square_matrix, const int &n, const double &i);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      index[j + i*p] = comm(i,j) >= minct;
    };
  };

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      int h = 0;
      vector<double> ratio(p, 0);
      for (int k = 0; k < p; k++) {
        if (index[i*p + k] && index[j*p + k]) {
          ratio[h] = static_cast<double>(comm(i,k)) / comm(j,k);
          h++;
        };
      };
      if (h >= inter_n) {
        sort(ratio.begin(), ratio.begin()+h);
        if (h % 2 == 1) {
          factor[i*n + j] = ratio[h / 2];
          factor[j*n + i] = 1 / (ratio[h / 2]);
        }
        else {
          factor[i*n + j] = (ratio[h / 2 - 1] + ratio[h / 2]) / 2;
          factor[j*n + i] = ( (1/ratio[h / 2 - 1]) + (1/ratio[h / 2]) ) / 2;
        }
      };
    };
  };
  GMPR::diag(factor, n, 1);
}

void GMPR::diag(vector<double> &square_matrix, const int &n, const double &i) {
  for (int j = 0; j < n; j++) square_matrix[n*j + j] = i;
};

// [[Rcpp::export]]
NumericVector gmpr(IntegerMatrix x, int min_ct, int intersect_no){
  /*vector<vector<int> > table { { 0,0,0,0,2,0,11,0,0,1,0,0,11,0,0,0,126,0,4,4,12,1,0,0,0,12,1,8,66,183,8,5,3,5,27,2,0,1,0,20,0,1,2,2,2,76,16,2,2,31 },{ 0,0,0,0,0,0,1,0,1,1,0,1,1,0,0,1,12,0,2,4,13,0,0,0,3,1,1,11,17,41,1,5,0,1,2,2,0,0,0,40,0,3,3,1,4,7,3,11,0,0 },
                                 { 0,2,0,0,0,0,9,0,2,1,0,0,1,0,1,0,88,0,2,1,5,0,0,0,0,82,1,7,6,21,2,2,3,0,3,16,0,2,0,3,0,2,3,0,1,10,15,3,0,2 },{ 1,2,1,0,2,0,10,0,2,1,0,0,7,0,0,0,8,0,1,1,4,0,0,0,0,4,0,0,2,12,2,1,1,0,0,0,0,1,0,10,0,1,0,0,0,2,7,1,0,0 },
  { 2,3,8,0,0,0,6,1,0,1,0,0,16,1,0,1,6,0,0,5,12,0,0,0,1,1,0,4,0,3,3,0,2,0,1,0,0,3,0,2,0,0,3,1,1,2,9,1,0,1 },{ 3,1,2,1,0,0,3,1,0,2,1,1,18,0,0,2,85,2,4,6,12,0,0,0,4,19,0,16,20,83,6,0,2,38,28,10,0,1,3,82,5,3,21,0,17,749,137,14,0,74 },
  { 0,3,34,0,0,0,1,0,1,0,0,0,12,0,0,0,30,0,0,0,6,0,0,0,1,1,1,0,8,60,6,0,3,1,1,5,0,0,0,9,5,0,1,0,5,10,14,0,0,1 },{ 0,0,2,0,0,0,0,1,0,0,1,0,2,0,0,0,82,6,22,15,54,1,0,0,4,10,0,7,1,40,19,6,18,8,66,9,0,3,0,92,0,0,6,2,2,10,10,1,1,33 },
  { 0,0,0,0,0,0,1,0,0,0,0,1,8,0,0,0,61,3,8,6,7,3,1,0,0,7,0,4,0,16,6,1,5,2,6,7,0,0,120,53,0,0,9,0,3,27,6,5,1,0 },{ 0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,714,0,103,0,18,0,0,0,1,75,16,2,90,219,1,4,3,4,46,116,0,35,0,59,0,1,2,2,3,19,4,1,1,3 },
  { 0,0,0,0,2,0,19,7,3,2,0,0,1,0,0,0,62,0,0,0,3,0,0,0,0,3,0,0,7,63,0,2,0,22,35,4,0,0,0,11,0,3,0,0,0,21,1,0,0,3 },{ 0,2,0,0,0,0,1,0,0,0,0,0,2,0,0,1,21,0,157,3,7,0,0,0,3,26,9,2,18,47,2,3,1,1,13,13,0,0,0,45,0,0,0,0,6,4,1,1,0,1 },
  { 0,1,0,0,0,0,2,0,0,0,0,0,0,0,0,0,82,0,2,0,5,0,0,0,0,116,1,1,4,15,3,0,0,0,6,13,0,0,0,1,0,0,0,0,0,2,2,0,0,0 },{ 8,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,9,1,0,0,4,0,0,0,2,0,0,1,1,7,2,0,1,0,1,0,0,0,0,8,0,0,2,0,0,5,1,0,0,0 },
  { 0,0,1,0,0,0,0,0,0,0,0,0,2,0,0,0,116,0,7,19,21,0,0,0,6,72,0,14,20,86,5,9,4,3,133,18,0,5,0,62,0,3,0,0,3,36,16,0,0,2 },{ 2,0,7,2,0,0,1,0,0,0,1,0,69,0,0,1,3,0,0,1,9,0,0,0,1,0,0,2,0,10,1,0,1,0,1,4,0,24,2,5,1,0,1,0,1,1,4,1,1,0 },
  { 0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,197,5,35,40,32,12,8,0,0,11,2,12,0,78,5,13,6,2,20,10,0,1,48,63,0,0,16,1,4,76,2,2,1,1 },{ 0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,97,2,0,12,35,15,0,0,0,1,0,41,1,78,27,1,12,35,41,7,0,8,0,31,0,2,1,0,2,16,60,9,0,6 },
  { 1,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,4,1,0,2,2,0,0,0,0,2,0,4,9,39,1,0,0,5,4,1,0,0,0,66,1,0,0,0,0,9,4,1,0,1 },{ 3,0,10,1,0,0,2,0,0,0,0,1,6,0,0,0,55,0,0,2,9,0,0,0,0,7,1,1,9,73,1,0,7,1,8,6,0,2,2,53,8,1,0,0,2,9,10,4,0,0 },
  { 0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,119,8,19,16,56,1,1,0,1,12,0,0,2,49,3,6,3,29,93,15,0,12,0,17,0,0,2,1,3,3,3,3,2,49 },{ 0,1,0,0,0,0,2,1,1,0,0,1,2,0,0,4,145,1,3,9,67,2,0,0,40,11,1,30,62,136,6,13,10,4,9,13,0,3,1,57,0,0,15,6,34,141,28,116,0,7 },
  { 0,0,0,1,0,0,4,0,0,1,0,0,1,0,0,2,69,2,0,5,29,1,0,0,1,12,0,13,48,125,15,15,6,5,12,6,0,4,0,74,2,1,12,1,16,178,23,35,2,20 },{ 0,0,0,0,0,0,1,1,0,1,0,0,2,0,0,0,68,2,3,6,87,0,0,0,20,9,1,32,38,72,9,32,5,4,16,5,0,5,0,76,0,0,4,12,26,30,3,45,0,10 },
  { 0,1,0,0,0,0,5,0,1,0,0,0,1,0,0,0,253,1,5,43,180,0,0,0,38,39,0,32,59,140,37,56,6,3,14,57,0,2,1,34,0,0,60,0,12,91,71,63,4,18 },{ 0,2,1,0,3,0,3,0,0,0,0,0,18,0,0,0,51,0,2,2,8,0,0,0,1,2,0,8,11,33,8,2,9,0,1,1,0,2,12,21,0,1,9,0,1,30,72,20,0,1 },
  { 0,0,0,8,0,0,3,1,0,0,0,2,18,0,2,7,68,1,1,40,99,1,0,0,30,9,0,21,44,76,15,15,14,1,12,13,0,0,13,64,1,0,16,2,12,130,63,14,3,15 },{ 1,0,1,2,1,0,1,0,0,0,3,0,9,1,1,5,35,1,1,7,61,1,0,1,9,3,0,10,2,10,25,9,17,1,1,2,0,5,8,5,0,2,17,0,9,20,25,13,0,19 },
  { 0,0,1,2,1,0,0,0,0,1,0,0,1,0,0,2,135,4,1,25,128,2,0,2,28,6,0,55,7,23,21,13,21,8,5,1,0,2,8,6,0,1,49,2,23,227,89,148,6,9 },{ 0,0,2,2,0,0,0,0,0,0,1,2,7,0,0,1,76,0,0,8,53,2,0,0,17,9,3,18,14,33,22,6,5,5,3,7,0,3,12,39,14,0,55,0,16,379,208,27,1,34 } };
  */
  GMPR A(x, x.nrow(), x.ncol(), min_ct, intersect_no);
  A.Factor();
  A.Size_factor();
  return Rcpp::wrap(A.size_factor);
  }
