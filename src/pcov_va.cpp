#include <Rcpp.h>
using namespace Rcpp;


NumericMatrix mmult(const NumericMatrix& a, const NumericMatrix& b){
  if (a.ncol() != b.nrow()) stop ("Incompatible matrix dimensions");
  NumericMatrix out(a.nrow(),b.ncol());
  NumericVector rm1, cm2;
  for (int i = 0; i < a.nrow(); ++i) {
    rm1 = a(i,_);
    for (int j = 0; j < b.ncol(); ++j) {
      cm2 = b(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);
    }
  }
  return out;
}


NumericMatrix mmult_ele(const NumericMatrix& a, const NumericMatrix& b){
  if (a.ncol() != b.ncol()) stop ("Incompatible matrix dimensions");
  if (a.nrow() != b.nrow()) stop ("Incompatible matrix dimensions");
  NumericMatrix out(a.nrow(),b.ncol());
  for (int i = 0; i < a.nrow(); ++i) {
    for (int j = 0; j < b.ncol(); ++j) {
      out(i,j) = a(i,j)*b(i,j);
    }
  }
  return out;
}


NumericMatrix mtrans(const NumericMatrix& a){
  NumericMatrix out(a.ncol(),a.nrow());
  for (int i = 0; i < a.nrow(); ++i) {
    for (int j = 0; j < a.ncol(); ++j) {
      out(j,i) = a(i,j);
    }
  }
  return out;
}


NumericMatrix extract_mat(const NumericMatrix& a, int rr){
  NumericMatrix out(a.nrow()-1,a.ncol());
  for (int i = 0; i < rr; ++i) {
    for (int j = 0; j < a.ncol(); ++j) {
      out(i,j) = a(i,j)-a(rr,j);
    }
  }
  for (int i = a.nrow()-1; i > rr; --i) {
    for (int j = 0; j < a.ncol(); ++j) {
      out(i-1,j) = a(i,j)-a(rr,j);
    }
  }
  return out;
}


NumericMatrix acos_mat(const NumericMatrix& a){
  NumericMatrix out(a.nrow(),a.ncol());
  for (int i = 0; i < a.nrow(); ++i) {
    for (int j = 0; j < a.ncol(); ++j) {
      if(a(i,j) > 1){
        out(i,j) = acos(1);
      }else if(a(i,j) < -1){
        out(i,j) = acos(-1);
      }else{
        out(i,j) = acos(a(i,j));
      }
    }
  }
  return out;
}


NumericVector get_norm(const NumericMatrix& a)
{
  NumericVector rm(a.nrow());
  double ss;

  for (int i = 0; i < a.nrow(); i++) {
    ss = 0;
    for (int j = 0; j < a.ncol(); j++) {
      ss += std::pow(a(i,j),2);
    }
    rm(i) = std::sqrt(ss);
  }

  return rm;
}


NumericMatrix norm_mat(const NumericMatrix& a, NumericVector b)
{
  NumericMatrix nm = Rcpp::clone(a);

  for (int i = 0; i < a.nrow(); i++) {
    if(b(i) > 0){
      for (int j = 0; j < a.ncol(); j++) {
        nm(i,j) = a(i,j) / b(i);
      }
    }
  }
  return nm;
}


NumericMatrix zero_mat(const NumericMatrix& a, NumericVector b){
  NumericMatrix out = Rcpp::clone(a);

  for (int i = 0; i < a.nrow(); ++i) {
    if(b[i]==0){
      for(int j = 0; j < a.ncol(); ++j){
        out(i,j) = 0;
        out(j,i) = 0;
      }
    }
  }
  return out;
}


double sum_offd_multiply(const NumericMatrix& a)
{
  double sum = 0;
  for (int i = 0; i < a.nrow(); i++) {
    for (int j = (i+1); j < a.nrow(); j++) {
      sum += a(i,i) * a(i,j);
    }
  }
  for (int i = 1; i < a.nrow(); i++) {
    for (int j = 0; j < i; j++) {
      sum += a(i,i) * a(i,j);
    }
  }
  return sum;
}


double sum_offd(const NumericMatrix& a)
{
  double sum=0;
  for (int i = 0; i < a.nrow(); i++) {
    for (int j = (i+1); j < a.nrow(); j++) {
      sum += a(i,j);
    }
  }
  for (int i = 1; i < a.nrow(); i++) {
    for (int j = 0; j < i; j++) {
      sum += a(i,j);
    }
  }
  return sum;
}


NumericMatrix subtractr(const NumericMatrix& a, int rr)
{
  NumericVector b(a.nrow());
  NumericMatrix a_ex(a.nrow()-1, a.ncol());
  NumericMatrix a_exnorm(a.nrow()-1, a.ncol());
  NumericMatrix a_exnormt(a.ncol(), a.nrow()-1);
  NumericMatrix a_square(a.nrow()-1, a.nrow()-1);
  NumericMatrix a_acos(a.nrow()-1, a.nrow()-1);
  NumericMatrix a_fillna(a.nrow()-1, a.nrow()-1);

  a_ex = extract_mat(a, rr);
  b = get_norm(a_ex);
  a_exnorm = norm_mat(a_ex, b);
  a_exnormt = mtrans(a_exnorm);
  a_square = mmult(a_exnorm, a_exnormt);
  a_acos = acos_mat(a_square);
  a_fillna = zero_mat(a_acos, b);

  return a_fillna;
}


NumericVector insa(const NumericMatrix& matx, const NumericMatrix& maty)
{
  NumericMatrix matxym(matx.nrow(), maty.ncol());
  NumericMatrix matxye(matx.nrow(), maty.ncol());
  NumericVector out(2);

  double sr, matxysum, xsum, ysum;
  double matxyesum, matxsum, matysum;
  double sa_value, in_value;
  int r = matx.nrow();

  matxym = mmult(matx, maty);
  matxysum = sum_offd(matxym);
  xsum = sum_offd_multiply(matx);
  ysum = sum_offd_multiply(maty);
  sr = (matxysum - xsum - ysum)/(r*(r-1)*(r-2));

  matxye = mmult_ele(matx, maty);
  matxyesum = sum_offd(matxye);
  matxsum = sum_offd(matx);
  matysum = sum_offd(maty);

  sa_value = matxsum*matysum/(r*r*(r-1)*(r-1));
  in_value = matxyesum/(r*(r-1))+sa_value-2*sr;

  out[0] = sa_value;
  out[1] = in_value;

  return out;
}


NumericVector one_pcor(const NumericMatrix& a, const NumericMatrix& b, int rr)
{
  NumericMatrix matx(a.nrow()-1, a.ncol());
  NumericMatrix maty(b.nrow()-1, b.ncol());
  NumericVector out(2);

  matx = subtractr(a, rr);
  maty = subtractr(b, rr);

  out = insa(matx, maty);

  return out;
}


NumericVector invsa_va(const NumericMatrix& a, const NumericMatrix& b)
{
  double sa = 0, inv = 0;
  NumericVector out(2);
  NumericMatrix res(a.nrow(), 2);

  for (int i = 0; i < a.nrow(); i++) {
    res(i,_) = one_pcor(a,b,i);
  }

  for (int i = 0; i < a.nrow(); i++) {
    sa += res(i,0);
  }

  for (int i = 0; i <a.nrow(); i++) {
    inv += res(i,1);
  }

  out[0] = inv/a.nrow();
  out[1] = sa/a.nrow();

  return out;
}


// [[Rcpp::export]]
double pcov_va(const NumericMatrix& a, const NumericMatrix& b){
  double out;
  NumericVector pcs;

  pcs = invsa_va(a, b);
  out = std::sqrt(pcs[0]);

  return out;
}


// [[Rcpp::export]]
double pcor_va(const NumericMatrix& a, const NumericMatrix& b){
  double out, pca, pcb;

  pca = std::sqrt(pcov_va(a,a));
  pcb = std::sqrt(pcov_va(b,b));
  out = pcov_va(a,b)/(pca*pcb);

  return out;
}


// [[Rcpp::export]]
double t_va(const NumericMatrix& a, const NumericMatrix& b){
  double out;
  NumericVector pcs;

  pcs = invsa_va(a,b);
  out = pcs[0]/(3.14159265358979323846*3.14159265358979323846 -pcs[1]);

  return out;
}



