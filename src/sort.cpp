#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <stdio.h>
using namespace Rcpp;



//' @title Compute how many days in month M, year Y.
//' @description Compute how many days in month M, year Y, as a subfunction for water supply prediction.
//' @param Y Year
//' @param M Month
//' @return Days of Month M in Year Y.
//' @examples
//' \dontrun{
//' dates(2022,12)
//' }
//' @export
// [[Rcpp::export]]
int dates(int Y,int M){
  int a[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
  if (Y % 4 == 0){
    if (Y % 100 == 0)
    {
      if (Y % 400 == 0)
        a[2]++;
    }
    else
      a[2]++;
  }
  return a[M];
}