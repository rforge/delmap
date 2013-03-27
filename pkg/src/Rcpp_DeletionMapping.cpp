#include "Rcpp_DeletionMapping.h"
#include "R_ext/Rdynload.h"
#include "Rcpp.h"

RcppExport SEXP Rcpp_DeletionMapping(){
    using namespace Rcpp ;

    Environment Envdelmap("package:delmap");
    Function testf = Envdelmap["testf"]; 
    
    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    NumericVector z ;
    z = testf();

    return z ;
}



