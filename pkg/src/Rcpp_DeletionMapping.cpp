#include "Rcpp_DeletionMapping.h"
#include "R_ext/Rdynload.h"
#include "Rcpp.h"

RcppExport SEXP RcppDeletionMapping(SEXP Arg1, SEXP Arg2, 
                          SEXP Arg3, SEXP Arg4, SEXP Arg5,
                          SEXP Arg6, SEXP Arg7, 
                          SEXP Arg8, SEXP Arg9,
                          SEXP Arg10, SEXP Arg11){
    using namespace Rcpp ;

    Environment delmap("package:delmap");
    Function iDeletionMapping = delmap["iDeletionMapping"]; 
   
    NumericMatrix dmap(Arg1);
    NumericMatrix imissing(Arg2);
    NumericVector psampled(Arg3);
    NumericMatrix odmap(Arg4);
    NumericVector temp(Arg5);
    CharacterVector  method(Arg6);
    NumericVector otlength(Arg7);
    NumericVector besttlength(Arg8);
    NumericMatrix bestmap(Arg9);
    int outer = as<int>(Arg10);
    int inner = as<int>(Arg11);

    Rprintf( " %i ", inner);

    List newz = List::create(Named("bestmap") = bestmap,
                             Named("besttlength") = besttlength,
                             Named("odmap") = odmap,
                             Named("imissing") = imissing,
                             Named("otlength") = otlength);
    for(int ii=0; ii <outer; ii++)
    {
      for( int jj=0; jj < inner; jj++)
     {
        List z = iDeletionMapping(dmap, 
                               newz("imissing"), 
                              psampled,
                              newz("odmap"), 
                              temp, 
                              method, 
                              newz("otlength"), 
                              newz("besttlength"), 
                              newz("bestmap"));
    
        newz("bestmap") = z("bestmap");
        newz("besttlength") = z("besttlength");
        newz("odmap") = z("odmap");
        newz("imissing") = z("imissing");
        newz("otlength") = z("otlength");
     }
   }

    return wrap(newz) ;
}



