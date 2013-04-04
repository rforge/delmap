#include "Rcpp_DeletionMapping.h"
#include "R_ext/Rdynload.h"
#include "Rcpp.h"

RcppExport SEXP RcppDeletionMapping(SEXP Arg1,  
                          SEXP Arg3, SEXP Arg4, SEXP Arg5,
                          SEXP Arg6, SEXP Arg7, 
                          SEXP Arg8, SEXP Arg9,
                          SEXP Arg10, SEXP Arg11){
    using namespace Rcpp ;

    Environment delmap("package:delmap");
    Function iDeletionMapping = delmap["iDeletionMapping"]; 

    Environment base("package:base");
    Function print = base["print.default"];
 
    NumericMatrix dmap(Arg1);
    NumericVector psampled(Arg3);
    NumericMatrix odmap(Arg4);
    double cooling =  as<double>(Arg5);
    CharacterVector  method(Arg6);
    NumericVector otlength(Arg7);
    NumericVector besttlength(Arg8);
    NumericMatrix bestmap(Arg9);
    int outer = as<int>(Arg10);
    int inner = as<int>(Arg11);

    double temp;
    int counter=0;

    
    NumericVector  t(outer*inner);



    temp = 0.5 / pow(cooling, outer);

    List newz = List::create(Named("dmap") = dmap,
                             Named("bestmap") = bestmap,
                             Named("besttlength") = besttlength,
                             Named("odmap") = odmap,
                             Named("otlength") = otlength,
                             Named("t") = t );
    for(int ii=0; ii <outer; ii++)
    {
      for( int jj=0; jj < inner; jj++)
     {
        List z = iDeletionMapping(newz("dmap"), 
                              psampled,
                              newz("odmap"), 
                              temp, 
                              method, 
                              newz("otlength"), 
                              newz("besttlength"), 
                              newz("bestmap"));
   
 
        newz("dmap") = z("dmap");
        newz("bestmap") = z("bestmap");
        newz("besttlength") = z("besttlength");
        newz("odmap") = z("odmap");
        newz("otlength") = z("otlength");

          t(counter) = as<int>(z("besttlength")) ;
        counter++;
        
     }
      temp = temp * cooling;
   }
    return wrap(newz) ;
}



