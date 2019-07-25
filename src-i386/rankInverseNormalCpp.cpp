/* RankInverseNormalCpp: utility is part of FRESA.CAD

   This program is free software under the terms of the 
   GPL Lesser General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
   
   Jose Tamez and Israel Alanis
  
*/
#include <Rcpp.h>

using namespace Rcpp;

extern "C" SEXP rankInverseNormalCpp(SEXP _nrows,SEXP _obs, SEXP _minvalue, SEXP _maxvalue, SEXP _dataframe)
{  //R_CStackLimit=(uintptr_t)-1;
	
	double zMax=8.0;
	double zout=0;
	double ranmin = R::pnorm(-zMax, 0.0, 1.0, 1, 0);
	double ranmax = R::pnorm(zMax, 0.0, 1.0, 1, 0);
	unsigned int nrows= Rcpp::as<unsigned int>(_nrows);
	NumericVector obs(_obs);
	double min_value= Rcpp::as<double>(_minvalue);
	double max_value= Rcpp::as<double>(_maxvalue);
	NumericVector dataframe(_dataframe);
	int size=dataframe.size();
	NumericVector out(nrows);
	double range = (max_value-min_value);
	double zmin = R::qnorm(1.0/(1.0+size), 0.0, 1.0, 1, 0);
	double zmx = R::qnorm(static_cast<double>(size)/(1.0+size), 0.0, 1.0, 1, 0);

//	int i=1;
        for (unsigned int j = 0; j < nrows; ++j)
        {
				double mrank=0;
		        int i=1;
				zout=0.0;
		        while ((obs[j] >= dataframe[i-1]) & (i <= size))
		        {
		          i=i+1; 
		        }
		        i = i-1;
		        if (i==0) 
		        {
				  zout = zmin*(1.0+2.0*(dataframe[0]-obs[j])/range);
		        }
		        else
				{
					if (i == size)
					{
						zout = zmx*(1.0+2.0*(obs[j]-dataframe[size-1])/range);
					}
					else
					{
						mrank = (i+(obs[j]-dataframe[i-1])/(dataframe[i]-dataframe[i-1]))/(1.0+size);
					}
				}
				if (zout == 0.0)
				{
					if (mrank<ranmin) mrank=R::pnorm(-zMax, 0.0, 1.0, 1, 0);
					if (mrank>ranmax) mrank=R::pnorm(zMax, 0.0, 1.0, 1, 0);
					zout = R::qnorm(mrank, 0.0, 1.0, 1, 0);
				}
				out[j] = zout;
		}
		 	
		return out;

}
