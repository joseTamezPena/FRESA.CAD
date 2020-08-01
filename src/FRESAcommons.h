/* FRESA.CAD: utilities for building and testing formula-based  
	models (linear, logistic or COX) for Computer Aided Diagnosis/Prognosis  
	applications.  Utilities include data adjustment, univariate analysis,  
	model building, model-validation, longitudinal analysis, reporting and visualization..  
 
   This program is free software under the terms of the  
   GPL Lesser General Public License as published by 
   the Free Software Foundation, either version 2 of the License, or 
   (at your option) any later version. 
   
   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of 
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   
    
   Jose Tamez and Israel Alanis 
   
*/ 
 
 
#ifndef Fresa_commons_h 
#define Fresa_commons_h 
 
 
#include <RcppArmadillo.h> 
#include <R.h> 
#include <Rmath.h> 
#ifdef _OPENMP 
#include <omp.h> 
#endif 
#include <R_ext/Utils.h> 
#include <exception> 
#include <iostream> 
#include <iomanip> 
#include <map> 
#include <cmath> 
 
// #define HAVE_UINTPTR_T 
#define CSTACK_DEFNS 7 
#define ARMA_NO_DEBUG 
#define ARMA_USE_BLAS 
 
using namespace arma; 
using namespace Rcpp; 
static const double THRESH = 36.0436534112975; 
static const double MTHRESH = -36.0436534112975; 
static const double INVEPS = 1.0/2.220446e-16; 
static const double DOUBLEEPS = 2.220446e-16; 
static const double MAXEXP = 4503599727262010.0; 
 
struct getVReclass{ 
  	 vec z_IDIs; 
	 vec z_NRIs; 
	 vec IDIs; 
	 vec NRIs; 
	 vec tz_IDIs; 
	 vec tz_NRIs; 
	 vec tIDIs; 
	 vec tNRIs; 
 
  }; 
   
 struct gvarNeRI{ 
	 vec tP_value; 
	 vec BinP_value; 
	 vec WilcoxP_value; 
	 vec FP_value; 
	 vec NeRIs; 
	 vec testData_tP_value; 
	 vec testData_BinP_value; 
	 vec testData_WilcoxP_value; 
	 vec testData_FP_value; 
	 vec testData_NeRIs; 
	}; 
	 
 struct improvedRes { 
	double p1; 
	double p2; 
	double NeRI; 
	double pvalue; 
	double binom_pValue; 
	double wilcox_pValue; 
	double t_test_pValue; 
	double F_test_pValue; 
  }; 
 
 
 
double qnorm(double p, double mu, double sigma); 
void chinv2(const mat &matrix , int n); 
int cholesky2(const mat &matrix, int n, double toler); 
void chsolve2(const mat &matrix, int n, vec &y); 
vec coxfit(int  maxiter,const  vec &time,const  vec &status,mat &covar,const vec &offset,const vec &weights,vec &strata,   int method, double eps, double toler,vec &beta,    int doscale); 
vec logit_link(const vec &mu); 
vec logit_linkinv(const vec &eta); 
vec logit_mu_eta(const vec &eta); 
vec binomial_dev_resids(const vec &y,const vec &mu,const vec &wt); 
vec modelFittingFunc(const mat &ymat,const mat &XP,const std::string & type,const vec &weights); 
vec modelFittingFunc(const mat &ymat,const mat &XP,const std::string & type); 
vec predictForFresaFunc(const vec &cf,const mat &newdata,const std::string & typ, const std::string & opc);		 
vec improveProbFuncSamples(const vec &x1,const vec &x2,const vec &y,unsigned int samples,double se_nri=0, double se_idi=0);  
vec improveProbFunc(const vec &x1,const vec &x2,const vec &y,double se_nri=0, double se_idi=0);  
getVReclass getVarBinFunc(const mat &dataframe,const std::string & type, const mat &independentFrame,const mat &bestFrame,const mat &bestTestFrame); 
double rocAUC(const vec &controls,const vec &cases, const std::string & direction); 
vec Fresarank(const vec &xi); 
vec residualForFRESAFunc(const vec &cf,const mat &newdata,const std::string & typ, const std::string & type,const mat &outcome); 
improvedRes improvedResidualsFunc(const vec &oldResiduals,const vec &newResiduals, const std::string & testType,unsigned int samples); 
improvedRes improvedResidualsFunc(const vec &oldResiduals,const vec &newResiduals, const std::string & testType); 
double pttest(const vec &xt, const vec &y,const std::string &tail); 
double ttest(const vec &x, const vec &y , double mu, bool paired, bool var_equal, const std::string & tail); 
double wilcoxtest(const vec &xt,const vec &y , double mu, bool paired, const std::string & tail,bool correct); 
double binomtest(double x, double n, double p , const std::string & tail); 
gvarNeRI getVarResFunc(const mat &dataframe, const std::string & type,const mat &testdataP,int testsamples,const mat &fullFrameP,const mat &fulltestFrameP); 
uvec equSamples(const mat &inputsample, unsigned int sort_indx=1,int histbins=5,double omin=0,double range=0); 
 
#endif 
// END 
