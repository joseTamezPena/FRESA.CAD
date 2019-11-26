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
 
 
#include <R.h> 
#include <Rinternals.h> 
#include <stdlib.h> // for NULL 
#include <R_ext/Rdynload.h> 
 
 
/* .Call calls */ 
 
 
extern "C" SEXP bootstrapValidationBinCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
extern "C" SEXP bootstrapValidationResCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
extern "C" SEXP equalizedSampling(SEXP, SEXP, SEXP); 
extern "C" SEXP ForwardResidualModelCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
extern "C" SEXP improvedResidualsCpp(SEXP, SEXP, SEXP, SEXP); 
extern "C" SEXP improveProbCpp(SEXP, SEXP, SEXP); 
extern "C" SEXP modelFittingCpp(SEXP, SEXP, SEXP); 
extern "C" SEXP wtmodelFittingCpp(SEXP, SEXP, SEXP, SEXP); 
extern "C" SEXP predictForFresa(SEXP, SEXP, SEXP, SEXP); 
extern "C" SEXP rankInverseNormalCpp(SEXP, SEXP, SEXP, SEXP, SEXP); 
extern "C" SEXP ReclassificationFRESAModelCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
extern "C" SEXP logRank(SEXP, SEXP, SEXP); 
extern "C" SEXP SLRNullDistribution(SEXP ,SEXP ,SEXP ,SEXP ,SEXP); 
extern "C" SEXP SLRDistribution(SEXP ,SEXP ,SEXP ,SEXP ,SEXP); 
 
static const R_CallMethodDef CallEntries[] = { 
    {"bootstrapValidationBinCpp",     (DL_FUNC) &bootstrapValidationBinCpp,      6}, 
    {"bootstrapValidationResCpp",     (DL_FUNC) &bootstrapValidationResCpp,      6}, 
    {"equalizedSampling",             (DL_FUNC) &equalizedSampling,              3}, 
    {"ForwardResidualModelCpp",       (DL_FUNC) &ForwardResidualModelCpp,       15}, 
    {"improvedResidualsCpp",          (DL_FUNC) &improvedResidualsCpp,           4}, 
    {"improveProbCpp",                (DL_FUNC) &improveProbCpp,                 3}, 
    {"modelFittingCpp",               (DL_FUNC) &modelFittingCpp,                3}, 
    {"wtmodelFittingCpp",             (DL_FUNC) &wtmodelFittingCpp,              3}, 
    {"predictForFresa",               (DL_FUNC) &predictForFresa,                4}, 
    {"rankInverseNormalCpp",          (DL_FUNC) &rankInverseNormalCpp,           5}, 
    {"ReclassificationFRESAModelCpp", (DL_FUNC) &ReclassificationFRESAModelCpp, 15}, 
    {"logRank",                       (DL_FUNC) &logRank,                        3}, 
    {"SLRNullDistribution",           (DL_FUNC) &SLRNullDistribution,            5}, 
    {"SLRDistribution",               (DL_FUNC) &SLRDistribution,                5}, 
    {NULL, NULL, 0} 
}; 
 
extern "C" void R_init_FRESA_CAD(DllInfo *dll) 
{ 
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL); 
    R_useDynamicSymbols(dll, FALSE); 
} 
 
extern "C" void R_init_FRESACAD(DllInfo *dll) 
{ 
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL); 
    R_useDynamicSymbols(dll, FALSE); 
} 
 
