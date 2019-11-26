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
 
 
#include "FRESAcommons.h" 
 
#define EPS 1e-4		 
 
#define RANDIN  GetRNGstate() 
#define RANDOUT PutRNGstate() 
#define UNIF unif_rand() 
 
 
#define MAX_TIES 1000 
 
 
 
 
extern "C" SEXP modelFittingCpp(SEXP _ymat,SEXP _xmat,SEXP _type) 
{ 
	std::string type = Rcpp::as<std::string>(_type); 
	Rcpp::NumericMatrix rymat(_ymat); 
	Rcpp::NumericMatrix rX(_xmat); 
	mat ymat(rymat.begin(), rymat.rows(), rymat.cols(),false); 
	mat X(rX.begin(), rX.rows(), rX.cols(),false); 
	vec betas=modelFittingFunc(ymat,X,type); 
	vec linearPredictors = predictForFresaFunc(betas,X,"linear",type); 
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("coefficients")=Rcpp::wrap(betas.t()), 
										   Rcpp::Named("linear.predictors")=Rcpp::wrap(linearPredictors));  
	return result; 
} 
 
extern "C" SEXP wtmodelFittingCpp(SEXP _ymat,SEXP _xmat,SEXP _type,SEXP _wts) 
{ 
	std::string type = Rcpp::as<std::string>(_type); 
	Rcpp::NumericMatrix rymat(_ymat); 
	Rcpp::NumericMatrix rX(_xmat); 
	Rcpp::NumericVector wts(_wts); 
	mat ymat(rymat.begin(), rymat.rows(), rymat.cols(),false); 
	mat X(rX.begin(), rX.rows(), rX.cols(),false); 
	vec betas=modelFittingFunc(ymat,X,type,wts); 
	vec linearPredictors = predictForFresaFunc(betas,X,"linear",type); 
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("coefficients")=Rcpp::wrap(betas.t()), 
										   Rcpp::Named("linear.predictors")=Rcpp::wrap(linearPredictors));  
	return result; 
} 
 
extern "C" SEXP improvedResidualsCpp(SEXP _oldResiduals,SEXP _newResiduals,SEXP _testType,SEXP _samples) 
{ 
	std::string testType = Rcpp::as<std::string>(_testType); 
	Rcpp::NumericVector roldResiduals(_oldResiduals); 
	Rcpp::NumericVector rnewResiduals(_newResiduals); 
	unsigned int samples = Rcpp::as<unsigned int>(_samples); 
	vec oldResiduals(roldResiduals.begin(), roldResiduals.size(),false); 
	vec newResiduals(rnewResiduals.begin(), rnewResiduals.size(),false); 
	improvedRes impred = improvedResidualsFunc(oldResiduals,newResiduals,testType,samples); 
 
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("p1")=Rcpp::wrap(impred.p1), 
										   Rcpp::Named("p2")=Rcpp::wrap(impred.p2), 
										   Rcpp::Named("NeRI")=Rcpp::wrap(impred.NeRI), 
										   Rcpp::Named("p.value")=Rcpp::wrap(impred.pvalue), 
										   Rcpp::Named("BinP.value")=Rcpp::wrap(impred.binom_pValue), 
										   Rcpp::Named("WilcoxP.value")=Rcpp::wrap(impred.wilcox_pValue), 
										   Rcpp::Named("tP.value")=Rcpp::wrap(impred.t_test_pValue), 
										   Rcpp::Named("FP.value")=Rcpp::wrap(impred.F_test_pValue) 
										   ); 
	return result; 
} 
 
 
 extern "C" SEXP improveProbCpp(SEXP _x1,SEXP _x2,SEXP _y) 
{ 
	Rcpp::NumericVector rx1(_x1); 
	Rcpp::NumericVector rx2(_x2); 
	Rcpp::NumericVector ry(_y); 
 
	vec x1(rx1.begin(), rx1.size(),false); 
	vec x2(rx2.begin(), rx2.size(),false); 
	vec y(ry.begin(), ry.size(),false); 
	vec imp=improveProbFunc(x1,x2,y); 
	 
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("z.idi")=Rcpp::wrap(imp[0]), 
										   Rcpp::Named("z.nri")=Rcpp::wrap(imp[1]), 
										   Rcpp::Named("idi")=Rcpp::wrap(imp[2]), 
										   Rcpp::Named("nri")=Rcpp::wrap(imp[3]) 
										   ); 
	return result; 
} 
 
 
extern "C" SEXP predictForFresa(SEXP _cf,SEXP _newdata,SEXP _typ,SEXP _opc) 
{ 
	std::string typ = Rcpp::as<std::string>(_typ); 
	std::string opc = Rcpp::as<std::string>(_opc); 
	Rcpp::NumericVector rcf(_cf); 
	Rcpp::NumericMatrix rnewdata(_newdata); 
	vec cf(rcf.begin(), rcf.size(),false); 
	mat newdata(rnewdata.begin(), rnewdata.rows(), rnewdata.cols(), false); 
	vec prediction = predictForFresaFunc(cf,newdata,typ,opc); 
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("prediction")=Rcpp::wrap(prediction)); 
	return result; 
} 
 
extern "C" SEXP equalizedSampling(SEXP _thedata,SEXP _indx,SEXP _breakWidth) 
{ 
	Rcpp::NumericMatrix rnewdata(_thedata); 
	unsigned int indx = Rcpp::as<unsigned int>(_indx); 
	unsigned int breakWidth = Rcpp::as<unsigned int>(_breakWidth); 
	mat newdata(rnewdata.begin(), rnewdata.rows(), rnewdata.cols(), false); 
	uvec uidex = equSamples(newdata,indx,breakWidth); 
	mat equmatrix = newdata.rows(uidex); 
	Rcpp::NumericMatrix result = Rcpp::wrap(equmatrix); 
	return result; 
} 
 
extern "C" SEXP logRank(SEXP _times,SEXP _status,SEXP _groups) 
{ 
	Rcpp::NumericVector rtimes(_times); 
	Rcpp::NumericVector rstatus(_status); 
	Rcpp::NumericVector rgroups(_groups); 
	vec times(rtimes.begin(), rtimes.size(),false); 
	vec status(rstatus.begin(), rstatus.size(),false); 
	vec groups(rgroups.begin(), rgroups.size(),false); 
 
	 
	uvec uidex = sort_index(times); 
	vec ntimes = times(uidex); 
	vec ngroups = groups(uidex); 
	vec nstatus = status(uidex); 
	int nt=times.size(); 
	int nt1 = nt-1; 
	int eqt=0; 
	int repeats=1; 
	for (int i=0; i<nt1; ++i) 
	{ 
		eqt += (ntimes[i]==ntimes[i+1])&&(ngroups[i]!=ngroups[i+1]); 
	} 
	if (eqt > 0) repeats = 10; //repeat if there is subjects with the same event time 
	 
	double LR=0; 
	double VR=0; 
	int n1=sum(groups); 
	for (int n=0;n<repeats;n++) 
	{ 
		if (repeats>1) 
		{ 
			ntimes = times+10*DOUBLEEPS*randu<vec>(nt); 
			uidex = sort_index(ntimes); 
			ngroups = groups(uidex); 
			nstatus = status(uidex); 
		} 
		double atRisk=nt; 
		double atRiskO1=n1; 
		for (int i=0; i<nt1; ++i) 
		{ 
			double expected=nstatus[i]*atRiskO1/atRisk; 
			LR += nstatus[i]*ngroups[i]-expected; 
			VR += expected*(atRisk-atRiskO1)*(atRisk-nstatus[i])/(atRisk*(atRisk-1.0)); 
			atRiskO1 -= ngroups[i]; 
			--atRisk; 
		}	 
		LR += nstatus[nt1]*ngroups[nt1]-nstatus[nt1]*atRiskO1/atRisk; 
	} 
	LR = LR/repeats; 
	VR = VR/repeats; 
	double SLR = LR/sqrt(VR); 
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("LR")=Rcpp::wrap(LR), 
										   Rcpp::Named("VR")=Rcpp::wrap(VR), 
										   Rcpp::Named("SLR")=Rcpp::wrap(SLR) 
										   ); 
	return result; 
} 
 
 
double standarizedLogRank(const vec &itimes,const vec &istatus,const vec &igroups,unsigned int type=0) 
{ 
	uvec uidex = sort_index(itimes); 
	vec groups = igroups(uidex); 
	vec status = istatus(uidex); 
	double LR=0; 
	double VR=0; 
	int n1=sum(groups); 
	int nt=itimes.size(); 
	double atRisk=nt; 
	double atRiskO1=n1; 
	int nt1 = nt-1; 
	int i; 
	for (i=0; i<nt1; ++i) 
	{ 
		double expected=status[i]*atRiskO1/atRisk; 
		LR += status[i]*groups[i]-expected; 
		VR += expected*(atRisk-atRiskO1)*(atRisk-status[i])/(atRisk*(atRisk-1.0)); 
		atRiskO1 -= groups[i]; 
		--atRisk; 
	} 
	LR += status[i]*groups[i]-status[i]*atRiskO1/atRisk; 
	if (VR>0) 
	{ 
		if (type==1) LR=LR/sqrt(VR); 
		if (type==2) LR=LR*LR/VR; 
	} 
	return LR; 
} 
 
extern "C" SEXP SLRNullDistribution(SEXP _times,SEXP _status,SEXP _groups,SEXP _samples,SEXP _type) 
{ 
	Rcpp::NumericVector rtimes(_times); 
	Rcpp::NumericVector rstatus(_status); 
	Rcpp::NumericVector rgroups(_groups); 
	unsigned int samples = Rcpp::as<unsigned int>(_samples); 
	unsigned int type = Rcpp::as<unsigned int>(_type); 
	const vec times(rtimes.begin(), rtimes.size(),false); 
	const vec status(rstatus.begin(), rstatus.size(),false); 
	const vec groups(rgroups.begin(), rgroups.size(),false); 
	vec logranksamples(samples); 
	unsigned int nt = times.size(); 
 
 
	#pragma omp parallel for schedule(dynamic)  
	for (unsigned int i = 0; i<samples; i++) 
	{ 
		logranksamples[i] = standarizedLogRank(randu<vec>(nt),shuffle(status),groups,type); 
    } 
	Rcpp::NumericVector result = Rcpp::wrap(logranksamples); 
	return result; 
} 
 
/***** Bootstapped distribution of the observed LR 
 
*/ 
extern "C" SEXP SLRDistribution(SEXP _times,SEXP _status,SEXP _groups,SEXP _samples,SEXP _type) 
{ 
	Rcpp::NumericVector rtimes(_times); 
	Rcpp::NumericVector rstatus(_status); 
	Rcpp::NumericVector rgroups(_groups); 
	unsigned int samples = Rcpp::as<unsigned int>(_samples); 
	unsigned int type = Rcpp::as<unsigned int>(_type); 
	const vec times(rtimes.begin(), rtimes.size(),false); 
	const vec status(rstatus.begin(), rstatus.size(),false); 
	const vec groups(rgroups.begin(), rgroups.size(),false); 
	vec logranksamples(samples); 
	uvec g1 = find(groups==1); 
	uvec g0 = find(groups==0); 
	const vec timesg1 = times(g1); 
	const vec timesg0 = times(g0); 
	const vec statusg1 = status(g1); 
	const vec statusg0 = status(g0); 
	const vec groupsg1 = groups(g1); 
	const vec groupsg0 = groups(g0); 
 
	unsigned int nt = times.size(); 
	unsigned int ng0 = timesg0.size(); 
	unsigned int ng0t1 = ng0-1; 
	unsigned int ng1 = timesg1.size(); 
	unsigned int ng1t1 = ng1-1; 
	 
	vec ngroups(nt); 
	ngroups.head(ng0) =  groupsg0; 
	ngroups.tail(ng1) =  groupsg1; 
 
	#pragma omp parallel for schedule(dynamic)  
	for (unsigned int i = 0; i<samples; i++) 
	{ 
		vec ntimes(nt); 
		vec nstatus(nt); 
		uvec bootsamplesg0 = randi<uvec>(ng0, distr_param(0,ng0t1)); 
		uvec bootsamplesg1 = randi<uvec>(ng1, distr_param(0,ng1t1)); 
		ntimes.head(ng0) = timesg0(bootsamplesg0)+10.0*DOUBLEEPS*randu<vec>(ng0); 
		ntimes.tail(ng1) = timesg1(bootsamplesg1)+10.0*DOUBLEEPS*randu<vec>(ng1); 
		nstatus.head(ng0) = statusg0(bootsamplesg0); 
		nstatus.tail(ng1) = statusg1(bootsamplesg1); 
		logranksamples[i] = standarizedLogRank(ntimes,nstatus,ngroups,type); 
    } 
	Rcpp::NumericVector result = Rcpp::wrap(logranksamples); 
	return result; 
} 
 
                                                           
double qnorm(double p, double mu, double sigma) 
{ 
	return R::qnorm(p, mu, sigma, 1, 0); 
} 
 
 
/* 
** This code is based on the chinv2 function from survival package  
** Reference: Therneau T (2014). A Package for Survival Analysis in S. R package version 2.37-7, http://CRAN.R-project.org/package=survival. 
*/ 
void chinv2(mat &matrix , int n) 
{ 
	double temp; 
	int i,j,k; 
	for (i=0; i<n; i++) 
	{ 
		if (matrix(i,i) >0)  
		{ 
		  	matrix(i,i) = 1/matrix(i,i);   
		 	for (j= (i+1); j<n; j++)  
		 	{ 
		   		matrix(i,j) = -matrix(i,j); 
		   		for (k=0; k<i; k++)    
					matrix(k,j) += matrix(i,j)*matrix(k,i); 
		   	} 
		} 
	} 
	for (i=0; i<n; i++)  
	{ 
		if (matrix(i,i)==0)  
		{  
			for (j=0; j<i; j++) matrix(i,j)=0; 
			for (j=i; j<n; j++) matrix(j,i)=0; 
		} 
		else  
		{ 
			for (j=(i+1); j<n; j++)  
			{ 
				temp = matrix(i,j)*matrix(j,j); 
				if (j!=i) matrix(j,i) = temp; 
				for (k=i; k<j; k++) 
				matrix(k,i) += temp*matrix(k,j); 
			} 
		} 
	} 
} 
 
/* 
** This code is based on the cholesky2 function from survival package  
** Reference: Therneau T (2014). A Package for Survival Analysis in S. R package version 2.37-7, http://CRAN.R-project.org/package=survival. 
*/ 
 
int cholesky2(mat &matrix, int n, double toler) 
{ 
    double temp; 
    double eps=0; 
    int  i,j,k; 
    for (i=0; i<n; i++)  
    { 
		if (matrix(i,i)> eps)  eps = matrix(i,i); 
		for (j=(i+1); j<n; j++)  matrix(i,j) = matrix(j,i); 
	} 
    eps *= toler; 
 
    int rank=0; 
    int nonneg=1; 
    for (i=0; i<n; i++)  
    { 
		double pivot = matrix(i,i); 
		if (pivot < eps)  
		{ 
		    matrix(i,i) =0; 
		    if (pivot < -8*eps) nonneg= -1; 
	    } 
		else   
		{ 
		    rank++; 
		    for (j=(i+1); j<n; j++)  
		    { 
				temp = matrix(i,j)/pivot; 
				matrix(i,j) = temp; 
				matrix(j,j) -= temp*temp*pivot; 
				for (k=(j+1); k<n; k++) matrix(j,k) -= temp*matrix(i,k); 
			} 
	    } 
	} 
    return(rank * nonneg); 
} 
 
/* 
** This code is based on the coxfit6 function from survival package  
** Reference: Therneau T (2014). A Package for Survival Analysis in S. R package version 2.37-7, http://CRAN.R-project.org/package=survival. 
** calls functions:  cholesky2, chsolve2, chinv2. 
*/ 
void chsolve2(const mat &matrix, int n, vec &y) 
{ 
	int i,j; 
	double temp; 
 
	for (i=0; i<n; i++)  
	{ 
		temp = y[i]; 
		for (j=0; j<i; j++) 
		   temp -= y[j] * matrix(j,i); 
		y(i) = temp ; 
	} 
	for (i=(n-1); i>=0; i--)  
	{ 
		if (matrix(i,i)==0)  y[i] =0; 
		else  
		{ 
			temp = y[i]/matrix(i,i); 
			for (j= i+1; j<n; j++) 
				temp -= y[j]*matrix(i,j); 
			y[i] = temp; 
	  	} 
	} 
} 
/* 
** This code is based on the coxfit6 function from survival package  
** Reference: Therneau T (2014). A Package for Survival Analysis in S. R package version 2.37-7, http://CRAN.R-project.org/package=survival. 
** calls functions:  cholesky2, chsolve2, chinv2. 
*/ 
vec coxfit(int  maxiter,const  vec &time,const  vec &status, mat &covar,const vec &offset,const vec &weights,vec &strata,   int method, double eps, double toler,vec &beta,    int doscale)  
{ 
    int i,j,k, obs;   
    double  wtave; 
    double  denom=0, zbeta, risk; 
    double  temp, temp2; 
    int     ndead;   
    double  tdeath=0;   
    double  newlk=0; 
    double  dtime, d2; 
    double  deadwt;   
    double  efronwt; 
    int     halving;    
    int     nrisk;    
    int     nobs=offset.n_elem,  nvar  = covar.n_cols; 
    vec a(nvar); 
    vec a2(nvar); 
    vec maxbeta(nvar); 
    vec newbeta(nvar); 
    vec scale(nvar); 
    mat cmat(nvar,nvar); 
    mat cmat2(nvar,nvar); 
    mat imat(nvar,nvar); 
    vec means(nvar); 
    vec u(nvar); 
    vec loglik(2); 
//    double sctest; 
//   int flag; 
    int iter; 
    tdeath=0; temp2=0; 
    for (i=0; i<nobs; i++)  
    { 
		temp2 += weights(i); 
		tdeath += weights(i) * status(i); 
    }	 
    for (i=0; i<nvar; i++)  
    { 
		temp=0; 
		for (obs=0; obs<nobs; obs++)  
	    	temp += weights(obs) * covar(obs,i); 
		temp /= temp2; 
		means(i) = temp; 
		for (obs=0; obs<nobs; obs++) covar(obs,i) -=temp; 
		if (doscale==1)  
		{  
	    	temp =0; 
	    	for (obs=0; obs<nobs; obs++)  
	    	{     
				temp += weights(obs) * std::abs(covar(obs,i)); 
	    	} 
	    	if (temp > 0) temp = temp2/temp;    
	    	else temp=1.0;  
	    	scale(i) = temp; 
	    	for (obs=0; obs<nobs; obs++)  covar(obs,i) *= temp; 
	    } 
	} 
    if (doscale==1)  
    { 
		for (i=0; i<nvar; i++) beta(i) /= scale(i);  
	} 
    else  
    { 
		for (i=0; i<nvar; i++) scale(i) = 1.0; 
	} 
    strata(nobs-1) =1; 
    loglik(1) =0; 
    for (i=0; i<nvar; i++) 
    { 
		u(i) =0; 
		a2(i) =0; 
		for (j=0; j<nvar; j++)  
		{ 
		    imat(i,j) =0 ; 
		    cmat2(i,j) =0; 
	    } 
	} 
    for (obs=nobs-1; obs>=0; )  
    { 
		if (strata(obs) == 1)  
		{ 
		    nrisk =0 ;   
		    denom = 0; 
		    for (i=0; i<nvar; i++)  
		    { 
				a(i) = 0; 
				for (j=0; j<nvar; j++) cmat(j,i) = 0; 
			} 
	    } 
		dtime = time(obs); 
		ndead =0;  
		deadwt =0;   
		efronwt=0;   
		while(obs >=0 &&time(obs)==dtime)  
		{ 
		    nrisk++; 
		    zbeta = offset(obs);     
		    for (i=0; i<nvar; i++) 
				zbeta += beta(i)*covar(obs,i); 
		    risk = exp(zbeta) * weights(obs); 
		    denom += risk; 
		    for (i=0; i<nvar; i++)  
		    { 
				a(i) += risk*covar(obs,i); 
				for (j=0; j<=i; j++) 
			    	cmat(j,i) += risk*covar(obs,i)*covar(obs,j); 
		    } 
 
		    if (status(obs)==1.0)  
		    { 
				ndead++; 
				deadwt += weights(obs); 
				efronwt += risk; 
				loglik(1)+= weights(obs)*zbeta; 
 
				for (i=0; i<nvar; i++)  
				    u(i) += weights(obs)*covar(obs,i); 
				if (method==1)  
				{  
				    for (i=0; i<nvar; i++)  
				    { 
						a2(i) +=  risk*covar(obs,i); 
						for (j=0; j<=i; j++) 
						    cmat2(j,i) += risk*covar(obs,i)*covar(obs,j); 
			        } 
				} 
		    } 
		    obs--; 
		    if (obs>=0) 
		    if (strata(obs)==1) break;   
 
	    } 
 
		if (ndead >0)  
		{  
		    if (method==0)  
		    {  
				loglik(1) -= deadwt* log(denom); 
		    
				for (i=0; i<nvar; i++)  
				{ 
				    temp2= a(i)/ denom;  
				    u(i) -=  deadwt* temp2; 
				    for (j=0; j<=i; j++) 
						imat(i,j) += deadwt*(cmat(j,i) - temp2*a(j))/denom; 
			    } 
			} 
		    else  
		    {  
				for (k=0; k<ndead; k++)  
				{ 
				    temp = static_cast<double>(k)/ ndead; 
				    wtave = deadwt/ndead; 
				    d2 = denom - temp*efronwt; 
				    loglik(1) -= wtave* log(d2); 
				    for (i=0; i<nvar; i++)  
				    { 
						temp2 = (a(i) - temp*a2(i))/ d2; 
						u(i) -= wtave *temp2; 
						for (j=0; j<=i; j++) 
						    imat(i,j) +=  (wtave/d2) * ((cmat(j,i) - temp*cmat2(j,i)) - temp2*(a(j)-temp*a2(j))); 
				    } 
			    } 
		 
				for (i=0; i<nvar; i++)  
				{ 
				    a2(i)=0; 
				    for (j=0; j<nvar; j++) cmat2(i,j)=0; 
		    	} 
			} 
	    } 
	}   
    loglik(0) = loglik(1); 
    for (i=0; i<nvar; i++)  
		maxbeta(i) = 20* std::sqrt(imat(i,i)/tdeath); 
    for (i=0; i<nvar; i++)  
		a(i) = u(i); 
 
//    flag= cholesky2(imat, nvar, toler); 
    cholesky2(imat, nvar, toler); 
    chsolve2(imat,nvar,a);        
    temp=0; 
    for (i=0; i<nvar; i++) temp +=  u(i)*a(i); 
//    sctest = temp;   
    for (i=0; i<nvar; i++) 
    { 
		newbeta(i) = beta(i) + a(i); 
	} 
    if (maxiter==0)  
    { 
		chinv2(imat,nvar); 
		for (i=0; i<nvar; i++)  
		{ 
		    beta(i) *= scale(i);  
		    u(i) /= scale(i); 
		    imat(i,i) *= scale(i)*scale(i); 
		    for (j=0; j<i; j++)  
		    { 
				imat(i,j) *= scale(i)*scale(j); 
				imat(j,i) = imat(i,j); 
			} 
    	} 
		goto finish; 
	} 
    halving =0 ;             
    for (iter=1; iter<= maxiter; iter++)  
    { 
		newlk =0; 
		for (i=0; i<nvar; i++)  
		{ 
		    u(i) =0; 
		    for (j=0; j<nvar; j++) 
			imat(j,i) =0; 
	    } 
		for (obs=nobs-1; obs>=0; )  
		{ 
		    if (strata(obs) == 1)  
		    {  
				denom = 0; 
				nrisk =0; 
				for (i=0; i<nvar; i++)  
				{ 
				    a(i) = 0; 
				    for (j=0; j<nvar; j++) cmat(j,i) = 0; 
			    } 
			} 
		    dtime = time(obs); 
		    deadwt =0; 
		    ndead =0; 
		    efronwt =0; 
		    while(obs>=0 && time(obs)==dtime)  
		    { 
				nrisk++; 
				zbeta = offset(obs); 
				for (i=0; i<nvar; i++) 
				    zbeta += newbeta(i)*covar(obs,i); 
				risk = exp(zbeta) * weights(obs); 
				denom += risk; 
 
				for (i=0; i<nvar; i++)  
				{ 
				    a(i) += risk*covar(obs,i); 
				    for (j=0; j<=i; j++) 
				    cmat(j,i) += risk*covar(obs,i)*covar(obs,j); 
			    } 
				if (status(obs)==1)  
				{ 
				    ndead++; 
				    deadwt += weights(obs); 
				    newlk += weights(obs) *zbeta; 
				    for (i=0; i<nvar; i++)  
						u(i) += weights(obs) *covar(obs,i); 
				    if (method==1)  
				    {  
						efronwt += risk; 
						for (i=0; i<nvar; i++)  
						{ 
						    a2(i) +=  risk*covar(obs,i); 
						    for (j=0; j<=i; j++) 
								cmat2(j,i) += risk*covar(obs,i)*covar(obs,j); 
						}    
			        } 
		  	    } 
				obs--; 
				if (obs>=0) 
				if (strata(obs)==1) break;  
		    } 
		    if (ndead >0)  
		    {   
				if (method==0)  
				{  
				    newlk -= deadwt* log(denom); 
				    for (i=0; i<nvar; i++)  
				    { 
						temp2= a(i)/ denom;   
						u(i) -= deadwt* temp2; 
						for (j=0; j<=i; j++) 
						  	imat(i,j) +=  (deadwt/denom)* (cmat(j,i) - temp2*a(j)); 
			        } 
			    } 
				else   
				{  
				    for (k=0; k<ndead; k++)  
				    { 
						temp = static_cast<double>(k) / ndead; 
						wtave= deadwt/ ndead; 
						d2= denom - temp* efronwt; 
						newlk -= wtave* log(d2); 
						for (i=0; i<nvar; i++)  
						{ 
						    temp2 = (a(i) - temp*a2(i))/ d2; 
						    u(i) -= wtave*temp2; 
						    for (j=0; j<=i; j++) 
								imat(i,j) +=  (wtave/d2)*((cmat(j,i) - temp*cmat2(j,i)) - temp2*(a(j)-temp*a2(j))); 
			            } 
			        } 
 
				    for (i=0; i<nvar; i++)  
				    {  
						a2(i) =0; 
						for (j=0; j<nvar; j++) cmat2(j,i) =0; 
			        } 
		        } 
			} 
	    }    
//		flag = cholesky2(imat, nvar, toler); 
		cholesky2(imat, nvar, toler); 
		if (std::abs(1-(loglik(1)/newlk))<= eps && halving==0)  
		{  
		    loglik(1) = newlk; 
		    chinv2(imat, nvar);      
		    for (i=0; i<nvar; i++) 
		    { 
				beta(i) = newbeta(i)*scale(i); 
				u(i) /= scale(i); 
				imat(i,i) *= scale(i)*scale(i); 
				for (j=0; j<i; j++)  
				{ 
				    imat(i,j) *= scale(i)*scale(j); 
				    imat(j,i) = imat(i,j); 
			    } 
		    } 
		    goto finish; 
		} 
		if (iter== maxiter) break;   
		if (newlk < loglik(1))    
		{    
			halving =1; 
			for (i=0; i<nvar; i++) 
			    newbeta(i) = (newbeta(i) + beta(i)) /2; 
		} 
		else  
		{ 
		    halving=0; 
		    loglik(1) = newlk; 
		    chsolve2(imat,nvar,u); 
		    j=0; 
		    for (i=0; i<nvar; i++)  
		    { 
				beta(i) = newbeta(i); 
				newbeta(i) = newbeta(i) +  u(i); 
				if (newbeta(i) > maxbeta(i)) newbeta(i) = maxbeta(i); 
				else if (newbeta(i) < -maxbeta(i)) newbeta(i) = -maxbeta(i); 
	    	} 
	    } 
	}  
    loglik(1) = newlk; 
    chinv2(imat, nvar); 
    for (i=0; i<nvar; i++)  
    { 
		beta(i) = newbeta(i)*scale(i); 
		u(i) /= scale(i); 
		imat(i,i) *= scale(i)*scale(i); 
		for (j=0; j<i; j++)  
		{ 
		    imat(i,j) *= scale(i)*scale(j); 
		    imat(j,i) = imat(i,j); 
	    } 
	} 
//    flag = 1000; 
	finish: 
    return join_cols(beta,means); 
} 
 
/*  
** This code is based on the logit_link function from  
** file src/library/stats/src/family.c 
** Part of the R package, http://www.R-project.org 
*/ 
 
vec logit_link(const vec &mu) 
{ 
 
 
    vec rans(mu.n_elem); 
	vec::iterator ransp = rans.begin(); 
    for (vec::const_iterator mup = mu.begin(); mup != mu.end(); ++mup,++ransp) 
    { 
		if (!std::isnan(*mup)) 
		{ 
			*ransp = (*mup < 1.0) ? ((*mup >= DOUBLEEPS) ? log(*mup/(1.0 - *mup)) : MTHRESH) : THRESH; 
		} 
		else 
		{ 
			*ransp = nan(""); 
		} 
    } 
    return rans; 
} 
/*  
** This code is based on the logit_linkinv function from  
** file src/library/stats/src/family.c 
** Part of the R package, http://www.R-project.org 
*/ 
 
vec logit_linkinv(const vec &eta) 
{ 
    vec rans(eta.n_elem); 
	vec::iterator ransp = rans.begin(); 
	double expeta=0; 
    for (vec::const_iterator etap = eta.begin(); etap != eta.end(); ++etap,++ransp) 
	{ 
		if (!std::isnan(*etap)) 
		{ 
			expeta = (*etap < THRESH) ? ( (*etap > MTHRESH) ? exp(*etap)  : 0.0) : MAXEXP; 
	//		if (expeta==MAXEXP)	*ransp = 1.0; 
	//		else *ransp = expeta/(1.0+expeta); 
			*ransp = (expeta==MAXEXP) ? 1.0 : expeta/(1.0+expeta); 
		} 
		else 
		{ 
			*ransp = nan(""); 
		} 
	} 
	return rans; 
 
 
} 
 
/*  
** This code is based on the binomial_dev_resids function from  
** file src/library/stats/src/family.c 
** Part of the R package, http://www.R-project.org 
*/ 
vec logit_mu_eta(const vec &eta) 
{ 
	double expeta,opexp; 
	 
    vec rans(eta.n_elem); 
	vec::iterator ransp = rans.begin(); 
    for (vec::const_iterator etap = eta.begin(); etap != eta.end(); ++etap,++ransp) 
    { 
		if (!std::isnan(*etap)) 
		{ 
			expeta = (*etap < THRESH) ? ( (*etap > MTHRESH) ? exp(*etap)  : 0.0) : MAXEXP; 
			if (expeta==MAXEXP) *ransp = 0.0; 
			else 
			{ 
				opexp = 1.0 + expeta; 
				*ransp = expeta/(opexp * opexp); 
			} 
		} 
		else 
		{ 
			*ransp = nan(""); 
		} 
    } 
	return rans; 
 
} 
 
/*  
** This code is based on the binomial_dev_resids function from  
** file src/library/stats/src/family.c 
** Part of the R package, http://www.R-project.org 
*/ 
vec binomial_dev_resids(const vec &y,const vec &mu,const vec &wt) 
{ 
    int i, n = y.n_elem; 
    double mui, yi;  
    vec rans = y; 
   	for (i = 0; i < n; i++) 
   	{ 
	      mui = mu[i]; 
	      yi = y[i]; 
	      rans[i] = 2.0 * wt[i] * ( ((yi) ? log(yi/mui) : log((1.0-yi)/(1.0-mui))  ) ); 
	} 
 
/* 
    int i, n = y.n_elem, lmu = mu.n_elem, lwt = wt.n_elem; 
    double mui, yi;  
    vec rans = y; 
    if(lmu > 1.0)  
    { 
    	for (i = 0; i < n; i++) 
    	{ 
		      mui = mu[i]; 
		      yi = y[i]; 
		      rans[i] = 2.0 * wt[lwt > 1.0 ? i : 0] * ( ((yi) ? (yi * log(yi/mui)) : 0.0 ) + ((1-yi) ? ( (1-yi) * log((1.0-yi)/(1.0-mui))) : 0.0 )); 
		} 
    }  
    else  
    { 
	  mui = mu[0]; 
	  for (i = 0; i < n; i++)  
	  { 
	    yi = y[i]; 
		rans[i] = 2.0 * wt[lwt > 1.0 ? i : 0] * ( ((yi) ? (yi * log(yi/mui)) : 0.0 ) + ((1-yi) ? ( (1-yi) * log((1.0-yi)/(1.0-mui))) : 0.0 )); 
	  } 
    } 
*/ 
    return rans; 
} 
        
	   
vec modelFittingFunc(const mat &ymat,const mat &XP,const std::string &type) 
{ 
	vec weights = ones<vec>(ymat.n_rows);  
	return modelFittingFunc(ymat,XP,type,weights); 
}	 
    
vec modelFittingFunc(const mat &ymat,const mat &XP,const std::string &type,const vec &weights) 
{ 
  int  nobs = ymat.n_rows; 
  mat X=XP; 
  mat WX = X; 
  vec start = zeros<vec>(X.n_cols); 
    if (X.n_cols==1)   
    { 
        if (X.n_elem==0)  
        { 
        	return start; 
    	} 
		else 
		{ 
			if (type != "COX") 
			{ 
				start = mean(ymat.col(1)); 
				if (type == "LOGIT") 
				{ 
					start = logit_link(start); 
				} 
				return start; 
			} 
		} 
    } 
  vec offset = zeros<vec>(nobs);  
//  vec weights = ones<vec>(nobs);  
  vec newstrat = zeros<vec>(nobs); 
 
  vec strata; 
if (type=="COX") 
{ 
	vec status = ymat.col(1); 
    vec stime = ymat.col(0); 
    uvec sorted; 
    if (strata.n_elem==0)  
    { 
       sorted = sort_index(stime); 
    } 
    else { 
          sorted=sort_index(stime); 
          sorted=sort_index(strata(sorted)); 
          strata=strata(sorted); 
         } 
         for (unsigned int i=0;i<X.n_cols;i++) 
         { 
           vec sor=X.col(i); 
            X.col(i)=sor(sorted); 
         } 
         stime=stime(sorted); 
         status=status(sorted); 
    vec betas=coxfit(20,stime,status,X,offset,weights, newstrat,1,1e-9,std::pow(2.220446e-16,0.75),start,1); 
    return betas; 
} 
  vec y = ymat.col(1); 
  vec mustart=ones<vec>(nobs);   
  vec mu=ones<vec>(nobs);   
  try{ 
  vec eta=ones<vec>(nobs);  
  int iter =0; 
  int maxit=25; 
  double acc=1e-08; 
  double  dev=0.0; 
  double tol = 1.0; 
  vec  n=ones<vec>(nobs); 
  if (type=="LOGIT") 
  { 
    mustart = (weights % y + 0.5)/(weights + 1.0); 
    eta=logit_link(mustart); 
    mu = logit_linkinv(eta + offset); 
    dev=sum(binomial_dev_resids(y, mu, weights)); 
  } 
 	vec varmu(mu.size()); 
	 vec mu_eta_val(eta.size()); 
  if (type=="LM") 
  { 
    maxit=1; 
    mustart=y; 
    eta=mustart; 
    mu=mustart; 
 	dev=sum(weights%(square(y-mu)));   
	varmu.ones(); 
	mu_eta_val.ones();   
   } 
	 double dev0; 
	 vec z=ones<vec>(nobs); 
	 vec W=ones<vec>(nobs); 
	 mat XTX; 
	 vec XTz; 
    while((tol>acc) && (iter<maxit)) 
	{ 
	    iter = iter + 1; 
	    dev0 = dev; 
	    
	    if (type=="LOGIT") 
	    { 
			varmu = mu % (1.0 - mu); 
			mu_eta_val = logit_mu_eta(eta); 
			for (int n=0;n<nobs;n++) 
			{ 
				W[n] = (varmu[n]>0) ?  ((weights[n]*mu_eta_val[n]*mu_eta_val[n])/varmu[n]) : 0.0 ; 
				if (mu_eta_val[n]<DOUBLEEPS) mu_eta_val[n] = DOUBLEEPS; 
			} 
			z =  W % ((eta - offset) + (y - mu)/mu_eta_val); 
			WX = X; 
			WX.each_col() %= W; 
	    } 
		else 
		{ 
			z=y; 
		} 
		 
		XTz = X.t()*z;   
		XTX = X.t()*WX;     
	 
		if ( XTX.has_nan() || XTz.has_nan() )   
		{ 
			for (unsigned int i=0;i<X.n_cols;i++) start[i]=nan(""); 
			return start; 
		} 
	    try 
		{  
			start = pinv(XTX)*XTz; 
		} 
		catch(std::exception const& e) 
		{ 
			Rcpp::Rcout<<"Exception::: "<<e.what()<<"\n"; 
			for (unsigned int i=0;i<X.n_cols;i++) start[i]=nan(""); 
			return start; 
		} 
		if (iter<maxit) 
		{ 
		    eta=X*start; 
		    if (type=="LOGIT") 
		    { 
		          mu = logit_linkinv(eta + offset); 
		          dev=sum(binomial_dev_resids(y, mu, weights)); 
		    } 
		    if (type=="LM") 
		    { 
		        mu=eta;  
		        dev=sum(weights%(square(y-mu))); 
		    } 
			tol = (std::abs(dev0-dev)/(std::abs(dev)+0.1)); 
		} 
  	} 
  return start; 
  } 
    catch(std::runtime_error const& e) 
{ 
    Rcpp::Rcout<<"Exception: "<<e.what()<<"\n"; 
	for (unsigned int i=0;i<X.n_cols;i++) start[i]=nan(""); 
   return start; 
}  
} 
 
vec predictForFresaFunc(const vec &cf,const mat &newdata,const std::string &typ,const std::string &opc)  
{ 
	vec out; 
	if(opc=="COX") 
	{ 
	    vec coef = cf.subvec(0,newdata.n_cols-1); 
	    vec means = cf.subvec(newdata.n_cols,cf.n_elem-1); 
		out = newdata*coef - sum(coef%means); 
	     if(typ =="prob")   
	    { 
//		  out = (1.0/(1.0+exp(-out))+1.0)-1.0; 
		  out = logit_linkinv(out); 
	    } 
	    if(typ =="risk")   
		{ 
			out = exp(out); 
		} 
	} 
	else 
	{	  
		out =newdata*cf;   
	    if(typ =="prob")   
		{ 
//		  out = (1.0/(1.0+exp(-out))+1.0)-1.0; 
		  out = logit_linkinv(out); 
		} 
	} 
	 	 
    return out; 
} 
 
vec improveProbFuncSamples(const vec &x1,const vec &x2,const vec &y, unsigned int samples,double se_nri, double se_idi) 
{ 
	if ((samples == 0) || (x1.n_elem >= samples)) 
	{ 
		return improveProbFunc(x1,x2,y); 
	} 
	else 
	{	 
		const unsigned int trials=6; // set to number of estimations 
		vec imp(6);			 
		mat imppart(imp.n_elem,trials); 
		vec vecx1(samples),vecx2(samples),vecy(samples); 
		vec impt(imp.n_elem); 
		unsigned int i,j,n,nj; 
		unsigned int elem=x1.n_elem; 
		for (i=0;i<trials;i++)	 
		{ 
			unsigned int off = (randi<uvec>(1))[0] % elem; 
			for (n=0;n<samples;n++)	 
			{ 
				nj = (n+off) % elem; 
				vecx1[n]=x1[nj];	 
				vecx2[n]=x2[nj];	 
				vecy[n]=y[nj]; 
			} 
			impt = improveProbFunc(vecx1,vecx2,vecy,se_nri,se_idi); 
			for (j=0;j<impt.n_elem;j++) imppart(j,i) = impt[j]; 
		} 
		imp = median(imppart,1); 
		return imp; 
	} 
} 
 
 
/*  
** This code is based on the improveProb function from Hmisc package  
** Reference: Frank E Harrell Jr. URL: http://biostat.mc.vanderbilt.edu/wiki/Main/Hmisc  
** R package version 2.37-7, http://CRAN.R-project.org/package=Hmisc. 
*/ 
 
vec improveProbFunc(const vec &x1,const vec &x2,const vec &y, double se_nri, double se_idi) 
{ 
	double idi = 0.0; 
	double nri = 0.0; 
	double z_idi = 0.0;	 
	double z_nri = 0.0; 
    vec d  = x2 - x1; 
	 
    vec da=d.elem(find(y==1)); 
    vec db=d.elem(find(y==0)); 
	double n_ev=da.size(); 
	double n_ne=db.size(); 
	 
	if ((n_ne>0)&&(n_ev>0)) 
	{ 
		double nup_ev = 0; 
		double nup_ne = 0;  
		double ndown_ev = 0; 
		double ndown_ne = 0; 
		for(unsigned int i=0;i<y.n_elem;i++) 
		{ 
			if (y[i]==1) 
			{ 
				nup_ev += (d[i]>0); 
				ndown_ev += (d[i]<0); 
			} 
			else 
			{ 
				nup_ne += (d[i]>0); 
				ndown_ne += (d[i]<0); 
			}	 
		} 
		 
		double pup_ev   = nup_ev/n_ev; 
		double pup_ne   = nup_ne/n_ne; 
		double pdown_ev = ndown_ev/n_ev; 
		double pdown_ne = ndown_ne/n_ne; 
		double v_ev = (nup_ev - ndown_ev)/n_ev; 
		double v_ne = (ndown_ne - nup_ne)/n_ne; 
		 
		nri = pup_ev - pdown_ev - (pup_ne - pdown_ne); 
 
		double v_nri_ev = (nup_ev + ndown_ev)/(n_ev*n_ev) - v_ev*v_ev/n_ev; 
		double v_nri_ne = (ndown_ne + nup_ne)/(n_ne*n_ne) - v_ne*v_ne/n_ne; 
		se_nri = std::max(se_nri,std::sqrt(v_nri_ev + v_nri_ne)); 
		if (se_nri < std::numeric_limits< double >::min()) se_nri= std::numeric_limits< double >::min(); 
		z_nri =  nri/se_nri; 
		idi = mean(da) - mean(db); 
		se_idi = std::max(se_idi,std::sqrt(var(da)/n_ev + var(db)/n_ne)); 
		if (se_idi < std::numeric_limits< double >::min()) se_idi = std::numeric_limits< double >::min(); 
		z_idi = idi/se_idi; 
 
	} 
    vec out(6); 
    	out[0]=z_idi; 
    	out[1]=z_nri; 
    	out[2]=idi; 
    	out[3]=nri; 
    	out[4]=se_idi; 
    	out[5]=se_nri; 
    return out; 
  } 
   
getVReclass getVarBinFunc(const mat &dataframe,const std::string &type,const mat &independentFrameP,const mat &bestFrameP,const mat &bestTestFrameP)  
{ 
	const int z_idi=0; 
	const int z_nri=1; 
	const int idi=2; 
	const int nri=3; 
	unsigned int i=3; 
	int nzero=0; 
	int n_var=dataframe.n_cols-3; 
	const mat *independentFrame = &independentFrameP; 
	const mat *bestFrame = &bestFrameP; 
	const mat *bestTestFrame= &bestTestFrameP; 
	if (type=="COX") 
	{ 
		n_var=dataframe.n_cols-2; 
		i=2; 
		nzero=1; 
	} 
	if (independentFrame->n_rows<=1) 
	{ 
		independentFrame = &dataframe; 
	} 
	if (bestFrame->n_rows<=1) 
	{ 
		bestFrame = &dataframe; 
		bestTestFrame = &dataframe; 
	} 
	vec model_zidi(n_var); 
	vec model_znri(n_var); 
	vec model_idi(n_var); 
	vec model_nri(n_var); 
	vec t_model_zidi(n_var); 
	vec t_model_znri(n_var); 
	vec t_model_idi(n_var); 
	vec t_model_nri(n_var); 
	unsigned int samples=dataframe.n_rows; 
	vec FullModel = modelFittingFunc(bestFrame->cols(0,1),bestFrame->cols(2,bestFrame->n_cols-1),type); 
	if (!FullModel.has_nan()) 
	{ 
		vec FullPredict_train = predictForFresaFunc(FullModel,bestFrame->cols(2,bestFrame->n_cols-1),"prob",type); 
		vec FullPredict = predictForFresaFunc(FullModel,bestTestFrame->cols(2,bestTestFrame->n_cols-1),"prob",type); 
		for (int  j=0; i<dataframe.n_cols; i++,j++) 
		{ 
			mat dataframe_sh = dataframe; 
			mat independentFrame_sh = *independentFrame; 
			if (n_var>nzero) 
			{ 
				dataframe_sh.shed_col(i); 
				independentFrame_sh.shed_col(i); 
			} 
			else 
			{ 
				dataframe_sh.fill(1); 
				independentFrame_sh.fill(1); 
			} 
			vec redModel= modelFittingFunc(dataframe_sh.cols(0,1),dataframe_sh.cols(2,dataframe_sh.n_cols-1),type); 
			vec iprob_t = improveProbFunc(predictForFresaFunc(redModel,dataframe_sh.cols(2,dataframe_sh.n_cols-1),"prob",type), FullPredict_train ,dataframe.col(1)); 
			vec iprob = improveProbFuncSamples(predictForFresaFunc(redModel,independentFrame_sh.cols(2,independentFrame_sh.n_cols-1),"prob",type),FullPredict,independentFrame->col(1),samples,iprob_t[5],iprob_t[4]); 
			model_zidi[j]=iprob[z_idi]; 
			model_idi[j]=iprob[idi]; 
			model_nri[j]=iprob[nri]; 
			model_znri[j]=iprob[z_nri]; 
			t_model_zidi[j]=iprob_t[z_idi]; 
			t_model_idi[j]=iprob_t[idi]; 
			t_model_nri[j]=iprob_t[nri]; 
			t_model_znri [j]=iprob_t[z_nri]; 
		} 
	} 
	getVReclass result; 
	 result.z_IDIs=model_zidi; 
	 result.z_NRIs=model_znri; 
	 result.IDIs=model_idi; 
	 result.NRIs=model_nri; 
	 result.tz_IDIs=t_model_zidi; 
	 result.tz_NRIs=t_model_znri; 
	 result.tIDIs=t_model_idi; 
	 result.tNRIs=t_model_nri; 
    return result; 
} 
 
double rocAUC(const vec &response,const vec &predictor,const std::string &direction)  
{ 
	std::string dir = "ascend"; 
	if (direction == "auto") 
	{ 
		vec controls=predictor(find(response==0)); 
		vec cases=predictor(find(response==1)); 
		if (median(controls) < median(cases)) 	dir = "descend"; 
	} 
	else  
	{ 
		if (direction=="<") dir = "descend"; 
	} 
    	 
	double ncases = sum(response); 
	double ncontrol = predictor.n_elem-ncases; 
	 
	uvec sindex = sort_index(predictor,dir.c_str()); 
	double auc=0.0; 
	double sen=0.0; 
	double spe=0.0; 
	double sena=0.0; 
	double spea=0.0; 
	for (unsigned int  i=0;i<predictor.n_elem;i++) 
	{ 
		sen = sen + response[sindex[i]]/ncases; 
		spe = spe + (response[sindex[i]]==0)/ncontrol; 
		auc = auc + (sen+(sen-sena)/2.0)*(spe-spea); 
		sena = sen; 
		spea = spe; 
	} 
 
	return(auc); 
} 
 
 
vec Fresarank(const vec &xi) 
{ 
    vec x(xi.n_elem); 
	vec so=sort(xi); 
	unsigned int n=xi.n_elem; 
    uvec si_x=sort_index(xi); 
    vec sortRankx(n+1); 
    unsigned int s; 
	for(unsigned int i=0;i<n;) 
	{ 
		s=i+1; 
		unsigned int h=1; 
		if ( (i+h) < n ) 
		{ 
			while(so(i)==so(i+h)) 
			{ 
				s=s+(i+1+h); 
				h++; 
				if (i+h == n) break; 
			} 
		} 
		double in = s/(double)(h); 
		unsigned int t=i+h; 
		for (;i<t;i++) sortRankx(i)=in; 
	} 
	for(unsigned int i=0;i<n;i++) x(si_x(i))=sortRankx(i); 
	return(x); 
} 
 
vec residualForFRESAFunc(const vec &cf,const mat &newdata,const std::string &typ,const std::string &type,const mat &outcome)  
{ 
	vec out; 
	if(type=="COX") // returns the log partial likelihood 
	{ 
		out=predictForFresaFunc(cf,newdata,"prob",type) - outcome.col(1); 
	} 
    if(type=="LM") 
	{ 
		out=predictForFresaFunc(cf,newdata,"linear",type) - outcome.col(1); 
	} 
    if(type=="LOGIT") 
    { 
		if (typ == "linear") 
		{ 
			out=predictForFresaFunc(cf,newdata,"linear",type) - outcome.col(1); 
		} 
		else 
		{ 
			out=predictForFresaFunc(cf,newdata,"prob",type) - outcome.col(1); 
		} 
	} 
	if (out.has_nan())  
	{ 
		Rcpp::Rcout<<"Warning NA predictFor NeRIs \n"; 
		out.elem(find_nonfinite(out)).fill(1.0e10); 
	} 
    return (out); 
} 
 
improvedRes improvedResidualsFunc(const vec &oldResiduals,const vec &newResiduals,const std::string &testType,unsigned int samples) 
{ 
	if (oldResiduals.n_elem >= samples) 
	{ 
		vec ored=oldResiduals; 
		vec nred=newResiduals; 
		if ((samples>0)&&(oldResiduals.n_elem > samples)) 
		{ 
			uvec samSamples = sort_index(randu(oldResiduals.n_elem));	 
			samSamples.resize(samples); 
			ored=oldResiduals(samSamples); 
			nred=newResiduals(samSamples); 
		} 
		improvedRes imp = improvedResidualsFunc(ored,nred,testType); 
		return imp; 
	} 
	else 
	{	 
		unsigned int trials=7; // set to number of estimations 
		improvedRes impt;			 
		mat imppart(8,trials); 
		vec vecx1(samples),vecx2(samples); 
		unsigned int i; 
		unsigned int elem=oldResiduals.n_elem; 
		uvec samSamples(elem); 
		for (i=0;i<trials;i++)	 
		{ 
 
			samSamples = sort_index(randu(elem)); 
			for (unsigned int n=0;n<samples;n++) 
			{ 
				unsigned int nj = n % elem; 
				vecx1[n]=oldResiduals[samSamples[nj]];	 
				vecx2[n]=newResiduals[samSamples[nj]];	 
			} 
 
			impt = improvedResidualsFunc(vecx1,vecx2,testType); 
			imppart(0,i)=impt.p1; 
			imppart(1,i)=impt.p2; 
			imppart(2,i)=impt.NeRI; 
			imppart(3,i)=impt.pvalue; 
			imppart(4,i)=impt.binom_pValue; 
			imppart(5,i)=impt.wilcox_pValue; 
			imppart(6,i)=impt.t_test_pValue; 
			imppart(7,i)=impt.F_test_pValue; 
		} 
		 
		vec medimp = median(imppart,1); 
		impt.p1=medimp[0]; 
		impt.p2=medimp[1]; 
		impt.NeRI=medimp[2]; 
		impt.pvalue = medimp[3]; 
		impt.binom_pValue = medimp[4]; 
		impt.wilcox_pValue = medimp[5]; 
		impt.t_test_pValue = medimp[6]; 
		impt.F_test_pValue = medimp[7]; 
		 
		return impt; 
	} 
} 
 
improvedRes improvedResidualsFunc(const vec &oldResiduals,const vec &newResiduals,const std::string &testType) 
{ 
    improvedRes result; 
 
 
	double pwil = 0.5; 
	double pbin = 0.5; 
	double ptstu = 0.5; 
	double f_test = 0.5; 
	double p1=0.0; 
	double p2=0.0; 
	double pvalue = 0.5; 
	double size = oldResiduals.n_elem; 
	double improved = 0.0; 
	if (size==0)  
	{ 
		Rcpp::Rcout<<"Zero Elements:ImproveResiduals \n"; 
	} 
	vec oldres = abs(oldResiduals); 
	vec newres = abs(newResiduals); 
	double reduction = 0.0; 
	double increase  = 0.0; 
	for (int i=0; i< size;i++) 
	{ 
		reduction += (newres[i]<oldres[i]); 
		increase += (newres[i]>oldres[i]); 
	} 
	p1 = reduction/size;															//	#proportion of subjects with improved residuals 
	p2 = increase/size;											//#proportion of subjects with worst residuals 
	improved = p1-p2; 									//			# the net improvement in residuals 
 
 
//					Rcout << oldResiduals; 
//					Rcout << "....\n"; 
//					Rcout << newResiduals; 
//					Rcout << "\n" << reduction << " : " << increase << "\n"; 
 
	double rss1 = sum(square(oldResiduals)); 
	 
	if (rss1 > 1000*DOUBLEEPS) 
	{ 
		pvalue = 0.5; 
		pwil = pvalue; 
		pbin =  pvalue; 
		ptstu =  pvalue; 
		f_test =  pvalue; 
		double rss2 = sum(square(newResiduals)); 
		if (rss2 > DOUBLEEPS) 
		{ 
//					Rcout << "\n 1:" << rss1 << " 2:" << rss2 << "\n"; 
			int sw=0; 
			std::string tail="greater"; 
 
			if (testType=="Binomial") sw=1; else 
			if (testType=="Ftest") sw=2; else 
			if (testType=="Wilcox") sw=3; else 
			if (testType=="tStudent") sw=4; 
			switch (sw) 
			{ 
				case 1: 
				{ 
					pbin = 0.5; 
					if (reduction>=increase) pbin = binomtest(reduction,size,0.5,tail); 
					pvalue = pbin;	 
					break; 
				} 
				case 2: 
				{ 
					rss2 = rss1/rss2; 
//					Rcout << "\n 1:" << rss2 << "\n"; 
					if (rss2>10000.0) rss2=10000.0; 
					pvalue = f_test = (1.0-R::pf(size*(rss2-1.0),1.0,size,1,0))*0.5; 
					break; 
				} 
				case 3: 
				{ 
					pvalue = pwil = wilcoxtest(oldres, newres,0.0,TRUE,tail,TRUE); 
					break; 
				} 
				case 4: 
				{ 
//					pvalue = ptstu = ttest(oldres, newres,0.0,TRUE,TRUE,tail); 
					pvalue = ptstu = pttest(oldres, newres,tail); 
					break; 
				} 
				default: 
				{ 
					rss2 = rss1/rss2; 
					if (rss2>10000.0) rss2=10000.0; 
//					Rcout << "\n 1:" << rss2 << "\n"; 
					f_test = (1.0-R::pf(size*(rss2-1.0),1.0,size,1,0))*0.5; 
					pvalue =  f_test; 
					pbin=0.5; 
					pwil=1.0; 
					if (reduction>=increase) pbin = binomtest(reduction,size,0.5,tail); 
					if (reduction>=increase) pwil = wilcoxtest(oldres, newres,0.0,TRUE,tail,TRUE); 
//					ptstu = ttest(oldres, newres,0.0,TRUE,TRUE,tail); 
					ptstu = pttest(oldres, newres,tail); 
					break; 
				} 
			} 
		} 
	} 
 
//	Rcout << "\n p:" << pvalue << " f: " << f_test << "\n"; 
 
	result.p1 = p1; 
	result.p2 = p2; 
	result.NeRI = improved; 
	result.pvalue = pvalue; 
	result.binom_pValue = pbin; 
	result.wilcox_pValue = pwil; 
	result.t_test_pValue = ptstu; 
	result.F_test_pValue = f_test; 
 
 
	 
	return (result); 
} 
 
double pttest(const vec &xt, const vec &y,const std::string &tail) 
{ 
    double pval=1.0; 
	vec x=xt-y; 
    double nx = x.n_elem; 
    double mx = mean(x); 
    double vy = std::min(var(x)/2,var(y)); 
    double df = nx-1.0; 
    double tstat; 
    double stderrr;  
	if(nx < 2.0)  
	{ 
		tstat = 0.0;	// let's make it zero 
	} 
	else 
	{ 
		stderrr = std::sqrt(vy/nx); 
		if(stderrr <= (10.0 * DOUBLEEPS)) 
		{ 
			tstat = 10.0*(mx > 0); // let's make it 10 
		} 
		else 
		{ 
			tstat = mx/stderrr; 
		} 
    } 
    if (tail == "less")  
    { 
		pval = R::pt(tstat, df,1,0); 
    } 
    else  
    	if (tail == "greater")  
    	{ 
			pval = R::pt(tstat, df, 0,0); 
    	} 
   		else  
			{ 
				pval = 2.0 * R::pt(-std::abs(tstat), df,1,0); 
			} 
    return(pval); 
}	 
 
 
/*  
** This code is based on the t.test function from  
** file src/library/stats/R/t.test.R 
** Part of the R package, http://www.R-project.org 
*/ 
 
double ttest(const vec &xt, const vec &y, double mu, bool paired, bool var_equal,const std::string &tail) 
{ 
    double pval=1.0; 
	vec x=xt; 
    if (paired)  
    { 
		x = x-y; 
    } 
    double nx = x.n_elem; 
    double mx = mean(x); 
    double vx = var(x); 
    double df=1.0; 
    double tstat; 
    double stderrr;  
    if(y.empty() || paired)  
    { 
        if(nx < 2.0)  
		{ 
//			stop("not enough 'x' observations"); 
			tstat = 0.0;	// let's make it zero 
		} 
		else 
		{ 
			df = nx-1.0; 
			stderrr = std::sqrt(vx/nx); 
				if(stderrr <= (10.0 * DOUBLEEPS)) 
				{ 
//		            stop("data are essentially constant"); 
					tstat = 10.0*(mx != mu); // if equal return 0, else 10 
				} 
			else 
			{ 
				tstat = (mx-mu)/stderrr; 
			} 
		} 
    }  
    else  
    { 
		double ny = y.n_elem; 
		double my = mean(y); 
		double vy = var(y); 
		 
		if(var_equal)  
		{ 
		    df = nx+ny-2.0; 
	            double v = 0.0; 
	            if(nx > 1) v = v + (nx-1.0)*vx; 
	            if(ny > 1) v = v + (ny-1.0)*vy; 
		    v = v/df; 
		    stderrr = std::sqrt(v*(1.0/nx+1.0/ny)); 
		}  
		else  
		{ 
		    double stderrx = std::sqrt(vx/nx); 
		    double stderry = std::sqrt(vy/ny); 
		    stderrr = std::sqrt(stderrx*stderrx + stderry*stderry); 
		    df = std::pow(stderrr,4)/(std::pow(stderrx,4)/(nx-1.0) + std::pow(stderry,4)/(ny-1.0)); 
		} 
	        if(stderrr <= 10.0 *DOUBLEEPS * std::max(std::abs(mx), std::abs(my))) 
			{ 
//	            stop("data are essentially constant"); 
				tstat=0.0; 
			} 
			else 
			{ 
				tstat = (mx - my - mu)/stderrr; 
			} 
    } 
    if (tail == "less")  
    { 
		pval = R::pt(tstat, df,1,0); 
    } 
    else  
    	if (tail == "greater")  
    	{ 
			pval = R::pt(tstat, df, 0,0); 
    	} 
   		else  
			{ 
				pval = 2.0 * R::pt(-std::abs(tstat), df,1,0); 
			} 
    return(pval); 
} 
/*  
** This code is based on the wilcox.test function from  
** file src/library/stats/R/wilcox.test.R 
** Part of the R package, http://www.R-project.org 
*/ 
 
double wilcoxtest(const vec &xt,const vec &y , double mu, bool paired,const std::string &tail,bool correct) 
{ 	 
    double pvalue=1.0; 
	vec x=xt; 
    double  correc = 0.0; 
    if(y.empty() || paired )  
    { 
        if (!y.empty()) x = x - y; 
        x = x - mu; 
        bool zeros = any(x == 0.0); 
        if (zeros) x = x.elem(find(x != 0.0)); 
        double n = x.n_elem; 
		if (n>1.0) 
		{ 
			vec r = Fresarank(abs(x)); 
			double stat = sum(r.elem(find(x > 0.0))); 
			{  
				std::map <double, double> rtiesm; 	 
                for(unsigned int i=0;i<r.n_elem;i++) 
                { 
                	rtiesm[r[i]]++; 
                } 
                vec rtiesv(rtiesm.size()); 
                std::map<double,double>::iterator ii=rtiesm.begin(); 
                 for(unsigned int i=0 ; ii!=rtiesm.end(); ++ii) 
				   { 
				   	   rtiesv[i]=(*ii).second; 
				   	   i++; 
				   } 
				double z = stat - n * (n + 1.0)/4.0; 
				double sigma = std::sqrt(n * (n + 1.0) * (2.0 * n + 1.0)/24.0 - sum(pow(rtiesv,3) - rtiesv)/48.0); 
				if(correct)  
				{ 
				   if(tail=="greater") correc=0.5; 
					else if(tail=="less") correc=-0.5; 
						else  
						{ 
							if (z!=0) correc=(z/std::abs(z)) * 0.5; 
								else correc=0.5; 
						} 
				} 
				if (sigma>0.0) z = (z - correc)/sigma; else z=10.0*(z != correc); 
		   
				if(tail=="less") pvalue= R::pnorm5(z,0.0,1.0,1,0); 
					else if(tail=="greater") pvalue = R::pnorm5(z,0.0,1.0,0,0); 
						else pvalue = 2.0 * std::min(R::pnorm5(z,0.0,1.0,1,0),R::pnorm5(z,0.0,1.0,0,0)); 
			} 
		} 
    } 
    else  
    {  
    //-------------------------- 2-sample case --------------------------- 
        vec r = Fresarank(x - y); 
		x.clear(); 
        double nx = x.n_elem; 
        double ny = y.n_elem; 
        bool exact= ((nx < 50) && (ny < 50)); 
 
        double stat =  sum(r) - nx * (nx + 1.0) / 2.0; 
        vec ur=unique(r); 
        bool ties=(r.n_elem != ur.n_elem); 
		ur.clear(); 
		 
        if(exact && !ties)  
        { 
        	 if(tail=="greater") 
                pvalue =R::pwilcox(stat - 1, nx, ny, 0,0); 
             else if(tail=="less") pvalue =R::pwilcox(stat, nx, ny,1,0); 
             else 
                   { 
                   	double p; 
                       if(stat > (nx * ny / 2.0)) 
                           p =R::pwilcox(stat - 1, nx, ny, 0,0); 
                       else 
                           p =R::pwilcox(stat, nx, ny,1,0); 
                        pvalue =std::min(2.0 * p, 1.0); 
                   } 
        } 
        else  
        { 
			std::map <double, double> rtiesm; 	 
            for(unsigned int i=0;i<r.n_elem;i++) 
            { 
                rtiesm[r[i]]++; 
            } 
            vec rtiesv(rtiesm.size()); 
            std::map<double,double>::iterator ii=rtiesm.begin(); 
            for(unsigned int i=0 ; ii!=rtiesm.end(); ++ii) 
			{ 
			   rtiesv[i]=(*ii).second; 
				i++; 
			} 
            double z = stat - nx * ny/2.0; 
            double sigma = std::sqrt((nx * ny/12.0) * ((nx + ny + 1)- sum(pow(rtiesv,3) - rtiesv)/((nx + ny) * (nx + ny - 1.0)))); 
            if(correct)  
            { 
               if (tail=="greater") correc=0.5; 
	           else if (tail=="less") correc=-0.5; 
	           else  
	           { 
		           	if (z!=0) 
		            	correc=(z/std::abs(z)) * 0.5; 
		            else 
		            	correc=0.5; 
		        } 
            } 
		    z = (z - correc) / sigma; 
		    if (tail=="less")pvalue= R::pnorm5(z,0.0,1.0,1,0); 
				else if (tail=="greater")pvalue = R::pnorm5(z,0.0,1.0,0,0); 
					else pvalue = 2 * std::min(R::pnorm5(z,0.0,1.0,1,0),R::pnorm5(z,0.0,1.0,0,0)); 
        } 
    } 
    return(pvalue); 
} 
 
double binomtest(double x, double n, double p ,const std::string &tail) 
{ 
//    x = R::fround(x,0); 
//    n = R::fround(n,0); 
    double pvalue=1.0; 
    if (tail=="less") pvalue =R::pbinom(x, n, p,1,0); 
    else if (tail=="greater") pvalue =R::pbinom(x-1, n, p, 0,0); 
	return(pvalue); 
} 
 
gvarNeRI getVarResFunc(const mat &dataframe,const std::string &type,const mat &testdataP,int testsamples,const mat &fullFrameP,const mat &fulltestFrameP)  
{ 
	gvarNeRI result; 
	const mat *testdata = &testdataP; 
	if (testdata->n_elem <= 1) 
	{ 
		testdata = &dataframe; 
	} 
	const mat *fulldataFrame = &fullFrameP; 
	const mat *fulltestFrame = &fulltestFrameP; 
	if (fulldataFrame->n_elem <=1 ) 
	{ 
		fulldataFrame = &dataframe; 
		fulltestFrame = &dataframe; 
	} 
	unsigned int i=3; 
	int nzero=0; 
	int nvar=dataframe.n_cols-3; 
	if (type=="COX") 
	{ 
		nvar=dataframe.n_cols-2; 
		nzero=1; 
		i=2; 
	} 
	 
//	Rcout<<"Train Rows"<<dataframe.n_rows<<" Test:"<<testdataP.n_rows<<" Train Rows: "<<fullFrameP.n_rows<<"Test Rows:"<<fulltestFrameP.n_rows<<"\n"; 
 
	vec FullModel = modelFittingFunc(fulldataFrame->cols(0,1),fulldataFrame->cols(2,fulldataFrame->n_cols-1),type); 
	vec FullResiduals = residualForFRESAFunc(FullModel,fulldataFrame->cols(2,fulldataFrame->n_cols-1),"",type,fulldataFrame->cols(0,1)); 
	vec testResiduals = residualForFRESAFunc(FullModel,fulltestFrame->cols(2,fulltestFrame->n_cols-1),"",type,fulltestFrame->cols(0,1)); 
	vec model_tpvalue(nvar); 
	vec model_bpvalue(nvar); 
	vec model_neri(nvar); 
	vec model_wpvalue(nvar); 
	vec model_fpvalue(nvar); 
	vec testmodel_tpvalue(nvar); 
	vec testmodel_bpvalue(nvar); 
	vec testmodel_neri(nvar); 
	vec testmodel_wpvalue(nvar); 
	vec testmodel_fpvalue(nvar); 
	unsigned int samples = dataframe.n_rows; 
	if (testsamples>0) samples=testsamples; 
	if (nvar>0) 
	{ 
		for (int  j=0; i<dataframe.n_cols; i++,j++) 
		{ 
			mat dataframe_sh = dataframe; 
			mat testdata_sh = *testdata; 
		 
			if (nvar>nzero) 
			{ 
				dataframe_sh.shed_col(i); 
				testdata_sh.shed_col(i); 
			} 
			else 
			{ 
				dataframe_sh.fill(1); 
				testdata_sh.fill(1); 
			} 
			vec redModel = modelFittingFunc(dataframe_sh.cols(0,1),dataframe_sh.cols(2,dataframe_sh.n_cols-1),type); 
			vec redResiduals = residualForFRESAFunc(redModel,dataframe_sh.cols(2,dataframe_sh.n_cols-1),"",type,dataframe_sh.cols(0,1)); 
			vec redTestResiduals = residualForFRESAFunc(redModel,testdata_sh.cols(2,testdata_sh.n_cols-1),"",type,testdata_sh.cols(0,1)); 
			improvedRes iprob = improvedResidualsFunc(redResiduals,FullResiduals," ",samples); 
			improvedRes testiprob = improvedResidualsFunc(redTestResiduals,testResiduals," ",samples); 
			model_tpvalue[j] = iprob.t_test_pValue; 
			model_bpvalue[j] = iprob.binom_pValue; 
			model_wpvalue[j] = iprob.wilcox_pValue; 
			model_fpvalue[j] = iprob.F_test_pValue; 
			model_neri[j] = iprob.NeRI; 
			testmodel_tpvalue[j] = testiprob.t_test_pValue; 
			testmodel_bpvalue[j] = testiprob.binom_pValue; 
			testmodel_wpvalue[j] = testiprob.wilcox_pValue; 
			testmodel_fpvalue[j] = testiprob.F_test_pValue; 
			testmodel_neri[j] = testiprob.NeRI; 
		} 
	} 
//	Rcout<<" t-pvalue: "<<model_tpvalue[0]<<"\n"; 
	 result.tP_value=model_tpvalue; 
	 result.BinP_value=model_bpvalue; 
	 result.WilcoxP_value=model_wpvalue; 
	 result.FP_value=model_fpvalue; 
	 result.NeRIs=model_neri; 
	 result.testData_tP_value=testmodel_tpvalue; 
	 result.testData_BinP_value=testmodel_bpvalue; 
	 result.testData_WilcoxP_value=testmodel_wpvalue; 
	 result.testData_FP_value=testmodel_fpvalue; 
	 result.testData_NeRIs=testmodel_neri; 
     return (result); 
} 
 
//**************** it will atempt to create a sampled data set with uniform desity distribution of the output 
//* a running window is used to estimate the density of the data. At regions with low desity the data will be duplicated 
// int sort_index the index of the feature to be equilized 
// int breakWidth the width of the sample window in the data. The number of samples per window will be (size/breakWidth) + 1 
// this funtion limits the maximun sample duplication to 4 
//************************** 
uvec equSamples_o(const mat &inputsample, unsigned int sort_indx, int breakWidth, double range) 
{ 
	vec outcome = inputsample.col(sort_indx); 
	unsigned int osize = outcome.size(); 
	double omin = min(outcome); 
	double omax = max(outcome); 
	if (range==0) range = omax-omin; 
	double delta = range/(10.0*osize); 
	vec histo=zeros<vec>(breakWidth); 
	double histowidth=(range+delta)/breakWidth; 
	 
	vec noiseoutcome = outcome + delta*(randu(osize)-0.5); // we add noise in order to break sample ties 
 
	delta = range/osize; 
	 
	uvec indices = sort_index(noiseoutcome); 
	uvec uindices(3*osize); // we just assume that all the samples will duplicated four times 
	int sdx=osize/(2*breakWidth); // there will be 2*sdx+1 samples at each density estimation 
	if (sdx<2) sdx=2; // the miniumun number of samples for desity estimation  
	 
	int indx = 0; 
	for (unsigned int i=0; i < osize;i++) 
	{ 
		int hind=(int)((outcome[i]-omin)/histowidth); 
		++histo[hind];  
	} 
	int maxhist=max(histo); 
	for (unsigned int i=0; i < osize;i++) 
	{ 
		int hind=(int)((outcome[i]-omin)/histowidth); 
		int ed1 = i; 
		int ed2 = osize-i-1; 
		int dsdx=sdx; 
		if (ed2<ed1) ed1=ed2; 
		if (dsdx>ed1) dsdx=ed1; 
		if (dsdx<2) dsdx=2;				// We limit the minimum width to 2 at sample extreme 
		int indx1 = (i-dsdx); 
		int indx2 = (i+dsdx); 
		if (indx1<0) indx1 = 0; 
		if (indx2 >= (int)osize) indx2 = osize-1; 
		int idx = indx2-indx1; 
		int nsamplesdx = (int)((noiseoutcome[indices[indx2]]-noiseoutcome[indices[indx1]])/(idx*delta)+0.5); 
		if (nsamplesdx>3) nsamplesdx=3;  // the maxumun number of replications is 2 
		if (nsamplesdx<1) nsamplesdx=1;   // we copy every sample regardless of high density 
		if (histo[hind]==maxhist) nsamplesdx=1; 
		for (int j=0; j< nsamplesdx;j++,indx++) 
		{ 
			uindices[indx]=indices[i]; 
		} 
// 		Rcout << i << " last:" << indx << " delta indx " << idx << " Samples "<< nsamplesdx << " data: " << outcome[i] << "\n"; 
	} 
	uindices.resize(indx); 
//	Rcout << " Range:" << range << " Delta:" << delta << " Histo[1]:" << histo[0] << "Histo[n]" << histo[breakWidth-1] << "\n";  
	if (((histo[0]+histo[breakWidth-1]) > 0.75*osize) && (std::abs(histo[0]-histo[breakWidth-1])/osize < 0.01)) uindices = indices; 
	return uindices; 
} 
 
//**************** it will atempt to create a sampled data set with uniform desity distribution of the output 
//************************** 
uvec equSamples(const mat &inputsample, unsigned int sort_indx, int histbins,double omin, double range) 
{ 
	vec outcome = inputsample.col(sort_indx); 
	unsigned int osize = outcome.size(); 
	if (range==0)  
	{ 
		omin = min(outcome); 
		double omax = max(outcome); 
		range = omax-omin; 
	} 
	double delta = range/(10.0*osize); 
	vec histo=zeros<vec>(histbins); 
	double histowidth=(range+delta)/histbins; 
	 
	delta = range/osize; 
	 
	uvec indices(osize); 
	uvec uindices(3*osize); // we just assume that all the samples will duplicated four times 
	 
	int indx = 0; 
	for (unsigned int i=0; i < osize;i++) 
	{ 
		int hind=(int)((outcome[i]-omin)/histowidth); 
		++histo[hind];  
		indices[i]=i; 
	} 
	double maxhist=max(histo); 
	vec prand = randu(osize); 
	for (unsigned int i=0; i < osize;i++) 
	{ 
		int hind=(int)((outcome[i]-omin)/histowidth); 
		double frac= maxhist/histo[hind]-1.0; 
		int nsamplesdx=1; 
		if (frac>0) 
		{ 
			if (frac>1)  
			{ 
				nsamplesdx = (int)(frac); 
				if (nsamplesdx>3) nsamplesdx=3; 
			} 
			else 
			{ 
				nsamplesdx=1+1*(prand[i]<frac); 
			} 
		} 
		for (int j=0; j< nsamplesdx;j++,indx++) 
		{ 
			uindices[indx]=indices[i]; 
		} 
//		Rcout  << i  << " Histo[i]:" << histo[hind] << " Frac:" <<frac<< " Samples:" << nsamplesdx << "\n";  
	} 
	uindices.resize(indx); 
	return uindices; 
} 
 
