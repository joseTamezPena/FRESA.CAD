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

struct bootVal{
  	 mat bootmodel;
  	 mat bootmodelmeans;
	 double BlindAccuracy;
	 double BlindSensitivity;
	 double BlindSpecificity;
	 std::vector<double> testoutcome;
	 std::vector<double> testprediction;
	 mat zNRI;
	 mat zIDI;
	 mat test_zIDI;
	 mat test_zNRI;
  }redBoot;

int redCurmodel_S_lastRemoved;


//********************************************************************************************
//**==========================================================================================
//********************************************************************************************
extern "C" SEXP bootstrapValidationBinCpp(SEXP _fraction,SEXP _loops,SEXP _dataframe,SEXP _type,SEXP _response,SEXP _bestdataframe)
{
	try 
	{  //R_CStackLimit=(uintptr_t)-1;
//		arma_rng::set_seed_random();
		double fraction = Rcpp::as<double>(_fraction);
		int loops = Rcpp::as<int>(_loops);
		std::string type = Rcpp::as<std::string>(_type);
		Rcpp::NumericMatrix dataf(_dataframe);
		Rcpp::NumericMatrix resp(_response);
		Rcpp::NumericMatrix bestX(_bestdataframe);

	    mat dataframe;
		mat bestXframe;

	    mat Response(resp.rows(), 2);
		mat dat(dataf.begin(), dataf.rows(), dataf.cols(), false);
		mat bestdatX(bestX.begin(), bestX.rows(), bestX.cols(), false);

		vec dtime=ones<vec>(Response.n_rows);
		if (resp.cols()==1) 
			{
				Response.col(0)=dtime;
				Response.col(1)=Rcpp::as<vec>(Rcpp::NumericVector(_response));
			}
			else 
			{
			 mat Res(resp.begin(), resp.rows(), resp.cols(), false);
			 Response=Res;
			}
		vec Outcome=Response.col(1);
		dtime=Response.col(0);
        int n_var2=dat.n_cols-1;
        if (type=="COX")
		{
			dataframe=dat.cols(1,dat.n_cols-1);
			bestXframe=bestdatX.cols(1,bestdatX.n_cols-1);
		}
		else 
		{
			dataframe=dat;
			bestXframe=bestdatX;
		}
		int n_var=dataframe.n_cols;

	    mat casesample=join_rows(Response.rows(find(Outcome==1)),dataframe.rows(find(Outcome==1)));
	    mat controlsample=join_rows(Response.rows(find(Outcome==0)),dataframe.rows(find(Outcome==0)));

	    mat bestCaseSample = join_rows(Response.rows(find(Outcome==1)),bestXframe.rows(find(Outcome==1)));
	    mat bestControlSample = join_rows(Response.rows(find(Outcome==0)),bestXframe.rows(find(Outcome==0)));
		
		
		int sizecases   =casesample.n_rows;
		int sizecontrol =controlsample.n_rows;
		int minsize = std::min(sizecases,sizecontrol);
		int totSamples = (int)(fraction*minsize);

		std::vector<double> testoutcome;
		std::vector<double> testprediction;

		vec trainRoc=ones<vec>(loops);
		vec accuracy=zeros<vec>(loops);
		vec sensitivity=zeros<vec>(loops);
		vec specificity=zeros<vec>(loops);
		vec taccuracy=ones<vec>(loops);
		vec tsensitivity=zeros<vec>(loops);
		vec tspecificity=zeros<vec>(loops);
		mat bcoef=zeros<mat>(loops,n_var); 	  
		mat zNRI=zeros<mat>(loops,n_var2); 
		mat zIDI=zeros<mat>(loops,n_var2); 
		mat test_zNRI=zeros<mat>(loops,n_var2); 
		mat test_zIDI=zeros<mat>(loops,n_var2); 
		mat NRI=zeros<mat>(loops,n_var2); 
		mat IDI=zeros<mat>(loops,n_var2); 
		double gacc = 0.0; 
		double gsen = 0.0; 
		double gspe = 0.0; 
		int smaptot = 0;
		int smapsen = 0;
		int smapspe = 0;
		int lastInserted=0;

		vec rawbetas = modelFittingFunc(Response,dataframe,type);

		
		#pragma omp parallel for schedule(dynamic) ordered shared(lastInserted,testoutcome,testprediction,gacc,gsen,gspe,smaptot,smapsen,smapspe,trainRoc,bcoef,accuracy,sensitivity,specificity,taccuracy,tsensitivity,tspecificity,zNRI,zIDI,test_zNRI,test_zIDI,NRI,IDI)  
		for (int doOver=0;doOver<loops;doOver++)
		{ 
			mat trainingSample;
			mat bestTrainingSample;
			mat myTestCases;
			mat myTestControl;
			mat bestTestCases;
			mat bestTestControl;
			uvec nsamCaseTest;
			uvec nsamControlTest;
			unsigned int minTest=2; // at least 2 samples for test evaluation
			uvec samCases;
			uvec samControl;
			do
			{
				samCases= randi<uvec>(totSamples, distr_param(0,sizecases-1));
				samControl=randi<uvec>(totSamples, distr_param(0,sizecontrol-1));
				vec auxcasesample=zeros<vec>(sizecases);
				vec auxcontrolsample=zeros<vec>(sizecontrol);
				for (int i =0;i<totSamples;i++)
				{
					auxcasesample[samCases(i)]=1;
					auxcontrolsample[samControl(i)]=1;
				}
				nsamCaseTest = find(auxcasesample==0);
				nsamControlTest = find(auxcontrolsample==0);
			} while ((nsamCaseTest.n_elem<minTest)||(nsamControlTest.n_elem<minTest));

			trainingSample=join_cols(casesample.rows(samCases),controlsample.rows(samControl));
			bestTrainingSample=join_cols(bestCaseSample.rows(samCases),bestControlSample.rows(samControl));
			
			
		    vec trainmodel = modelFittingFunc(trainingSample.cols(0,1),trainingSample.cols(2,trainingSample.n_cols-1)+100.0*DOUBLEEPS*randn(trainingSample.n_rows,trainingSample.n_cols-2),type);
			if (!trainmodel.has_nan())
			{
				// the model beta coefficients
				vec coef=trainmodel;
				if (type == "COX") 
				{
					  vec coefaux = trainmodel.subvec(0,(trainmodel.n_elem/2)-1);
   					  coef=coefaux;
				}


				myTestCases=casesample.rows(nsamCaseTest);
				myTestControl=controlsample.rows(nsamControlTest);
				bestTestCases=bestCaseSample.rows(nsamCaseTest);
				bestTestControl=bestControlSample.rows(nsamControlTest);


				int ncases = myTestCases.n_rows;
				int ncontrol = myTestControl.n_rows;
				int ntesmin = std::max(ncases,ncontrol);
				int randnum = (randi<uvec>(1))[0];
				
// forcing to have the same number of cases and controls
				mat testSample(2*ntesmin,myTestCases.n_cols);
				mat bestTestSample(2*ntesmin,bestTestCases.n_cols);
				for (int i=0;i<ntesmin;i++)
				{
					testSample.row(2*i) =  myTestCases.row((i + randnum) % ncases);
					testSample.row(2*i+1) = myTestControl.row((i + randnum) % ncontrol);
					bestTestSample.row(2*i) = bestTestCases.row((i + randnum) % ncases);
					bestTestSample.row(2*i+1) = bestTestControl.row((i + randnum) % ncontrol);
				}
				ncases = ntesmin;
				ncontrol = ntesmin;
				

				
				// predicting the test
				vec p = predictForFresaFunc(trainmodel,testSample.cols(2,testSample.n_cols-1),"prob",type);
				// predicting the training
				vec trainmodel_predictors = predictForFresaFunc(trainmodel,trainingSample.cols(2,trainingSample.n_cols-1),"prob",type);

				// performance AUC and NRI/IDI
				double auul = rocAUC(trainingSample.col(1), trainmodel_predictors,"auto");
				getVReclass modelReclas = getVarBinFunc(trainingSample,type,testSample,bestTrainingSample,bestTestSample);
				
				double acc = 0.0;
				double sen = 0.0;
				double spe = 0.0;
				for(unsigned int i=0;i<testSample.n_rows;i++)
				{
					if (testSample(i,1)>0)
					{
				    	if (p(i)>=0.5) sen++;
					}
					if (testSample(i,1)==0)
					{
				    	if (p(i)<0.5) spe++;
					}	
				}
				acc = sen+spe;

				double trsen=0;
				double trspe=0;
				double tracc=0;
				for(unsigned int i=0;i<trainingSample.n_rows;i++)
				{
					if (trainingSample(i,1)>0)
					{
					 	if (trainmodel_predictors(i)>=0.5) trsen++;
					}
					if (trainingSample(i,1)==0)
					{
				    	if (trainmodel_predictors(i)<0.5) trspe++;
					}	
				}
				tracc = trsen+trspe;
				tracc = tracc/trainingSample.n_rows;
				trsen = trsen/totSamples;
				trspe = trspe/totSamples;
				
				
#pragma omp critical
{
				for (unsigned int i=0;i<p.n_elem;i++)
				{
					testoutcome.push_back(testSample(i,1));
					testprediction.push_back(p(i));
				}
				gacc = gacc + acc;
				smaptot = smaptot + testSample.n_rows;
				gsen = gsen + sen; 
				smapsen = smapsen + ncases;
				gspe = gspe + spe; 
				smapspe = smapspe + ncontrol;
				acc = acc/testSample.n_rows;
				sen = sen/ncases;
				spe = spe/ncontrol;

				trainRoc(lastInserted) = auul;
				bcoef.row(lastInserted) = coef.t();	
				accuracy(lastInserted) = acc;
				sensitivity(lastInserted) = sen;
				specificity(lastInserted) = spe;
				taccuracy(lastInserted)=tracc;
				tsensitivity(lastInserted)=trsen; 
				tspecificity(lastInserted)=trspe;
				zNRI.row(lastInserted) = modelReclas.tz_NRIs.t();
				zIDI.row(lastInserted) = modelReclas.tz_IDIs.t();
				test_zNRI.row(lastInserted) = modelReclas.z_NRIs.t();
				test_zIDI.row(lastInserted) = modelReclas.z_IDIs.t();
				NRI.row(lastInserted) = modelReclas.NRIs.t();
				IDI.row(lastInserted) = modelReclas.IDIs.t();
				lastInserted=lastInserted+1;
}
			}
		}
		if ((lastInserted>0)&&(lastInserted<loops))
		{
			trainRoc.resize(lastInserted);
			accuracy.resize(lastInserted);
			sensitivity.resize(lastInserted);
			specificity.resize(lastInserted);
			taccuracy.resize(lastInserted);
			tsensitivity.resize(lastInserted);
			tspecificity.resize(lastInserted);
			bcoef.resize(lastInserted,n_var);
			NRI.resize(lastInserted,n_var2);
			IDI.resize(lastInserted,n_var2);
			zNRI.resize(lastInserted,n_var2);
			zIDI.resize(lastInserted,n_var2);
			test_zNRI.resize(lastInserted,n_var2);
			test_zIDI.resize(lastInserted,n_var2);
			
			if ((3*lastInserted)<loops) Rcout<<"Warning: only "<<lastInserted<<" samples of "<<loops<<" were inserted\n";
		}
		double blindAUC=0.5;
		double BlindAccuracy=0.5;
		double BlindSensitivity=0;
		double BlindSpecificity=0;
		if (lastInserted>0)
		{
			blindAUC = rocAUC(testoutcome,testprediction,"auto");
			BlindAccuracy = gacc/static_cast<double>(smaptot);
			BlindSensitivity = gsen/static_cast<double>(smapsen);
			BlindSpecificity = gspe/static_cast<double>(smapspe);
		}
//		Rcout<<"Blind AUC: "<<blindAUC<<"\n";
		Rcpp::List resul =Rcpp::List::create(
			Rcpp::Named("test_zNRI")=Rcpp::wrap(test_zNRI), 
			Rcpp::Named("test_zIDI")=Rcpp::wrap(test_zIDI),
			Rcpp::Named("rawbetas")=Rcpp::wrap(rawbetas), 
			Rcpp::Named("blindAUC")=Rcpp::wrap(blindAUC) 
			);
		Rcpp::List result = Rcpp::List::create(Rcpp::Named("BlindAccuracy")=Rcpp::wrap(BlindAccuracy), 
						Rcpp::Named("BlindSensitivity")=Rcpp::wrap(BlindSensitivity), 
						Rcpp::Named("BlindSpecificity")=Rcpp::wrap(BlindSpecificity), 
						Rcpp::Named("bcoef")=Rcpp::wrap(bcoef), 
						Rcpp::Named("zNRI")=Rcpp::wrap(zNRI), 
						Rcpp::Named("zIDI")=Rcpp::wrap(zIDI), 
						Rcpp::Named("trainRoc")=Rcpp::wrap(trainRoc), 
						Rcpp::Named("NRI")=Rcpp::wrap(NRI), 
						Rcpp::Named("IDI")=Rcpp::wrap(IDI), 
						Rcpp::Named("taccuracy")=Rcpp::wrap(taccuracy), 
						Rcpp::Named("tsensitivity")=Rcpp::wrap(tsensitivity), 
						Rcpp::Named("tspecificity")=Rcpp::wrap(tspecificity), 
						Rcpp::Named("accuracy")=Rcpp::wrap(accuracy), 
						Rcpp::Named("sensitivity")=Rcpp::wrap(sensitivity), 
						Rcpp::Named("specificity")=Rcpp::wrap(specificity), 
						Rcpp::Named("testprediction")=Rcpp::wrap(testprediction), 
						Rcpp::Named("testoutcome")=Rcpp::wrap(testoutcome),
						Rcpp::Named("resul")=resul
						); 
						
		return (result);
	} 
	catch( std::exception &ex ) 
	{ 
		forward_exception_to_r( ex );
		return 0;
	}	
	catch(...) 
	{ 
		::Rf_error( "c++ exception (unknown reason)" );
		return 0;
	}
}

std::string binaryFowardSelection(const unsigned int size,const int totSamples,const double zthr2,const std::string &selType,
								const int loops,const std::string &Outcome,const std::string &timeOutcome,const std::string &type,const unsigned int maxTrainModelSize,
								const mat &casesample,const mat &controlsample,const std::vector < std::string > &ovnames,
								std::map<std::string, int> &lockup,const std::string &baseForm,std::vector <int> &mynamesLoc, 
								const std::vector<std::string> &covari, vec &basevalues,bool isRandom,unsigned int featuresize, vec &zadjst)
{
//	arma_rng::set_seed_random();
	unsigned int inserted = 0;
	double zmin = 0;
	double maxrec = 0;
	arma::mat mysample;
	arma::mat myTestsample;
	arma::mat mysampleMat;
	arma::mat myTestsampleMat;
	arma::mat mysampleMatbase;
	arma::mat myTestsampleMatbase;
	std::string frm1 = baseForm;
	vec bestmodel;
	vec bestmodelCoef;
	vec newmodel;
	vec bestpredict_train;
	vec bestpredict;
	vec singleTrainPredict;
	vec singleTestPredict;
	std::vector < std::string > vnames = ovnames;
	std::string gfrm1;

  	const int z_idi=0;
  	const int z_nri=1;
	int sizecases = casesample.n_rows;
	int sizecontrol = controlsample.n_rows;
	int trainsize = totSamples*2;
	int outcomeidx =  lockup[Outcome];

	
	if (loops > 1)
	{
		uvec ntestcases;
		uvec ntestcontrol;
		unsigned int minnsamples = 2;
		mat myTestCases;
		mat myTestControl;
		vec auxcasesample=zeros<vec>(sizecases);
		vec auxcontrolsample=zeros<vec>(sizecontrol);
		uvec samCases;
		uvec samControl;
		do
		{
			auxcasesample=zeros<vec>(sizecases);
			auxcontrolsample=zeros<vec>(sizecontrol);
			samCases = randi<uvec>(totSamples, distr_param(0,sizecases-1));
			samControl = randi<uvec>(totSamples, distr_param(0,sizecontrol-1));
			for (int i =0;i<totSamples;i++)
			{
				auxcasesample[(samCases(i))]=1;
				auxcontrolsample[(samControl(i))]=1;
			}				
			ntestcases = find(auxcasesample==0);
			ntestcontrol = find(auxcontrolsample==0);
		}
		while ((ntestcases.n_elem<minnsamples) || (ntestcontrol.n_elem<minnsamples));

		mysample=join_cols(casesample.rows(samCases),controlsample.rows(samControl));
		myTestsample=join_cols(casesample.rows(ntestcases),controlsample.rows(ntestcontrol));
		vec outcome = mysample.col(outcomeidx);
		vec testoutcome = myTestsample.col(outcomeidx);
		if (isRandom) 
		{
			if ((size+1)<casesample.n_cols) // if true will randomy sample columns
			{
				vec outcome = mysample.col(outcomeidx);
				arma::mat randcases;
				arma::mat randcontrol;
				uvec rindex= randi<uvec>(1,distr_param(1,casesample.n_cols-1));
				if (size>1)
				{
					rindex = sort_index(randu<vec>(casesample.n_cols-1))+ones<uvec>(casesample.n_cols-1);
					rindex.resize(size);
				}
				randcases = casesample.cols(rindex);
				randcontrol = controlsample.cols(rindex);
				mysample=join_cols(randcases.rows(samCases),randcontrol.rows(samControl));
				mysample.insert_cols(outcomeidx,1);
				myTestsample=join_cols(randcases.rows(ntestcases),randcontrol.rows(ntestcontrol));
				myTestsample.insert_cols(outcomeidx,1);
//				Rcout << endl << "Co(" << myTestsample.n_rows << "," << myTestsample.n_cols << ")" << endl;
			}
			mysample.col(outcomeidx) = outcome(sort_index(randu(mysample.n_rows)));
			myTestsample.col(outcomeidx) = testoutcome(sort_index(randu(myTestsample.n_rows)));
		}
		else
		{
			myTestCases=casesample.rows(ntestcases);
			myTestControl=controlsample.rows(ntestcontrol);
			int ncases = myTestCases.n_rows;
			int ncontrol = myTestControl.n_rows;
			int ntesmin = std::max(ncases,ncontrol);
			int randnum = (randi<uvec>(1))[0];

			myTestsample = zeros<mat>(2*ntesmin,myTestCases.n_cols);
			for (int i=0;i<ntesmin;i++)
			{
				myTestsample.row(2*i) =  myTestCases.row((i + randnum) % ncases);
				myTestsample.row(2*i+1) = myTestControl.row((i + randnum) % ncontrol);
			}
		}
	}
	else
	{
		mysample=join_cols(casesample,controlsample);
  	  	myTestsample=join_cols(casesample,controlsample);
	}

	
    mat myoutcome(mysample.n_rows,2);
	myoutcome.col(1)=mysample.col(outcomeidx);
       
	mysampleMatbase.set_size(mysample.n_rows);
	myTestsampleMatbase.set_size(myTestsample.n_rows);
	mysampleMatbase.fill(1.0);
	myTestsampleMatbase.fill(1.0);
	if (covari[0]!="1")
	{
		if (type == "COX")
		{
			myoutcome.col(0)=mysample.col(lockup[timeOutcome]);
		}
		for (unsigned int i=0;i<covari.size();i++)
		{
			mysampleMatbase=join_rows(mysampleMatbase,mysample.col(lockup[covari[i]]));
			myTestsampleMatbase=join_rows(myTestsampleMatbase,myTestsample.col(lockup[covari[i]]));		  
		}
	}
	else
	{
		if (type == "COX")
		{
			myoutcome.col(0)=mysample.col(lockup[timeOutcome]);
		}
	}
    if((mysampleMatbase.n_cols>1)or(type=="COX"))
	{
        bestmodel =modelFittingFunc(myoutcome,mysampleMatbase,type);
	}
    else
    {
       bestmodel=mean(myoutcome.col(1));
       if (type=="LOGIT") 
	   {
		   if ((bestmodel(0) >= DOUBLEEPS)&&(bestmodel(0) < 1))  
		   {
			   bestmodel=log(bestmodel/(1.0-bestmodel));
		   }
		   else 
		   {
			   if (bestmodel(0) == 1 ) bestmodel(0) = THRESH;
			   else bestmodel(0) = MTHRESH;
		   }
	   }
	}
  	bestpredict_train=predictForFresaFunc(bestmodel,mysampleMatbase,"prob",type);
  	if (loops > 1) bestpredict=predictForFresaFunc(bestmodel,myTestsampleMatbase,"prob",type);
	int changes = 1;
	int jmax = -1;
  	unsigned int j;
	std::string inname = "Inserted";
	vec testoutvec = myTestsample.col(outcomeidx);
	vec trainoutvec = mysample.col(outcomeidx);
	double merror=1.0;
	double tol=1.0e-6; // error tolerance
	unsigned int maximumNoselectedCount=25;
	int cycle=0;
	vec iprob(6);
	vec iprob_t(6);
	uvec sindex=linspace<uvec>(0,size-1,size);
	uvec nideex=sindex;
	for ( j=0; j<size;j++) 
	{
		if (vnames[j] != inname)
		{
			nideex[j]=lockup[vnames[j]];
		}
	}
	while ((changes>0)&&(merror>tol))
	{
		changes = 0;
	  	maxrec = zthr2;
	  	jmax = -1;
		unsigned int lastpass=0;
		if ((cycle>0)&&(size > (maximumNoselectedCount/5)))
		{
			vec res = bestpredict_train-trainoutvec;
			vec correlations=zeros<vec>(size);
			for ( j=0; j<size;j++)
			{
				if (vnames[j] != inname)
				{
					vec cr = abs(cor(res ,mysample.col(nideex[j])));
					if (!cr.has_nan())
					{
						correlations[j]=cr[0];
					}
				}
			}
			sindex = sort_index(correlations,"descend");
		}
		unsigned int passtestcout=1;
		double ztmin=0;
		int nfors = 0;
		unsigned int  jn=0;
		unsigned int  ajnskip=0;
		unsigned int noselectstop = maximumNoselectedCount;
	  	for (jn=0; jn<size;jn++,nfors++)
	  	{
			j=sindex[jn];
			zmin=0.0;
			unsigned int jnskip = 0;
			if (jn > maximumNoselectedCount)
			{
				jnskip = randi<uvec>(1,distr_param(0,jn/maximumNoselectedCount))[0];
			}
	  		if (vnames[j] != inname)
	  		{
				unsigned int vidx = nideex[j];
				mysampleMat=join_rows(mysampleMatbase,mysample.col(vidx));
				if (loops > 1) myTestsampleMat=join_rows(myTestsampleMatbase,myTestsample.col(vidx));
	  			gfrm1 = frm1+" + "+vnames[j];
             	if((mysampleMat.n_cols>1)or(type=="COX"))
				{
					newmodel =modelFittingFunc(myoutcome,mysampleMat,type);
				}
				else
				{
					newmodel=mean(myoutcome.col(1));
					if (type=="LOGIT") 	
					{
						if ((newmodel[0] >= DOUBLEEPS)&&(newmodel[0] < 1))  
						{
							newmodel=log(newmodel/(1.0-newmodel));
						}
						else 
						{
						   if (newmodel[0] == 1) newmodel[0] = THRESH;
						   else newmodel[0] = MTHRESH;
						}
					}
				}
				if (!newmodel.has_nan()) 
	  			{
  					singleTrainPredict=predictForFresaFunc(newmodel,mysampleMat,"prob",type);
					iprob_t = improveProbFunc(bestpredict_train,singleTrainPredict,trainoutvec);            
					merror = mean(abs(singleTrainPredict-trainoutvec)); // train estimation error
					if (selType=="zIDI") 
					{
						ztmin=iprob_t[z_idi];
					}
					else
					{
						ztmin=iprob_t[z_nri];
					}
					if (loops > 1) 
					{
						singleTestPredict=predictForFresaFunc(newmodel,myTestsampleMat,"prob",type);
						iprob = improveProbFuncSamples(bestpredict,singleTestPredict,testoutvec,trainsize,iprob_t[5],iprob_t[4]);
						if (!iprob.has_nan() && !iprob_t.has_nan())
						{
							if (selType=="zIDI") 
							{
								zmin = std::min(iprob[z_idi],iprob_t[z_idi]);
							}
							else
							{
								zmin = std::min(iprob[z_nri],iprob_t[z_nri]);
							}
						}
					}
					else
					{
						zmin=ztmin;
					}
//					if ((zmin+ztmin)  > 2*zthr2)
					if (ztmin  > zthr2)
					{
						passtestcout += 1+ajnskip;
						lastpass = 0;
					}
					else
					{
						lastpass += 1+ajnskip;
					}
					if ( zmin > maxrec )
					{
						jmax = j;
						maxrec = zmin;
						if (cycle>0) noselectstop = maximumNoselectedCount/5;
					}
					if (cycle==0) 
					{
						basevalues[j] = ztmin;
					}
					if ((loops>1)||(cycle>0))
					{
						if ((merror<tol) || (lastpass>noselectstop))
						{
							jn=size; // lets stop.
						}
						jn += jnskip;
					}
				}
				else
				{
					lastpass += 1+ajnskip;
				}
				ajnskip=jnskip;
			}
	  	}
		double a_padj=0;
		double porg = R::pnorm(-zthr2, 0.0, 1.0, 1, 0);
		a_padj = 4.0*((double)(passtestcout)/(featuresize-cycle));
//		if (cycle == 0) a_padj = 2.0*((double)(passtestcout)/(featuresize-cycle)); // to avoid selecting a false feature from the beginning 

		if (a_padj>1) a_padj = 1.0;
		a_padj = a_padj*porg;
		if (!isRandom) zadjst[cycle] = -1.0*(qnorm(a_padj,0.0,1.0)); 
		else zadjst[cycle] = zthr2; 
		
//		Rcout  <<size<<" : "<< cycle<<" : "<< nfors <<" : "<<jn<<" : "<< zadjst[cycle]<< " : "<< featuresize << " : "<< passtestcout << " : "<< a_padj << " : " << lastpass<<"\n";	
		
		if ((jmax >= 0) && (maxrec > zadjst[cycle]) && (inserted<maxTrainModelSize))   
	  	{
			unsigned int vidx=nideex[jmax];
			
			if (vnames[jmax] != inname)
			{
				gfrm1 = frm1+" + ";
				gfrm1 = gfrm1+ovnames[jmax];
				mysampleMatbase=join_rows(mysampleMatbase,mysample.col(vidx));

				bestmodel = modelFittingFunc(myoutcome,mysampleMatbase,type);
				bestpredict_train = predictForFresaFunc(bestmodel,mysampleMatbase,"prob",type);

				if (loops > 1) 
				{
					myTestsampleMatbase=join_rows(myTestsampleMatbase,myTestsample.col(vidx));
					bestpredict = predictForFresaFunc(bestmodel,myTestsampleMatbase,"prob",type);
				}
				merror = mean(abs(bestpredict_train-trainoutvec)); // train estimation error
				frm1 = gfrm1;	
				mynamesLoc.push_back(jmax);
				inserted = inserted + 1;
				if (merror>tol) 
				{
					changes = changes + 1;
				}
			}
			vnames[jmax] = inname;
	  	}
		++cycle;  		
	}
	return frm1;
}

extern "C" SEXP ReclassificationFRESAModelCpp(SEXP _size, SEXP _fraction,SEXP _pvalue, SEXP _loops,  SEXP _covariates, SEXP _Outcome, SEXP _variableList, SEXP _maxTrainModelSize, SEXP _type, SEXP _timeOutcome, SEXP _selType, SEXP _dataframe,SEXP _colnames,SEXP _featuresize,SEXP _cores)
{
try {   // R_CStackLimit=(uintptr_t)-1;
 		#ifdef _OPENMP
//		Rcout<<"::============ Available CPU Cores(Threads) : "<<omp_get_num_procs()<<" ===Requested Threads : "<< as<unsigned int>(_cores) << endl;
	    omp_set_num_threads(as<unsigned int>(_cores));
//	    omp_set_num_threads(1);
		#endif	
//		srand (time(NULL));
//		arma_rng::set_seed_random();
		unsigned int size = as<unsigned int>(_size); 
		double fraction = as<double>(_fraction);
		double pvalue = as<double>(_pvalue);
		int loops = as<int>(_loops);
		std::string Outcome = as<std::string>(_Outcome); 
		CharacterVector colnames(_colnames);
	    unsigned int maxTrainModelSize = Rcpp::as<unsigned int>(_maxTrainModelSize);
		std::string type = as<std::string>(_type);
		std::string timeOutcome = as<std::string>(_timeOutcome);
		std::string selType = as<std::string>(_selType);
		Rcpp::NumericMatrix DF(_dataframe);
    	arma::mat dataframe(DF.begin(), DF.rows(), DF.cols(), false);
		std::vector<std::string>covari=as<std::vector<std::string> >(_covariates);
		std::string baseForm=Outcome;
		std::map<std::string, int> lockup;
		std::vector<std::string> formulas;
		std::vector<int> mynames;
		std::vector<std::string> ovnames=as<std::vector<std::string> >(CharacterVector(_variableList));
		unsigned int featuresize=Rcpp::as<unsigned int>(_featuresize);
		vec bbasevalues=zeros<vec>(size);
		arma::mat AZthr=zeros<arma::mat>(size,loops);

		double zthr = -1.0*(qnorm(pvalue,0.0,1.0)); //double zthr = abs(Rf_qnorm5(pvalue,0.0, 1.0, 1, 0));
//		FILE *arch=fopen("Jose.txt","a+");

		int len_outcome=Outcome.size();
		bool israndom=false;
		if (len_outcome>6)
		{
			std::string key=Outcome;
			key.resize(6);
			if (key=="RANDOM")
			{
				israndom=true;
				Outcome=Outcome.substr(6,len_outcome-6);
				if (size<ovnames.size())
				{
					for (unsigned int i=1; i<=size; i++) 
					{
						ovnames[i-1]=colnames[i];
					}
				}
			}
		}




		for(int i=0;i<CharacterVector(colnames).size();i++)
		{
			lockup[std::string(CharacterVector(colnames)[i])]=i;
		}
		if (type=="COX")
		{
		  baseForm = "Surv("+timeOutcome+","+Outcome+")";
		}
		baseForm = baseForm+" ~ ";
		for(unsigned int i=0;i<covari.size();i++)
		{
			baseForm = baseForm+"+"+covari[i];
		}
		

		int indoutcome=lockup[Outcome];
		if (israndom)
		{
			vec outcome = dataframe.col(indoutcome);
			dataframe.col(indoutcome)=outcome(sort_index(randu(dataframe.n_rows)));
		}
		mat casesample = dataframe.rows(find(dataframe.col(indoutcome)==1));
		mat controlsample = dataframe.rows(find(dataframe.col(indoutcome)==0));
		int sizecases   = casesample.n_rows;
		int sizecontrol = controlsample.n_rows;
		int minsize = std::min(sizecases,sizecontrol);
		int totSamples = (int)(fraction*minsize+0.49);

		double zthr2 = zthr;
		if (fraction<1) 
		{
			zthr2 = zthr*std::sqrt(fraction);
		}
		if (zthr2<std::abs(qnorm(0.1,0,1))) zthr2 = std::abs(qnorm(0.1,0,1)); 

//		Rcout << "\n" <<pvalue<<" <-pv : z-> "<<zthr2<<" Cases: "<<  sizecases <<" Controls: "<<  sizecontrol << "\n";
//		fprintf(arch,"Size=%d, Pvalue=%8.3f,Zthr=%8.3f, Cases= %d, Control= %d, Total=%d (%d,%d), f=%d\n",size,pvalue,zthr2,sizecases,sizecontrol,totSamples,casesample.n_cols,controlsample.n_cols,featuresize);
//		uvec samCases = randi<uvec>(totSamples, distr_param(0,sizecases-1));
//		fprintf(arch,"%d,%d,%d,%d,%d,%d \n",samCases[0],samCases[1],samCases[2],samCases[3],samCases[4],samCases[5]);

		if (size > ovnames.size()) size = ovnames.size();
		std::string inname = "Inserted";
		for(unsigned int j=0;j<covari.size();j++)
		    for(unsigned int i=0;i<ovnames.size();i++)
					 if (covari[j]==ovnames[i]) ovnames[i]=inname;  	
#pragma omp parallel for schedule(dynamic) ordered shared (mynames,formulas,bbasevalues,AZthr)  
    	for (int doOver=0; doOver<loops;doOver++)
	 	{
			vec zadjst=zeros<vec>(size);
			zadjst.fill(100);
			vec basevalues=zeros<vec>(size);
//			std::vector<std::string> oovnames=ovnames;
			std::vector <int> mynamesLoc;
//			int doBootstrapSample = 1+loops*((doOver % 20)!=0);

			std::string  frm = binaryFowardSelection(size,totSamples,zthr2,selType,loops,
			Outcome,timeOutcome,type,maxTrainModelSize,casesample,controlsample,
			ovnames,lockup,baseForm,mynamesLoc,covari,basevalues,israndom,featuresize,zadjst);

#pragma omp critical
{
				mynames.insert(mynames.end(), mynamesLoc.begin(), mynamesLoc.end());
				formulas.push_back(frm);
				for (unsigned int i=0;i<size;i++) 
				{
					bbasevalues[i] += basevalues[i];
					AZthr(i,doOver) = zadjst[i];
				}
//				fprintf(arch,"doOver %10d: %s\n",doOver,frm.c_str());
//				fflush(arch);

}

			if ((doOver%100)==99) 
				{
#pragma omp critical
					Rcout << ".";	
//					Rcout << doOver <<" : "<< frm << std::endl;	
				}

  		}
		vec medAZthr=zeros<vec>(size);
		for (unsigned int i=0;i<size;i++) 
		{
			bbasevalues[i] =  bbasevalues[i]/loops;
			rowvec drow = AZthr.row(i);
			uvec whogz = find(drow < 100);
			if (whogz.n_elem > 0) medAZthr[i] = median(drow(whogz));
			else
			{
				if (i>0) medAZthr[i] = medAZthr[i-1];
			}
		}


  		if (mynames.size() == 0) mynames.push_back(0);
//		fprintf(arch,"%d Var# %d\n",loops,mynames.size());
//		fprintf(arch,"Ended %d\n",formulas.size());
//		fclose(arch);
	    List result = List::create(Named("mynames")=wrap(mynames),Named("formula.list")=wrap(formulas),Named("Base.values")=wrap(bbasevalues),Named("Zthr")=wrap(medAZthr));
	    return (result);
	} 
	catch( std::exception &ex ) 
	{ 
		forward_exception_to_r( ex );
		return 0;
	}	
	catch(...) 
	{ 
		::Rf_error( "c++ exception (unknown reason)" );
		return 0;
	}
}
