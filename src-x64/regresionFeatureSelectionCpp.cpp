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


struct NeRIbootVal{
	mat bootmodel;
	mat bootmodelmeans;
	mat NeRi;
	mat tpvalue;
	mat wpvalue;
	mat spvalue;
	mat Fpvalue;
	mat test_tpvalue;
	mat test_wpvalue;
	mat test_spvalue;
	mat test_Fpvalue;
 }NeRIredBoot;

int redCurmodel_S_lastRemovedNeRI;


//********************************************************************************************
//**==========================================================================================
//********************************************************************************************
extern "C" SEXP bootstrapValidationResCpp(SEXP _fraction,SEXP _loops,SEXP _dataframe,SEXP _type,SEXP _response,SEXP _bestdataframe)
{
	try 
	{  //R_CStackLimit=(uintptr_t)-1;
//		omp_set_num_threads(4);
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
        if (type == "COX")
		{
			dataframe=dat.cols(1,dat.n_cols-1);
			bestXframe=bestdatX.cols(1,bestdatX.n_cols-1);
		}
		else 
		{
			dataframe=dat;
			bestXframe=bestdatX;
		}
	    mat ttcasesample = join_rows(Response,dataframe);
	    mat tbestSample = join_rows(Response,bestXframe);

		
		int tsizecases  = dataframe.n_rows;
		std::vector<double>trainResiduals;
		vec trainSampledRMSE(loops);
		std::vector<double>testOutcome;
		std::vector<double>testPrediction;
		std::vector<double>testResiduals;
		vec testSampledRMSE(loops);
		vec varoutcome(loops);
		mat bcoef(loops,dataframe.n_cols); 	  
		mat bmeans(loops,dataframe.n_cols); 
		mat	NeRi(loops,n_var2); 
		mat	tpvalue=ones<mat>(loops,n_var2);  
		mat	wpvalue=ones<mat>(loops,n_var2);  
		mat	spvalue=ones<mat>(loops,n_var2);  
		mat	Fpvalue=ones<mat>(loops,n_var2);  
		mat	test_tpvalue=ones<mat>(loops,n_var2); 
		mat	test_wpvalue=ones<mat>(loops,n_var2); 
		mat	test_spvalue=ones<mat>(loops,n_var2); 
		mat	test_Fpvalue=ones<mat>(loops,n_var2);
		int lastInserted=0;		
		int totSamples = (int)(fraction*tsizecases);
//		uvec eqindex = equSamples(tcasesample);
//		tcasesample = tcasesample.rows(eqindex);
//		bestSample = bestSample.rows(eqindex);
//		int sizecases = ttcasesample.n_rows;
//		double omin=min(Outcome);
//		double range = max(Outcome)-omin;

#pragma omp parallel for schedule(dynamic) ordered shared(lastInserted,trainResiduals,testOutcome,testPrediction,testResiduals,trainSampledRMSE,testSampledRMSE,bcoef,bmeans,NeRi,tpvalue,wpvalue,spvalue,Fpvalue,test_tpvalue,test_wpvalue,test_spvalue,test_Fpvalue,varoutcome)
		for (int doOver=0;doOver<loops;doOver++)
		{ 
			uvec eqindex = equSamples(ttcasesample);
			mat tcasesample = ttcasesample.rows(eqindex);
			mat bestSample = tbestSample.rows(eqindex);
			int sizecases = tcasesample.n_rows;
			mat trainingSample;
			mat bestTrainingSample;
			uvec ntestSample;
			uvec samCases;
			unsigned int cycle = 0;
			unsigned int mintestSample = 2;
			do
			{ 
				samCases = randi<uvec>(totSamples, distr_param(0,sizecases-1));
				vec auxcasesample=zeros<vec>(sizecases);
				for (int i=0;i<totSamples;i++)
				{
					auxcasesample[samCases[i]]=1;
				}
				ntestSample=find(auxcasesample==0);
				++cycle;
				if (cycle > 5) 
				{
					ntestSample=find(auxcasesample==1); //same as train if no enough data sample
				}
			} while (ntestSample.n_elem<mintestSample);	

			trainingSample = tcasesample.rows(samCases);
			bestTrainingSample = bestSample.rows(samCases);	
			varoutcome(doOver) = var(tcasesample.col(1));
			vec trainmodel = modelFittingFunc(trainingSample.cols(0,1),trainingSample.cols(2,trainingSample.n_cols-1),type);
			if (!trainmodel.has_nan())
			{
				mat testSample = tcasesample.rows(ntestSample);
				mat bestTestSample = bestSample.rows(ntestSample);
				vec coef = trainmodel;
				vec coxmean;
				if (type == "COX") 
				{
					 coef = trainmodel.subvec(0,(trainmodel.n_elem/2)-1);
					 coxmean = trainmodel.subvec((trainmodel.n_elem/2),trainmodel.n_elem-1);
				}

				vec residualTrain = residualForFRESAFunc(trainmodel,trainingSample.cols(2,trainingSample.n_cols-1),"",type,trainingSample.cols(0,1));
				vec residualTest = residualForFRESAFunc(trainmodel,testSample.cols(2,testSample.n_cols-1),"",type,testSample.cols(0,1));
				gvarNeRI modelReclas = getVarResFunc(trainingSample,type,testSample,totSamples,bestTrainingSample,bestTestSample);

#pragma omp critical
{
				for (unsigned int i=0;i<residualTrain.n_elem;i++)
				{
					trainResiduals.push_back(residualTrain(i));
				}

				for (unsigned int i=0;i<residualTest.n_elem;i++)
				{
					testOutcome.push_back(testSample(i,1));
					testPrediction.push_back(testSample(i,1)+residualTest(i));
					testResiduals.push_back(residualTest(i));
				}
				trainSampledRMSE(lastInserted) = std::sqrt(mean(residualTrain%residualTrain));
				testSampledRMSE(lastInserted) = std::sqrt(mean(residualTest%residualTest));
				bcoef.row(lastInserted) = coef.t();	
				if (type == "COX") bmeans.row(lastInserted) = coxmean.t();
				NeRi.row(lastInserted) = modelReclas.NeRIs.t();
				tpvalue.row(lastInserted) = modelReclas.tP_value.t();
				wpvalue.row(lastInserted) = modelReclas.WilcoxP_value.t();
				spvalue.row(lastInserted) = modelReclas.BinP_value.t();
				Fpvalue.row(lastInserted) = modelReclas.FP_value.t();
				test_tpvalue.row(lastInserted) = modelReclas.testData_tP_value.t();
				test_wpvalue.row(lastInserted) = modelReclas.testData_WilcoxP_value.t();
				test_spvalue.row(lastInserted) = modelReclas.testData_BinP_value.t();
				test_Fpvalue.row(lastInserted) = modelReclas.testData_FP_value.t();
				lastInserted = lastInserted + 1;
}
			}
		}
		if ((lastInserted>0)&&(lastInserted<loops))
		{
			trainSampledRMSE.resize(lastInserted);
			testSampledRMSE.resize(lastInserted);
			bcoef.resize(lastInserted,dataframe.n_cols);	
			if (type == "COX") bmeans.resize(lastInserted,dataframe.n_cols);
			NeRi.resize(lastInserted,n_var2);
			test_tpvalue.resize(lastInserted,n_var2);
			test_wpvalue.resize(lastInserted,n_var2);
			test_spvalue.resize(lastInserted,n_var2);
			test_Fpvalue.resize(lastInserted,n_var2);
			tpvalue.resize(lastInserted,n_var2);
			wpvalue.resize(lastInserted,n_var2);
			spvalue.resize(lastInserted,n_var2);
			Fpvalue.resize(lastInserted,n_var2);
		
			if ((3*lastInserted)<loops) Rcout<<"Warning: only "<<lastInserted<<" samples of "<<loops<<" were inserted\n";
		}
		double totsd =0;
		for (int doOver=0;doOver<loops;doOver++) totsd += varoutcome(doOver);
		totsd = std::sqrt(totsd/loops);
		
//		Rcout<<" First t-pvalue: "<<tpvalue(0,0)<<"\n";

	 	rowvec bootmodel = mean (bcoef);
		if (type == "COX") 
		{ 
			rowvec CoxMmodel = mean (bmeans);
			rowvec aux(bootmodel.n_elem*2);
			for (unsigned int i=0;i<bootmodel.n_elem;i++)
			{
				aux[i] = bootmodel[i];
				aux[bootmodel.n_elem+i] = CoxMmodel[i];
			}
			bootmodel.resize(aux.n_elem);
			bootmodel = aux;
		}
	 	vec basemodel = modelFittingFunc(ttcasesample.cols(0,1),ttcasesample.cols(2,ttcasesample.n_cols-1),type);
		double startRMSE = std::sqrt(mean(square(residualForFRESAFunc(basemodel,ttcasesample.cols(2,ttcasesample.n_cols-1),"",type,ttcasesample.cols(0,1)))));
		double bootRMSE = std::sqrt(mean(square(residualForFRESAFunc(bootmodel.t(),ttcasesample.cols(2,ttcasesample.n_cols-1),"",type,ttcasesample.cols(0,1)))));
		Rcpp::List result = Rcpp::List::create(Rcpp::Named("NeRi")=Rcpp::wrap(NeRi), \
						Rcpp::Named("bcoef")=Rcpp::wrap(bcoef), \
						Rcpp::Named("tpvalue")=Rcpp::wrap(tpvalue), \
						Rcpp::Named("wpvalue")=Rcpp::wrap(wpvalue), \
						Rcpp::Named("spvalue")=Rcpp::wrap(spvalue), \
						Rcpp::Named("Fpvalue")=Rcpp::wrap(Fpvalue), \
						Rcpp::Named("test_tpvalue")=Rcpp::wrap(test_tpvalue), \
						Rcpp::Named("test_wpvalue")=Rcpp::wrap(test_wpvalue), \
						Rcpp::Named("test_spvalue")=Rcpp::wrap(test_spvalue), \
						Rcpp::Named("test_Fpvalue")=Rcpp::wrap(test_Fpvalue), \
						Rcpp::Named("trainSampledRMSE")=Rcpp::wrap(trainSampledRMSE), \
						Rcpp::Named("testOutcome")=Rcpp::wrap(testOutcome), \
						Rcpp::Named("testPrediction")=Rcpp::wrap(testPrediction), \
						Rcpp::Named("testResiduals")=Rcpp::wrap(testResiduals), \
						Rcpp::Named("testSampledRMSE")=Rcpp::wrap(testSampledRMSE),\
						Rcpp::Named("startRMSE")=Rcpp::wrap(startRMSE), \
						Rcpp::Named("bootRMSE")=Rcpp::wrap(bootRMSE), \
						Rcpp::Named("trainResiduals")=Rcpp::wrap(trainResiduals), \
						Rcpp::Named("OutcomeSD")=Rcpp::wrap(totsd)
						
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

std::string residualFowardSelection(const unsigned int size,const int totSamples,const double pthr2,const std::string &testType,
									const int loops,const std::string &Outcome,const std::string &timeOutcome,const std::string &type,const unsigned int maxTrainModelSize,
									const mat &dataframe, const std::vector < std::string > &ovnames,
									std::map<std::string, int> &lockup,const std::string &baseForm,std::vector <int> &mynamesLoc, 
									const std::vector<std::string> &covari,vec &basevalues,bool isRandom,unsigned int featuresize, vec &pthrs)
{
	std::string frm1 = baseForm;
	unsigned int inserted = 0;
	double iprob=0;
	unsigned int sizecases = dataframe.n_rows;
//	unsigned int featuresize=std::max((unsigned int)dataframe.n_cols-1,(unsigned int)ovnames.size());

//	Rcout <<" Hello \n";	

	arma::mat mysample = dataframe;
	arma::mat myTestsample = dataframe;

	arma::mat mysampleMatbase;
	arma::mat myTestsampleMatbase;
	vec bestmodel;
	vec newmodel;
	vec bestTrainResiduals;
	vec bestResiduals;
	std::vector < std::string > vnames = ovnames;
	std::string inname = "Inserted";

	int outinx = lockup[Outcome];
	vec randOutput=dataframe.col(outinx);
	double 	istdoutcome=stddev(randOutput);
	if (istdoutcome>0) istdoutcome=1.0/istdoutcome;
	if (isRandom)
	{
		randOutput =  randOutput(sort_index(randu<vec>(dataframe.n_rows)));
	}
	
	int timeinx	= 0;
	if (type=="COX")
	{
		timeinx	= lockup[timeOutcome];
	}

	for(unsigned  int j=0;j<covari.size();j++)
		for(unsigned int i=0;i<vnames.size();i++)
			if (covari[j]==vnames[i]) vnames[i]=inname;
	unsigned int cycle = 0;
	if (loops > 1)
	{
		unsigned int mitestsize=2;
		uvec notinserted;
		uvec samCases;
		vec auxcasesample;
		double omin=min(randOutput);
		double range=max(randOutput)-omin;
		cycle = 0;
		do
		{
			samCases= randi<uvec>(totSamples, distr_param(0,sizecases-1));
			auxcasesample=zeros<vec>(sizecases);
			for (int i = 0;i<totSamples;i++)
			{
				auxcasesample[samCases(i)]=1;
			}
			notinserted = find(auxcasesample==0);
			++cycle;
			if (cycle > 5) return frm1;	 //exit no models
		}	
		while (notinserted.n_elem<mitestsize);
		mysample = dataframe.rows(samCases);
		myTestsample = dataframe.rows(notinserted);
//		Rcout <<sizecases<<":"<<loops<<":"<< myTestsample.n_rows << ":"<< mysample.n_rows <<"\n";	

		if (isRandom) 
		{		
			outinx=0;
			if ((size+1)<dataframe.n_cols) // if true will randomy sample columns
			{
				arma::mat randcols;
				uvec rindex= randi<uvec>(1,distr_param(1,dataframe.n_cols-1));
				if (size>1)
				{
					rindex= sort_index(randu<vec>(dataframe.n_cols-1))+ones<uvec>(dataframe.n_cols-1);
					rindex.resize(size);
				}
				randcols = dataframe.cols(rindex);
				randcols.insert_cols(outinx,1);
				mysample = randcols.rows(samCases);
				myTestsample = randcols.rows(notinserted);
			}
			mysample.col(outinx) = randOutput(samCases);
		}
		if (isRandom) 
		{
			myTestsample.col(outinx) = randOutput(notinserted);
		}
		else
		{
			myTestsample = myTestsample.rows(equSamples(myTestsample,outinx,5,omin,range));
			mysample = mysample.rows(equSamples(mysample,outinx,5,omin,range));
		}
	}

	mat train_outcome = join_rows(mysample.col(timeinx),mysample.col(outinx));
	mat test_outcome = join_rows(myTestsample.col(timeinx),myTestsample.col(outinx));

//	Rcout <<sizecases<<":"<<loops<<":"<< myTestsample.n_rows << ":"<< mysample.n_rows <<"\n";	

    mat myoutcome(mysample.n_rows,2);
	myoutcome.col(1) = train_outcome.col(1);
	mysampleMatbase.set_size(mysample.n_rows);
	myTestsampleMatbase.set_size(myTestsample.n_rows);
	mysampleMatbase.fill(1.0);
	myTestsampleMatbase.fill(1.0);
	if (type=="COX")
	{
		myoutcome.col(0)=mysample.col(lockup[timeOutcome]);
	}
	if (covari[0]!="1")
	{
		for (unsigned int i=0;i<covari.size();i++)
	  	{
			mysampleMatbase=join_rows(mysampleMatbase,mysample.col(lockup[covari[i]]));
	  		myTestsampleMatbase=join_rows(myTestsampleMatbase,myTestsample.col(lockup[covari[i]]));		  
  	    }
	}
	if((mysampleMatbase.n_cols>1)or(type=="COX"))
		bestmodel = modelFittingFunc(myoutcome,mysampleMatbase,type);
	else
    {
		bestmodel = mean(myoutcome.col(1));
       	if (type=="LOGIT")
		{
			if ((bestmodel(0) >= DOUBLEEPS)&&(bestmodel(0) < 1.0))
			{
				bestmodel = log(bestmodel/(1.0-bestmodel));
			}
			else
			{
				if (bestmodel(0) == 1) bestmodel(0) = THRESH;
				else bestmodel(0) = MTHRESH;
			}
		}
    }
  	bestTrainResiduals=residualForFRESAFunc(bestmodel,mysampleMatbase,"",type,train_outcome);
//	Rcout << "Best:\n";
//	Rcout << bestmodel;
//	Rcout << "....Residuals\n";
//	Rcout << sum(square(bestTrainResiduals)) << "\n";

  	if (loops > 1) bestResiduals=residualForFRESAFunc(bestmodel,myTestsampleMatbase,"",type,test_outcome);
	int changes = 1;
	double minpiri = pthr2;
	int jmax = -1;
  	std::string gfrm1;
	mat myTestsampleMat;
	mat mysampleMat;
	vec testResiduals;
	vec trainResiduals;
	double error=1.0;
	double tol=1.0e-6;
	unsigned int maximumNoselectedCount=25;
	double padjs=1.0;
	uvec sindex=linspace<uvec>(0,size-1,size);
	uvec nideex=sindex;
	for (unsigned int j=0; j<size;j++) 
	{
		if (vnames[j] != inname)
		{
			nideex[j]=lockup[vnames[j]];
		}
	}
	cycle=0;
	while ((changes>0)&&(error>tol)&&(cycle < 2*featuresize))
	{
		changes = 0;
	  	minpiri = pthr2;
	  	jmax = -1;
		unsigned int lastpass=0;
		unsigned int passtestcout=1;
		unsigned int noselectstop = maximumNoselectedCount;
		if ((cycle>0)&&(size > (maximumNoselectedCount/5)))
		{
			vec correlations=zeros<vec>(size);
			for (unsigned int  j=0; j<size;j++)
			{
				if (vnames[j] != inname)
				{
					vec cr = abs(cor(bestTrainResiduals ,mysample.col(nideex[j])));
					if (!cr.has_nan())
					{
						correlations[j]=cr[0];
					}
				}
			}
			sindex = sort_index(correlations,"descend");
		}
	  	for (unsigned int jn=0; jn<size;jn++)
	  	{
			unsigned int j=sindex[jn];
			if (vnames[j] != inname)
	  		{
				unsigned int vidx = nideex[j];

				mysampleMat=join_rows(mysampleMatbase,mysample.col(vidx));
				if (loops > 1) myTestsampleMat=join_rows(myTestsampleMatbase,myTestsample.col(vidx));

	  			gfrm1 = frm1+" + "+vnames[j];
             	if((mysampleMat.n_cols>1)or(type=="COX"))
				{
					newmodel = modelFittingFunc(myoutcome,mysampleMat,type);
				}
			    else
			    {
					newmodel=mean(myoutcome.col(1));
        			if (type=="LOGIT")
					{
						if ((newmodel[0] >= DOUBLEEPS)&&(newmodel[0]<1))
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

					trainResiduals = residualForFRESAFunc(newmodel,mysampleMat,"",type,train_outcome);

//					Rcout << "Model :\n";
//					Rcout << newmodel;
//					Rcout << "....Residuals\n";
//					Rcout << sum(square(trainResiduals)) << "\n";
					

					error = mean(abs(trainResiduals))*istdoutcome;
			   		double mpiri,piri = improvedResidualsFunc(bestTrainResiduals,trainResiduals,testType,totSamples).pvalue;
//					Rcout << piri << " Pvalue \n";
					mpiri = piri;					
					if (loops > 1) 
					{
						testResiduals = residualForFRESAFunc(newmodel,myTestsampleMat,"",type,test_outcome);
						iprob = improvedResidualsFunc(bestResiduals,testResiduals,testType,totSamples).pvalue;
						if (!std::isnan(iprob) && !std::isnan(piri))
						{
							mpiri =std::max(iprob,piri);
						}
					}
//					if ((mpiri+piri) < 2*pthr2)
					if (piri<pthr2)
					{
						passtestcout += 1;
						lastpass = 0;
					}
					else
					{
						lastpass += 1;
					}
					if (!std::isnan(mpiri) && (mpiri  < minpiri))
					{
						jmax = j;
						minpiri = mpiri;
//						Rcout <<cycle<<":"<< error << ":"<< iprob <<":"<< piri <<"\n";	
						if (cycle>0) noselectstop = maximumNoselectedCount/5;
					}
					if (cycle==0) 
					{
						basevalues[j]=piri;
					}
					if ((loops > 1)||(cycle>0))
					{
						if ((error<tol)||(lastpass>noselectstop))
						{
							jn=size; // exit
						}
					}
	  			}
				else
				{
					lastpass += 1;
				}	
	  		}
	  	}
		padjs=1.0;
		padjs = 4.0*((double)(passtestcout)/(featuresize-cycle));
//		if (cycle == 0) padjs = 2.0*((double)(passtestcout)/(featuresize-cycle)); // to avoid selecting a false feature from the beginning 
		if (padjs > 1.0) padjs = 1.0;
		pthrs[cycle] = pthr2*padjs;
		if (pthrs[cycle]>pthr2) pthrs[cycle]=pthr2;
		if (isRandom) pthrs[cycle]=pthr2;
		
		if ((jmax >= 0) && (minpiri < pthrs[cycle]) && (inserted<maxTrainModelSize) )   
	  	{
//			if (!isRandom) Rcout <<cycle<<" : "<< featuresize <<" : "<<pthr2<<" : "<<minpiri<< " : "<< passtestcout << " : "<< padjs <<"\n";	
			if (vnames[jmax] != inname)
			{
				unsigned int vidx = lockup[vnames[jmax]];
				gfrm1 = frm1+" + ";
				gfrm1 = gfrm1+ovnames[jmax];
				mysampleMatbase = join_rows(mysampleMatbase,mysample.col(vidx));
				bestmodel = modelFittingFunc(myoutcome,mysampleMatbase,type);
				bestTrainResiduals = residualForFRESAFunc(bestmodel,mysampleMatbase,"",type,train_outcome);
				if (loops > 1) 
				{
					myTestsampleMatbase = join_rows(myTestsampleMatbase,myTestsample.col(vidx));
					bestResiduals = residualForFRESAFunc(bestmodel,myTestsampleMatbase,"",type,test_outcome);
				}
				error = mean(abs(bestTrainResiduals))*istdoutcome;
				frm1 = gfrm1;	
				mynamesLoc.push_back(jmax);
				inserted = inserted + 1;
				if (error>tol) 
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

extern "C" SEXP ForwardResidualModelCpp(SEXP _size, SEXP _fraction,SEXP _pvalue, SEXP _loops,  SEXP _covariates, SEXP _Outcome, SEXP _variableList, SEXP _maxTrainModelSize, SEXP _type, SEXP _timeOutcome, SEXP _testType, SEXP _dataframe,SEXP _colnames,SEXP _featuresize,SEXP _cores)
{
try {   //R_CStackLimit=(uintptr_t)-1;
		std::string type = as<std::string>(_type);
		std::string testType = as<std::string>(_testType);
#ifdef _OPENMP
//		Rcout<<"::============ Available CPU Cores(Threads) : "<<omp_get_num_procs()<<" ===Requested Threads : " << as<unsigned int>(_cores) <<endl;
		omp_set_num_threads(as<unsigned int>(_cores));
//		omp_set_num_threads(1);
#endif	

		unsigned int size = as<unsigned int>(_size); 
		double fraction = as<double>(_fraction);
		double pvalue = as<double>(_pvalue);
		int loops = as<int>(_loops);
		std::string Outcome = as<std::string>(_Outcome); 
		CharacterVector colnames(_colnames);
	    int maxTrainModelSize = Rcpp::as<int>(_maxTrainModelSize);
		std::string timeOutcome = as<std::string>(_timeOutcome);
		Rcpp::NumericMatrix DF(_dataframe);
    	arma::mat dataframe(DF.begin(), DF.rows(), DF.cols(), false);
		std::vector<std::string>covari=as<std::vector<std::string> >(_covariates);
		std::string baseForm=Outcome;
		std::map<std::string, int> lockup;
		std::vector<std::string> formulas;
		std::vector<int> mynames;
		vec bbasevalues=zeros<vec>(size);
		mat pthrshold=ones<mat>(size,loops);;
		unsigned int featuresize=Rcpp::as<unsigned int>(_featuresize);
		std::vector<std::string>ovnames=as<std::vector<std::string> >(CharacterVector(_variableList)); 
		double pthr = pvalue;

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
		baseForm = baseForm + " ~";
		for(unsigned int i=0;i<covari.size();i++)
		{
			baseForm = baseForm+" + "+covari[i];
		}
		int sizecases   =dataframe.n_rows;
		int totSamples = (int)(fraction*sizecases);
		double pthr2 = 1.0-R::pnorm(std::sqrt(fraction)*std::abs(qnorm(pthr,0.0,1.0)),0.0,1.0,1,0);
	    if (pthr2>0.1) pthr2 = 0.1;
		if (size > ovnames.size()) size = ovnames.size();
//		Rcout<<pvalue<<"pv: "<<pthr2<<" : "<<baseForm <<" : "<< size <<" : " << featuresize <<"\n";


		
		
#pragma omp parallel for schedule(dynamic) ordered shared(mynames,formulas,bbasevalues,pthrshold)
    	for (int doOver=0; doOver<loops;doOver++)
	 	{   
			std::vector <int> mynamesLoc;
			std::string  frm;
			vec basevalues=ones<vec>(size);
			vec pthrs=zeros<vec>(size);
			mynamesLoc.clear();
//			int doBootstrapSample = 1+loops*((doOver % 20)!=0);
			frm = residualFowardSelection(size,totSamples,pthr2,testType,loops,
			Outcome,timeOutcome,type,maxTrainModelSize,dataframe,
			ovnames,lockup,baseForm,mynamesLoc,covari,basevalues,israndom,featuresize,pthrs);
			if ((doOver%100)==99)
			{				
#pragma omp critical
				Rcout << ".";	
//				Rcout<<doOver<<" : "<<frm<<std::endl;	
			}
#pragma omp critical
{
				mynames.insert(mynames.end(), mynamesLoc.begin(), mynamesLoc.end());
				formulas.push_back(frm);
				for (unsigned int i=0;i<size;i++) 
				{
					bbasevalues[i] += basevalues[i];
					pthrshold(i,doOver) = pthrs[i];
				}

}
  		}
		vec med_pthrshold=zeros<vec>(size);
		for (unsigned int i=0;i<size;i++) 
		{
			bbasevalues[i] = bbasevalues[i]/loops;
			rowvec drow = pthrshold.row(i);
			uvec whogz = find(drow);
			if (whogz.n_elem > 0 ) med_pthrshold[i] = median(drow(whogz));
			else
			{
				if ( i >0 ) med_pthrshold[i]=med_pthrshold[i-1];
			}
		}
		
  		if (mynames.size() == 0) mynames.push_back(0);
	    List result = List::create(Named("mynames")=wrap(mynames),Named("formula.list")=wrap(formulas),Named("Base.values")=wrap(bbasevalues),Named("p.thresholds")=wrap(med_pthrshold));
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
