#include "CHydro.h"

CHydro::CHydro(parameterMap* pM) {
	pMap = pM;

	mDataRoot = parameter::getS(*pMap,"HYDRO_OUTPUT_DATAROOT","./");

	int sLength = mDataRoot.size();
	if (mDataRoot[sLength-1] != '/') mDataRoot.append("/");
  
	if (!checkDirectoryExistence(mDataRoot)){
		printf("\nHYDRO_OUTPUT_DATAROOT\n (((%s))) does not exist! Aborting.\n",
			   mDataRoot.c_str());
		exit(1);
	}
	
	mDebug = parameter::getB(*pMap,"HYDRO_DEBUG",false);
	mIoIntegrals = parameter::getB(*pMap,"HYDRO_IO_INTEGRALS", false);
	mIoSlices = parameter::getB(*pMap,"HYDRO_IO_SLICES",false);
	mIoSpots = parameter::getB(*pMap,"HYDRO_IO_SPOTS",false);
	mIoTechqm = parameter::getB(*pMap,"HYDRO_IO_TECHQM",false);
	mIoOscarFull = parameter::getB(*pMap,"HYDRO_IO_OSCARFULL",false);
	mIoOscarHyper = parameter::getB(*pMap,"HYDRO_IO_OSCARHYPER",false);
	mIoFull = parameter::getB(*pMap,"HYDRO_IO_FULL",false);
	mIoSpectra = parameter::getB(*pMap,"HYDRO_IO_SPECTRA",false);

	mSawtooth = parameter::getB(*pMap,"HYDRO_SAWTOOTH",false);
	mRK2  = parameter::getB(*pMap,"HYDRO_RK2",true);
	mRK4  = parameter::getB(*pMap,"HYDRO_RK4",false);
	mLinT = parameter::getB(*pMap,"HYDRO_LINT",false);
	mLogT = parameter::getB(*pMap,"HYDRO_LOGT",false);
	mLogSinhT = parameter::getB(*pMap,"HYDRO_LOGSINHT",true);

	mTStep = parameter::getI(*pMap,"HYDRO_TSTEP",-1);
	mXSize = parameter::getI(*pMap,"HYDRO_XSIZE",60);
	mYSize = parameter::getI(*pMap,"HYDRO_YSIZE",60);
	mNSize = parameter::getI(*pMap,"HYDRO_NSIZE",20);
	mPrintX = parameter::getI(*pMap,"HYDRO_PRINTX",0);
	mPrintY = parameter::getI(*pMap,"HYDRO_PRINTY",0);
	mPrintN = parameter::getI(*pMap,"HYDRO_PRINTN",0);

	mIoSliceTStep  = parameter::getI(*pMap,"HYDRO_IO_SLICES_TSTEP",10);
	mIoTechqmTStep = parameter::getI(*pMap,"HYDRO_IO_TECHQM_TSTEP",10);
	mIoOscarTStep  = parameter::getI(*pMap,"HYDRO_IO_OSCAR_TSTEP",10);
	mIoSpotsTStep  = parameter::getI(*pMap,"HYDRO_IO_SPOTS_TSTEP",1);
	mIoFullTStep   = parameter::getI(*pMap,"HYDRO_IO_FULL_TSTEP",100);
	
	mBjorken = parameter::getB(*pMap,"HYDRO_BJORKEN",true);
	mPureBjorken = parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false);
	mHalving = parameter::getB(*pMap, "HYDRO_HALVING", false);
	mOctant = parameter::getI(*pMap, "HYDRO_OCTANT", 3);

	mInitNS = parameter::getD(*pMap,"HYDRO_INIT_NS",0.);
	mInitNSTxxP = parameter::getD(*pMap,"HYDRO_INIT_NS_TXXP",1.);
	
	mInitFlow = parameter::getD(*pMap,"HYDRO_INIT_FLOW",0.);

	mOldFile = parameter::getB(*pMap,"HYDRO_OLD_FILE",false);
	mOldFileName = parameter::getS(*pMap,"HYDRO_OLD_FILENAME","");

	mFoTemp = parameter::getD(*pMap,"HYDRO_FOTEMP",0.165);
	mT0 = parameter::getD(*pMap,"HYDRO_T0",1.0);
	mDt = parameter::getD(*pMap,"HYDRO_DT",0.03);
	mE0 = parameter::getD(*pMap,"HYDRO_E0",1.0);

	
	nPrintSlice = 0; nPrintFull = 0;
	
		//	zeroPointers();
}

CHydro::~CHydro(){
	delete onMesh;
	delete offMesh;
	delete tempMesh;
	delete mEos;
	
	if (mRK4) {
		delete k1;
		delete k2;
	}
}

void CHydro::zeroPointers(){
	fEX = NULL;
	fEN = NULL;
}

int CHydro::runHydro() {
	time(&start);

		// Open IO Files
	if (mIoTechqm) {
		string temp = mDataRoot;
		temp.append("fTQMX.txt");
		fTQMX = fopen(temp.c_str(),"w");
		
		temp = mDataRoot;
		temp.append("fTQMY.txt");
		fTQMY = fopen(temp.c_str(),"w");
		
		temp = mDataRoot;
		temp.append("fTQMR.txt");
		fTQMR = fopen(temp.c_str(),"w");
	}
  
		// make main mesh
		//	printf("making main mesh....\n");
	if (mOldFile) {
		if (mOctant == 3)
			onMesh = new octMesh(pMap,mOldFileName.c_str());
		else 
			onMesh = new fullMesh(pMap,mOldFileName.c_str());
	}
	else {
		if (mOctant == 3)
			onMesh = new octMesh(pMap); 
		else 
			onMesh = new fullMesh(pMap);
	}
	
	mEos = new CEos(pMap);
    onMesh->setEos(mEos);
	
	if (parameter::getB(*pMap,"HYDRO_KLN_INPUT",false))
		onMesh->fixEntropyMesh();

		//	printf("filling active cells...\n");
	onMesh->fillActiveCells();
	
		// alter initial conditions
	if (mInitFlow != 0.) 
		onMesh->addInitialFlow();
	
	if (mInitNS != 0. || mInitNSTxxP != 0.)   
		onMesh->initNS();
	
		// make additional meshes
		//	printf("making additional meshes...\n");
	if (mOctant == 3){
		offMesh  = new octMesh(onMesh);
		tempMesh = new octMesh(onMesh);
		if (mRK4) {
			k1 = new octMesh(onMesh); 
			k2 = new octMesh(onMesh);
		}
	}
	else {
		offMesh  = new fullMesh(onMesh);
		tempMesh = new fullMesh(onMesh);
		if (mRK4) {
			k1 = new fullMesh(onMesh); 
			k2 = new fullMesh(onMesh);
		}
	}
	
		//	printf("active matching between meshes...\n");
	offMesh->fillActiveCells();
	tempMesh->fillActiveCells();
	onMesh->copyActive(offMesh,tempMesh);
	
	if (mRK4) {
		k1->fillActiveCells();
		k2->fillActiveCells();
		onMesh->copyActive(k1,k2);
	}
	
#ifdef HYDRO_BOOST_THREADS
	offMesh->setEos(mEos);
	tempMesh->setEos(mEos);
	CCellHelper* mHelper = new CCellHelper();
	tempMesh->setHelper(mHelper);
	onMesh->setHelper(mHelper);
	offMesh->setHelper(mHelper);
#endif
	
	printActive();
	
#ifdef HYDRO_XDMF_OUTPUT
	StartXdmf();
	WriteXdmf();
#endif
	
	sLoss=0.; 
	eLoss=0.;
	
	CMesh *oldMesh;
	FILE *fOSU;
	
		// first IO dump
		//printf("first i/o dump...\n");
	if (mIoIntegrals) {
		openFileIntegrals();
		onMesh->calcIntegrals();
		double intS = onMesh->integralS();
		intS0 = intS;
		double intE = onMesh->integralE();
		intE0 = intE;

		fprintf(fIS,"%0.9g 0.\n", mT0);
		fprintf(fIE,"%0.9g 0.\n", mT0);
		
		fprintf(fIEP,"0. %0.6g\n", onMesh->integralEp());
		fprintf(fIEX,"0. %0.6g\n", onMesh->integralEx());
		fprintf(fIVR,"0. %0.6g\n", onMesh->averageVr());
	}
	if (mIoSpots) {
		openFileSpots();
		printSpots();
	}
	if (mIoOscarFull) 
		openOscarFull();
	if (mIoOscarHyper) {
		openOscarHyper();
		if (mOctant == 3)
			oldMesh = new octMesh(tempMesh);
		else 
			oldMesh = new fullMesh(tempMesh);

		string osuFN("/osuFOS.dat");
		osuFN = mDataRoot + osuFN;
		fOSU = fopen(osuFN.c_str(),"w");
		oldMesh->setOsuFos(fOSU);
	}
	if (mIoSlices) printSlices();
	if (mIoFull && false)
		printMesh();
	
    time(&now);
	
		// time to complete estimate
	if (!mPureBjorken) 
		printf("Expected Run Time < %0.3g sec.  \n", 
			   (mXSize*mYSize*mNSize*150.)/2.7E5);
	else 
		printf("Expected Run Time < %0.3g sec.  \n",
			   (mXSize*mYSize*150)/5.E5);
	fflush(stdout);
    printf("Time Step 0 [%0.6g](elapsed time %0.6g sec) .",
		   onMesh->getTau(),difftime(now,start)); fflush(stdout);
	
	nPrintUpdate=1;
	mTemp0 = tempMesh->getT(0,0,0);
	int t=1;
	
	if (mTemp0 < mFoTemp) 
		t = mTStep+1;
		
    for (;(t<=mTStep || mTStep==-1);t++) {
		onMesh->deaden();

		if (mSawtooth) 
			onMesh->forward(tempMesh,offMesh,mDt);
		else if (mRK2) {
			onMesh->copyActive(offMesh,tempMesh);

			onMesh->forward(offMesh,mDt/2.);

/*
			CMesh* m = tempMesh;
			tempMesh = offMesh;
			printSlices();
			tempMesh = m;
*/
			
			onMesh->forward(tempMesh,offMesh,mDt);
		}
		else if (mRK4) {
			onMesh->forward(offMesh,mDt/2.); 
			onMesh->forward(k1,offMesh,mDt/2.);
			onMesh->forward(k2,k1,mDt/2.);
			k2->setTau(k1->getTau());     //k[2] mesh is actually at half step
			onMesh->forward(tempMesh,k2,mDt);
			onMesh->forward(offMesh,k1,k2,tempMesh,mDt);
		}

		if (!tempMesh->detectCrash()) {
			printf("\n\n********** WARNING -- isHydro3 resulted in crash at t=%d.************\n Dumping pre-crash Mesh...... \n",t);
			if (mIoSlices)
				printSlices();
			
			if (mIoSpots)
				closeFileSpots();
			
			if (mIoIntegrals)
				closeFileIntegrals();
			
			if (mIoOscarHyper)
				fclose(fOscarHyper);
			
			if (mIoOscarFull)
				fclose(fOscarFull);
			
			return 1;
		}

		if (tempMesh->getT(0,0,0) < mFoTemp) 
			if (!tempMesh->anyActiveCells())	
				break;

		if (mIoIntegrals) 
			printIntegrals();

		if (mIoSlices && t%mIoSliceTStep==0) 
			printSlices();

		/*
		if (tempMesh->getTau() > 1.6 && tempMesh->getTau() < 1.6 + mDt) {
			printf("\n1.6 print is number %d (%g)\n",nPrintSlice,tempMesh->getTau());
			printSlices();
		}
		*/

		if (mIoFull && t%mIoFullTStep==0) 
			printMesh();
		
		if (mIoSpots && t%mIoSpotsTStep == 0) 
			printSpots();

		if (mIoTechqm && t%mIoTechqmTStep == 0)
			printTQM();

		if (mIoOscarFull && t%mIoOscarTStep == 0)
			printOscarFull(t);

		if (mIoOscarHyper && t%mIoOscarTStep == 0){
			tempMesh->genFOS(oldMesh);
			oldMesh->copyMesh(tempMesh);
			
			if (!printOscarHyper(t)){
				printf("\nILLICIT FREEZEOUT SURFACE at t=%g\n",tempMesh->getTau());
				mIoOscarHyper=false;  // quit trying to print FOS
			}
		}
		
		if (mHalving)
			if (tempMesh->getTau() - onMesh->getTau() > 2.*mDt*mT0) {
				printf("!");
				mDt /= 2.;
			}

		if ( nPrintUpdate < 50.*(pow(tempMesh->getT(0,0,0),-2)-pow(mTemp0,-2))/(pow(mFoTemp,-2)-pow(mTemp0,-2))) {
			if (nPrintUpdate%10==0){
				time(&now);
				printf("\nAt Time %d: [%g] -- elapsed time %g sec .",t,tempMesh->getTau(),difftime(now,start)); 
				fflush(stdout);
			}
			else {
				printf(".");
				fflush(stdout);
			}
			nPrintUpdate++;
		}

		if (mSawtooth) {
			deadMesh = onMesh;
			onMesh = offMesh;
			offMesh = tempMesh;
			tempMesh = deadMesh;
		}
		else if (mRK2) {
			deadMesh = onMesh;
			onMesh = tempMesh;
			tempMesh = deadMesh;
		}
	  // no reassignments for RK4
		
	}
		// if we quit early, undo mesh reassignments
	if (mRK2 && tempMesh->getT(0,0,0) > mFoTemp) {
		deadMesh = tempMesh;
		tempMesh = onMesh;
		onMesh = deadMesh;
	}

	if (mIoSlices)
		printSlices();
	
	if (mIoFull){
		printMesh();
#ifdef HYDRO_XDMF_OUTPUT
		StopXdmf();
#endif
	}
	
	if (mIoTechqm)
		printTQM();
	
	if (mIoSpots) {
		printSpots();
		closeFileSpots();
	}

	if (mIoIntegrals)
		closeFileIntegrals();

	if (mIoOscarHyper)	{
		fclose(fOscarHyper);
		if (tempMesh->getT(0,0,0) < mFoTemp)
			tempMesh->genFOS(oldMesh);
		fclose(fOSU);
		delete oldMesh;
	}

	if (mIoOscarFull) {
		printOscarFull(t);
		fclose(fOscarFull);
	}
	
	time(&now);
	printf("\nHydro Finished %d [%g]...... (Total Time: %0.3g sec)\n",
		   t,tempMesh->getTau(),difftime(now,start)); 
	fflush(stdout);
	
	if (mIoSpectra) {
		time(&start);
		printf("\nCalculating Spectra (Exp Time %0.3g s)....",mXSize*mYSize*mNSize*1500./7.E6); fflush(stdout);
		
		string fn = "fdNdY.txt";
		fn = mDataRoot + fn;
		FILE *fdNdY = fopen( fn.c_str(),"w");
		
		const int nRap = 10;
		double dRap = 0.4;   
		double dNdY[nRap];
		
		tempMesh->calcdNdY(0.1396, 1.0, dNdY, dRap, nRap);
		
		for (int i=0;i<nRap;i++)
			fprintf(fdNdY,"%g %g\n",dRap*i,dNdY[i]);
		
		time(&now);
		printf("Finished! [%g s]\n",difftime(now,start));
	}
	
	tempMesh->writeIEPvEta();
	
#ifdef HYDRO_BOOST_THREADS
	delete mHelper;
#endif
	return 0;
}

void CHydro::spectraFromFile() {
	
	if (mOctant == 3)
		tempMesh = new octMesh(pMap,mOldFileName.c_str());
	else 
		tempMesh = new fullMesh(pMap,mOldFileName.c_str());
	
	mEos = new CEos(pMap);
    onMesh->setEos(mEos);
	
	time(&start);
	printf("\nCalculating Spectra (Exp Time %0.3g s)....",mXSize*mYSize*mNSize*1500./7.E6); fflush(stdout);
	
	string fn = "fdNdY.txt";
	fn = mDataRoot + fn;
	FILE *fdNdY = fopen( fn.c_str(),"w");
	
	const int nRap = 10;
	double dRap = 0.4;   
	double dNdY[nRap];
	
	tempMesh->calcdNdY(0.1396, 1.0, dNdY, dRap, nRap);
	
	for (int i=0;i<nRap;i++)
		fprintf(fdNdY,"%g %g\n",dRap*i,dNdY[i]);
	
	time(&now);
	printf("Finished! [%g s]\n",difftime(now,start));
}

void CHydro::openOscarHyper() {
		// 1
		//  fOscarHyper = fopen("fOscarHyper.OSCAR2008H","w");
	string fn = parameter::getS(*pMap,"HYDRO_IO_OSCARHYPER_FN","freezeout.dat");
	fn = mDataRoot + fn;
	fOscarHyper = fopen( fn.c_str(),"w");
	fprintf(fOscarHyper,"OSCAR2008H  ");
	
	if ( parameter::getD(*pMap,"HYDRO_SVRATIO",0.0) > 0. || parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
		fprintf(fOscarHyper,"viscous     ");
	else 
		fprintf(fOscarHyper,"ideal       ");
	fprintf(fOscarHyper,"final_hs     \n");
	
		// 2
	fprintf(fOscarHyper,"INIT: ");
	int mA = parameter::getD(*pMap,"GLAUBER_A",197.); int mB = parameter::getD(*pMap,"GLAUBER_B",0.);
	fprintf(fOscarHyper," Glauber - rho = %g Gev/fm^3, xi = %g, sigma = %fb, A = %f, b = %f fm\n",
			parameter::getD(*pMap,"GLAUBER_Rho0",0.16)*1.,parameter::getD(*pMap,"GLAUBER_Xi",0.3)*1.,
			parameter::getD(*pMap,"GLAUBER_Sigma",0.4)*10.,parameter::getD(*pMap,"GLAUBER_A",197.),parameter::getD(*pMap,"GLAUBER_B",0.0));
	
		// 3
	fprintf(fOscarHyper,"EOS: ");
	if ( parameter::getB(*pMap,"EQOFST_LATEOS",true))
		fprintf(fOscarHyper," Lat-HRG interpolated, 0 Mev Shift, p4 action\n");
	else 
		fprintf(fOscarHyper," ideal gas of massless pions\n");
	
		// 4
	fprintf(fOscarHyper,"CHARGES: ");
	fprintf(fOscarHyper,"none\n");
	
		// 5
	double temp = parameter::getD(*pMap,"HYDRO_FOTEMP",0.13);
	fprintf(fOscarHyper,"HYPER:  hyper T=%g MeV isotherm \n",temp);
	
		// 6
	fprintf(fOscarHyper,"GEOM: ");
	if ( parameter::getB(*pMap,"HYDRO_BJORKEN",true)) 
		if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false)) 
			fprintf(fOscarHyper,"scaling2d\n");
		else 
			fprintf(fOscarHyper,"3d\n");
		else 
			if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false)) 
				fprintf(fOscarHyper,"2d\n");
			else 
				fprintf(fOscarHyper,"3d-cart");
	
		// 7
	fprintf(fOscarHyper,"GRID:  Euler\n");
	if (mOctant == 0) {
		fprintf(fOscarHyper,"%d %d %d %d 0 5 3\n",
				parameter::getI(*pMap,"HYDRO_TSTEP",100),2*parameter::getI(*pMap,"HYDRO_XSIZE",60)+1,
				2*parameter::getI(*pMap,"HYDRO_YSIZE",60)+1,2*parameter::getI(*pMap,"HYDRO_NSIZE",20)+1);
		fprintf(fOscarHyper,"%f ",parameter::getD(*pMap,"HYDRO_T0",1.0));
		if (parameter::getB(*pMap,"HYDRO_LINT",false)) 
			fprintf(fOscarHyper,"%f ",parameter::getI(*pMap,"HYDRO_TSTEP",100)
					*parameter::getD(*pMap,"HYDRO_DT",0.1)+parameter::getD(*pMap,"HYDRO_T0",1.0));
			// else FIXME
		
		fprintf(fOscarHyper,"%f %f %f %f %f %f\n",
				-parameter::getD(*pMap,"HYDRO_DX",0.1)*parameter::getI(*pMap,"HYDRO_XSIZE",60),
				parameter::getD(*pMap,"HYDRO_DX",0.1)*parameter::getI(*pMap,"HYDRO_XSIZE",60),
				-parameter::getD(*pMap,"HYDRO_DY",0.1)*parameter::getI(*pMap,"HYDRO_YSIZE",60),
				parameter::getD(*pMap,"HYDRO_DY",0.1)*parameter::getI(*pMap,"HYDRO_YSIZE",60),
				-parameter::getD(*pMap,"HYDRO_DN",0.1)*parameter::getI(*pMap,"HYDRO_NSIZE",20),
				parameter::getD(*pMap,"HYDRO_DN",0.1)*parameter::getI(*pMap,"HYDRO_NSIZE",20));
	}
	else {
		fprintf(fOscarHyper,"%d %d %d %d 0 5 3\n",
				parameter::getI(*pMap,"HYDRO_TSTEP",100),parameter::getI(*pMap,"HYDRO_XSIZE",60)+1,
				parameter::getI(*pMap,"HYDRO_YSIZE",60)+1,parameter::getI(*pMap,"HYDRO_NSIZE",20)+1);
		fprintf(fOscarHyper,"%f ",parameter::getD(*pMap,"HYDRO_T0",1.0));
		if (parameter::getB(*pMap,"HYDRO_LINT",false)) 
			fprintf(fOscarHyper,"%f ",parameter::getI(*pMap,"HYDRO_TSTEP",100)
					*parameter::getD(*pMap,"HYDRO_DT",0.1)+parameter::getD(*pMap,"HYDRO_T0",1.0));
			// else FIXME
		
		fprintf(fOscarHyper,"0. %f 0. %f 0. %f\n",parameter::getD(*pMap,"HYDRO_DX",0.1)*parameter::getI(*pMap,"HYDRO_XSIZE",60),
				parameter::getD(*pMap,"HYDRO_DY",0.1)*parameter::getI(*pMap,"HYDRO_YSIZE",60),
				parameter::getD(*pMap,"HYDRO_DN",0.1)*parameter::getI(*pMap,"HYDRO_NSIZE",20));
	}
	
		// 8
	fprintf(fOscarHyper,"VISCOSITY:  ");
	if ( parameter::getD(*pMap,"HYDRO_SVRATIO",0.0) > 0.)
		if (parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
			fprintf(fOscarHyper,"shear and bulk viscosity\n");
		else 
			fprintf(fOscarHyper,"shear viscosity only\n");
		else 
			if (parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
				fprintf(fOscarHyper,"bulk viscosity only\n");
			else 
				fprintf(fOscarHyper,"none");
	
		//  fprintf(fOscarHyper,"COMM:  Dissipative Quantities: eta, tau_pi, s\n");
	fprintf(fOscarHyper,"COMM:  Reflective symmetry in x,y, and eta assumed\n");
	
	fprintf(fOscarHyper,"END_OF_HEADER\n");
	
	printOscarHyper(0);
}

void CHydro::openOscarFull() {
	
		// 1
	fOscarFull = fopen("fOscarFull.OSCAR2008H","w");
	fprintf(fOscarFull,"OSCAR2008H  ");
	if ( parameter::getD(*pMap,"HYDRO_SVRATIO",0.0) > 0. || parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
		fprintf(fOscarFull,"viscous     ");
	else 
		fprintf(fOscarFull,"ideal       ");
	fprintf(fOscarFull,"history     \n");
	
		// 2
	fprintf(fOscarFull,"INIT: ");
	int mA = parameter::getD(*pMap,"GLAUBER_A",197.); int mB = parameter::getD(*pMap,"GLAUBER_B",197.);
	fprintf(fOscarFull," Glauber - rho = %g Gev/fm^3, xi = %g, sigma = %fb, A = %f, b = %f fm\n",
			parameter::getD(*pMap,"GLAUBER_Rho0",0.16)*1.,parameter::getD(*pMap,"GLAUBER_Xi",0.3)*1.,
			parameter::getD(*pMap,"GLAUBER_Sigma",0.4)*10.,parameter::getD(*pMap,"GLAUBER_A",197.),parameter::getD(*pMap,"GLAUBER_B",0.0));
	
		// 3
	fprintf(fOscarFull,"EOS: ");
	if (parameter::getB(*pMap,"EQOFST_LATEOS",true))
		fprintf(fOscarFull," Lat-HRG interpolated, 0 Mev Shift, p4 action\n");
	else 
		fprintf(fOscarFull," ideal gas of quarks and gluons\n");
	
		// 4
	fprintf(fOscarFull,"CHARGES: ");
	fprintf(fOscarFull,"none\n");
  	
		// 5
	fprintf(fOscarFull,"HYPER:  full evolution\n");
	
		// 6
	fprintf(fOscarFull,"GEOM: ");
	if (parameter::getB(*pMap,"HYDRO_BJORKEN",true)) 
		if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false)) 
			fprintf(fOscarFull,"scaling2d\n");
		else 
			fprintf(fOscarFull,"3d\n");
		else 
			if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false)) 
				fprintf(fOscarFull,"2d\n");
			else fprintf(fOscarFull,"3d-cart");
	
		// 7
	fprintf(fOscarFull,"GRID:  Euler\n");
	if (mOctant == 3) {
		fprintf(fOscarFull,"%d %d %d %d 0 5 3\n",
				parameter::getI(*pMap,"HYDRO_TSTEP",100),parameter::getI(*pMap,"HYDRO_XSIZE",60)+1,
				parameter::getI(*pMap,"HYDRO_YSIZE",60)+1,parameter::getI(*pMap,"HYDRO_NSIZE",20)+1);
		fprintf(fOscarFull,"%f ",parameter::getD(*pMap,"HYDRO_T0",1.0));
		if (parameter::getB(*pMap,"HYDRO_LINT",false)) 
			fprintf(fOscarFull,"%f ",parameter::getI(*pMap,"HYDRO_TSTEP",100)*parameter::getD(*pMap,"HYDRO_DT",0.1)+parameter::getD(*pMap,"HYDRO_T0",1.0));
			// else FIXME
		
		fprintf(fOscarFull,"0. %f 0. %f 0. %f\n",parameter::getD(*pMap,"HYDRO_DX",0.1)*parameter::getI(*pMap,"HYDRO_XSIZE",60),
				parameter::getD(*pMap,"HYDRO_DY",0.1)*parameter::getI(*pMap,"HYDRO_YSIZE",60),
				parameter::getD(*pMap,"HYDRO_DN",0.1)*parameter::getI(*pMap,"HYDRO_NSIZE",20));
	}
	else {printf("!!!!!!!!!! OSCAR OUTPUT FAILING !!!!!!!!!!"); return;}
	
		// 8
	fprintf(fOscarFull,"VISCOSITY:  ");
	if (parameter::getD(*pMap,"HYDRO_SVRATIO",0.0) > 0.)
		if (parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
			fprintf(fOscarFull,"shear and bulk viscosity\n");
		else fprintf(fOscarFull,"shear viscosity only\n");
		else 
			if (parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
				fprintf(fOscarFull,"bulk viscosity only\n");
			else fprintf(fOscarFull,"none");
	
	fprintf(fOscarFull,"COMM:  Dissipative Quantities: a_i (1-5, if SV>0), b (if BV>0)\n");
	fprintf(fOscarFull,"COMM:  Dissipative Qunatities: above will be rescaled to have units GeV/fm^3 (if necessary)\n");
	fprintf(fOscarFull,"COMM:  Transport Coefficients: eta, tau_pi, s\n");
	fprintf(fOscarFull,"COMM:  Reflective symmetry in x,y, and eta assumed\n");
	
	fprintf(fOscarFull,"END_OF_HEADER\n");
	
	printOscarFull(0);
}

void CHydro::openSliceFiles() {
	string temp;
	char numC[2];
	sprintf(numC,"%d",nPrintSlice);
	
	string num;
	if (nPrintSlice > 9){
		num.append(1,numC[0]);
		num.append(1,numC[1]);
	}
	else {
		num.append(1,numC[0]);
	}

	if (nPrintSlice == 99) return;
	
	nPrintSlice++;
	
	temp = mDataRoot;
	temp.append("fEx");  
	temp.append(num);
	temp.append(".txt");
	fEX = fopen(temp.c_str(),"w");
	  
	
	temp = mDataRoot;
	temp.append("fEy");
	temp.append(num);
	temp.append(".txt");
	fEY = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fEz");
	temp.append(num);
	temp.append(".txt");
	fEZ = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fS");
	temp.append(num);
	temp.append(".txt");
	fS = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fT");
	temp.append(num);
	temp.append(".txt");
	fT = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fP");
	temp.append(num);
	temp.append(".txt");
	fP = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fEN");
	temp.append(num);
	temp.append(".txt");
	fEN = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fU");
	temp.append(num);
	temp.append(".txt");
	fU0 = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fUx");
	temp.append(num);
	temp.append(".txt");
	fUx = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fUy");
	temp.append(num);
	temp.append(".txt");
	fUy = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fUz");
	temp.append(num);
	temp.append(".txt");
	fUz = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fUzz");
	temp.append(num);
	temp.append(".txt");
	fUzz = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fULz");
	temp.append(num);
	temp.append(".txt");
	fULz = fopen(temp.c_str(),"w");
			  
	temp = mDataRoot;
	temp.append("fA1x-");
	temp.append(num);
	temp.append(".txt");
	fA1x = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fA2x-");
	temp.append(num);
	temp.append(".txt");
	fA2x = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fA3x-");
	temp.append(num);
	temp.append(".txt");
	fA3x = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fA4x-");
	temp.append(num);
	temp.append(".txt");
	fA4x = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fA5x-");
	temp.append(num);
	temp.append(".txt");
	fA5x = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fB");
	temp.append(num);
	temp.append("x.txt");
	fBx = fopen(temp.c_str(),"w");

	temp = mDataRoot;
	temp.append("fA1y-");
	temp.append(num);
	temp.append(".txt");
	fA1y = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA2y-");
	temp.append(num);
	temp.append(".txt");
	fA2y = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA3y-");
	temp.append(num);
	temp.append(".txt");
	fA3y = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA4y-");
	temp.append(num);
	temp.append(".txt");
	fA4y = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA5y-");
	temp.append(num);
	temp.append(".txt");
	fA5y = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fB");
	temp.append(num);
	temp.append("y.txt");
	fBy = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA1z-");
	temp.append(num);
	temp.append(".txt");
	fA1z = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA2z-");
	temp.append(num);
	temp.append(".txt");
	fA2z = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA3z-");
	temp.append(num);
	temp.append(".txt");
	fA3z = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA4z-");
	temp.append(num);
	temp.append(".txt");
	fA4z = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA5z-");
	temp.append(num);
	temp.append(".txt");
	fA5z = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fB");
	temp.append(num);
	temp.append("z.txt");
	fBz = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fTxx");
	temp.append(num);
	temp.append(".txt");
	fTxx = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fTxxNS");
	temp.append(num);
	temp.append(".txt");
	fTxxNS = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fTyy");
	temp.append(num);
	temp.append(".txt");
	fTyy = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fTzz");
	temp.append(num);
	temp.append(".txt");
	fTzz = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fTxy");
	temp.append(num);
	temp.append(".txt");
	fTxy = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxz");
	temp.append(num);
	temp.append(".txt");
	fTxz = fopen(temp.c_str(),"w");
		
	temp = mDataRoot;
	temp.append("fTyz");
	temp.append(num);
	temp.append(".txt");
	fTyz = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fTrr");
	temp.append(num);
	temp.append(".txt");
	fTrr = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fTpp");
	temp.append(num);
	temp.append(".txt");
	fTpp = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fTpr");
	temp.append(num);
	temp.append(".txt");
	fTpr = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTIS");
	temp.append(num);
	temp.append(".txt");
	fTIS = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fTISB");
	temp.append(num);
	temp.append(".txt");
	fTISB = fopen(temp.c_str(),"w");
	  
	temp = mDataRoot;
	temp.append("fFOS");
	temp.append(num);
	temp.append(".txt");
	fFOS = fopen(temp.c_str(),"w");
}

void CHydro::openFileIntegrals() {
    string temp = mDataRoot;
	temp.append("fIE.txt");
    fIE = fopen( temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fIS.txt");
	fIS = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fSLoss.txt");
    fSLoss = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fSLoss2.txt");
    fSLoss2 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fELoss.txt");
    fELoss = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fEp.txt");
	fIEP = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fEx.txt");
	fIEX = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fVr.txt");
	fIVR = fopen(temp.c_str(),"w");
}

void CHydro::openFileSpots() {
	string temp = mDataRoot;
	temp.append("fE.txt");
    fE0 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fT.txt");
    fT0  = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxx.txt");
    fTxx0 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTzz.txt");
    fTzz0 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fUL.txt");
	fUL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fEL.txt");
	fEL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA1L.txt");
	fA1L = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA2L.txt");
	fA2L = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fBL.txt");
	fBL  = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fDzz.txt");
	fDzz = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fFOSL.txt");
	fFOSL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fFOSSST.txt");
	fFOSSST = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxxL.txt");
	fTxxL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxxNSL.txt");
	fTxxNSL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxxML.txt");
	fTxxML = fopen(temp.c_str(),"w");
	
}

void CHydro::closeSliceFiles() {
	fclose(fEX); fclose(fEY); fclose(fEN); fclose(fEZ); fclose(fS); fclose(fT); fclose(fP); 
	fclose(fU0); fclose(fUx); fclose(fUy); fclose(fUz); fclose(fULz); fclose(fUzz);
	fclose(fA1x); fclose(fA2x); fclose(fA3x); fclose(fA4x); fclose(fA5x); fclose(fBx);
	fclose(fA1y); fclose(fA2y); fclose(fA3y); fclose(fA4y); fclose(fA5y); fclose(fBy);
	fclose(fA1z); fclose(fA2z); fclose(fA3z); fclose(fA4z); fclose(fA5z); fclose(fBz);
	fclose(fTxx); fclose(fTxxNS); fclose(fTyy); fclose(fTzz); fclose(fTxy); fclose(fTxz); fclose(fTyz); 
	fclose(fTrr); fclose(fTpp); fclose(fTpr);
	fclose(fTIS); fclose(fTISB); fclose(fFOS);
}

void CHydro::closeFileIntegrals() {
    fclose(fIE);
	fclose(fIS);
	fclose(fSLoss);
    fclose(fELoss);
	fclose(fIEP);
	fclose(fIEX);
	fclose(fIVR);
}

void CHydro::closeFileSpots() {
    fclose(fE0);
    fclose(fT0);
    fclose(fTxx0);
    fclose(fTzz0);
	fclose(fUL);
	fclose(fEL);
	fclose(fA1L);
	fclose(fA2L);
	fclose(fBL);
	fclose(fDzz);
	fclose(fFOSL);
	fclose(fFOSSST);
	fclose(fTxxL);
	fclose(fTxxNSL);
	fclose(fTxxML);
}

void CHydro::printSlices() { 
	openSliceFiles();
	
	int nPrintX = min(mPrintX,tempMesh->getXSize()-1);
	int nPrintY = min(mPrintY,tempMesh->getYSize()-1);
	int lPrintN = min(mPrintN,tempMesh->getNSize()-1);
	
	if (mPureBjorken) lPrintN=0;

    if (!mPureBjorken) {
		int i;
		if (mOctant==3) i=0; else i=-tempMesh->getNSize();
		for (; i<tempMesh->getNSize(); i++) {
			if (!tempMesh->getActive(i,nPrintX,nPrintY)) continue;
			double localN = tempMesh->getX(i,nPrintX,nPrintY,3);			
			fprintf(fEN,"%0.9g %0.9g\n", localN, tempMesh->getE(i,nPrintX,nPrintY));
			
			fprintf(fEZ, "%0.9g %0.9g\n", localN, tempMesh->getE(i,nPrintX,nPrintY));
			fprintf(fUz, "%0.9g %0.9g\n", localN, tempMesh->getTau()*tempMesh->getS(i,nPrintX,nPrintY,3));
			fprintf(fUzz,"%0.9g %0.9g\n", localN, tempMesh->getS(i,nPrintX,nPrintY,3));
			
			fprintf(fA1z,"%0.9g %0.9g\n", localN, tempMesh->getS(i,nPrintX,nPrintY,4));
			fprintf(fA2z,"%0.9g %0.9g\n", localN, tempMesh->getS(i,nPrintX,nPrintY,5));
			fprintf(fA3z,"%0.9g %0.9g\n", localN, tempMesh->getS(i,nPrintX,nPrintY,6));
			fprintf(fA4z,"%0.9g %0.9g\n", localN, tempMesh->getS(i,nPrintX,nPrintY,7));
			fprintf(fA5z,"%0.9g %0.9g\n", localN, tempMesh->getS(i,nPrintX,nPrintY,8));
			fprintf(fBz, "%0.9g %0.9g\n", localN, tempMesh->getS(i,nPrintX,nPrintY,10));
			
			fprintf(fT,"%0.9g %0.9g\n", localN, tempMesh->getT(i,nPrintX,nPrintY));
			
			double P = tempMesh->getP(i,nPrintX,nPrintY);
				//fprintf(fP,"%0.9g %0.9g\n", localN, P);
			
			fprintf(fTxx,"%0.9g %0.9g\n", localN, tempMesh->getPixy(i,nPrintX,nPrintY,1,1)/P);
			fprintf(fTyy,"%0.9g %0.9g\n", localN, tempMesh->getPixy(i,nPrintX,nPrintY,2,2)/P);
			fprintf(fTzz,"%0.9g %0.9g\n", localN, pow(tempMesh->getTau(),2)*tempMesh->getPixy(i,nPrintX,nPrintY,3,3)/P);
			fprintf(fTxy,"%0.9g %0.9g\n", localN, tempMesh->getPixy(i,nPrintX,nPrintY,1,2)/P);
			fprintf(fTyz,"%0.9g %0.9g\n", localN, tempMesh->getTau()*tempMesh->getPixy(i,nPrintX,nPrintY,2,3)/P);
		}
	}
	
	double mFOS[XSIZE+YSIZE][2];
	double mFOSigma[XSIZE+YSIZE][2];
	
	/*
	int fosSize;
	tempMesh->getFOS(mFOS,mFOSigma,fosSize);
	for (int i=0;i<fosSize;i++){
		//printf("%g %g\n",mFOS[i][0], mFOS[i][1]);
		fprintf(fFOS,"%0.6g %0.6g\n", mFOS[i][0], mFOS[i][1]);
	}
	 */
	
	int lPrintY = nPrintY;
	for (int i=0; i<tempMesh->getXSize(); i++) {
		if (!tempMesh->getActive(lPrintN,i,lPrintY)) break;
		tempMesh->selfUpdate(lPrintN,i,lPrintY);
		
		double localX = tempMesh->getX(1);
		fprintf(fEX,"%0.9g %0.9g\n", localX, tempMesh->getE());

		fprintf(fUx,"%0.9g %0.9g\n", localX, tempMesh->getS(1));
		fprintf(fUy,"%0.9g %0.9g\n", localX, tempMesh->getS(2));
		
		fprintf(fA1x,"%0.9g %0.9g\n", localX, tempMesh->getS(4));
		fprintf(fA2x,"%0.9g %0.9g\n", localX, tempMesh->getS(5));
		fprintf(fA3x,"%0.9g %0.9g\n", localX, tempMesh->getS(6));
		fprintf(fA4x,"%0.9g %0.9g\n", localX, tempMesh->getS(7));
		fprintf(fA5x,"%0.9g %0.9g\n", localX, tempMesh->getS(8));
		
		fprintf(fBx, "%0.9g %0.9g\n", localX, tempMesh->getS(10));
		
		double P = tempMesh->getP();
		fprintf(fT,"%0.9g %0.9g\n", localX, tempMesh->getT());
			//		fprintf(fP,"%0.9g %0.9g\n", localX, P);
		
		fprintf(fTxx,"%0.9g %0.9g\n", localX, tempMesh->getPixy(1,1)/P);
		fprintf(fTrr,"%0.9g %0.9g\n", localX, tempMesh->getTRR()/P);
		fprintf(fTpp,"%0.9g %0.9g\n", localX, tempMesh->getTPhiPhi()/P);
		fprintf(fTpr,"%0.9g %0.9g\n", localX, tempMesh->getTRPhi()/P);
		
		fprintf(fTxxNS,"%0.9g %0.9g\n",localX, tempMesh->getPixyNS(1,1)/P);
		fprintf(fTyy,"%0.9g %0.9g\n",  localX, tempMesh->getPixy(2,2)/P);
		fprintf(fTzz,"%0.9g %0.9g\n",  localX, pow(tempMesh->getTau(),2)*tempMesh->getPixy(3,3)/P);
		
		fprintf(fTxy,"%0.9g %0.9g\n", localX, tempMesh->getPixy(1,2)/P);
		fprintf(fTxz,"%0.9g %0.9g\n", localX, tempMesh->getTau()*tempMesh->getPixy(1,3)/P);
		fprintf(fTyz,"%0.9g %0.9g\n", localX, tempMesh->getTau()*tempMesh->getPixy(2,3)/P);
		
		fprintf(fTIS, "%0.9g %0.9g\n", localX, tempMesh->getTIS());
		fprintf(fTISB,"%0.9g %0.9g\n", localX, tempMesh->getTISB());
		
		if (mPureBjorken) 
			fprintf(fUz,"%0.9g %0.9g\n", tempMesh->getX(2), tempMesh->getS(3));
	}
	
	if (mOctant==0)
		for (int i=0; i<tempMesh->getXSize(); i++) {
			if (!tempMesh->getActive(lPrintN,-i,lPrintY)) break;
			tempMesh->selfUpdate(lPrintN,-i,lPrintY);
			
			double localX = tempMesh->getX(1);
			fprintf(fEX,"%0.9g %0.9g\n", localX, tempMesh->getE());
			
			fprintf(fUx,"%0.9g %0.9g\n", localX, tempMesh->getS(1));
			fprintf(fUy,"%0.9g %0.9g\n", localX, tempMesh->getS(2));
			
			fprintf(fA1x,"%0.9g %0.9g\n", localX, tempMesh->getS(5));
			fprintf(fA2x,"%0.9g %0.9g\n", localX, tempMesh->getS(6));
			fprintf(fA3x,"%0.9g %0.9g\n", localX, tempMesh->getS(7));
			fprintf(fA4x,"%0.9g %0.9g\n", localX, tempMesh->getS(8));
			fprintf(fA5x,"%0.9g %0.9g\n", localX, tempMesh->getS(9));
			
			fprintf(fBx, "%0.9g %0.9g\n", localX, tempMesh->getS(10));
			
			double P = tempMesh->getP();
			fprintf(fT,"%0.9g %0.9g\n", localX, tempMesh->getT());
				//		fprintf(fP,"%0.9g %0.9g\n", localX, P);
			
			fprintf(fTxx,"%0.9g %0.9g\n", localX, tempMesh->getPixy(1,1)/P);
			fprintf(fTrr,"%0.9g %0.9g\n", localX, tempMesh->getTRR()/P);
			fprintf(fTpp,"%0.9g %0.9g\n", localX, tempMesh->getTPhiPhi()/P);
			fprintf(fTpr,"%0.9g %0.9g\n", localX, tempMesh->getTRPhi()/P);
			
			fprintf(fTxxNS,"%0.9g %0.9g\n",localX, tempMesh->getPixyNS(1,1)/P);
			fprintf(fTyy,"%0.9g %0.9g\n",  localX, tempMesh->getPixy(2,2)/P);
			fprintf(fTzz,"%0.9g %0.9g\n",  localX, pow(tempMesh->getTau(),2)*tempMesh->getPixy(3,3)/P);
			
			fprintf(fTxy,"%0.9g %0.9g\n", localX, tempMesh->getPixy(1,2)/P);
			fprintf(fTxz,"%0.9g %0.9g\n", localX, tempMesh->getTau()*tempMesh->getPixy(1,3)/P);
			fprintf(fTyz,"%0.9g %0.9g\n", localX, tempMesh->getTau()*tempMesh->getPixy(2,3)/P);
			
			fprintf(fTIS, "%0.9g %0.9g\n", localX, tempMesh->getTIS());
			fprintf(fTISB,"%0.9g %0.9g\n", localX, tempMesh->getTISB());
			
			if (mPureBjorken) 
				fprintf(fUz,"%0.9g %0.9g\n", tempMesh->getX(2), tempMesh->getS(3));
		}
	
	int iPrintN = lPrintN, iPrintX = nPrintX;
	for (int i=0; i<tempMesh->getYSize(); i++) {
		if (!tempMesh->getActive(iPrintN,iPrintX,i)) break;
		tempMesh->selfUpdate(iPrintN,iPrintX,i);
		double localY = tempMesh->getX(2);
		
		fprintf(fEY, "%0.9g %0.9g\n", localY, tempMesh->getE());
		fprintf(fUy, "%0.9g %0.9g\n", localY, tempMesh->getS(2));
		fprintf(fUx, "%0.9g %0.9g\n", localY, tempMesh->getS(1));
		
		fprintf(fA1y,"%0.9g %0.9g\n", localY, tempMesh->getS(5));
		fprintf(fA2y,"%0.9g %0.9g\n", localY, tempMesh->getS(6));
		fprintf(fA3y,"%0.9g %0.9g\n", localY, tempMesh->getS(7));
		fprintf(fA4y,"%0.9g %0.9g\n", localY, tempMesh->getS(8));
		fprintf(fA5y,"%0.9g %0.9g\n", localY, tempMesh->getS(9));
		
		fprintf(fTxy,"%0.9g %0.9g\n", localY, tempMesh->getPixy(1,2)/tempMesh->getP());
		
		fprintf(fBy, "%0.9g %0.9g\n", localY, tempMesh->getS(10));
	}
	
	if (mOctant==0)
		for (int i=0; i<tempMesh->getYSize(); i++) {
		if (!tempMesh->getActive(iPrintN,iPrintX,-i)) break;
		tempMesh->selfUpdate(iPrintN,iPrintX,-i);
		double localY = tempMesh->getX(2);
		
		fprintf(fEY, "%0.9g %0.9g\n", localY, tempMesh->getE());
		fprintf(fUy, "%0.9g %0.9g\n", localY, tempMesh->getS(2));
		fprintf(fUx, "%0.9g %0.9g\n", localY, tempMesh->getS(1));
		
		fprintf(fA1y,"%0.9g %0.9g\n", localY, tempMesh->getS(5));
		fprintf(fA2y,"%0.9g %0.9g\n", localY, tempMesh->getS(6));
		fprintf(fA3y,"%0.9g %0.9g\n", localY, tempMesh->getS(7));
		fprintf(fA4y,"%0.9g %0.9g\n", localY, tempMesh->getS(8));
		fprintf(fA5y,"%0.9g %0.9g\n", localY, tempMesh->getS(9));
		
		fprintf(fTxy,"%0.9g %0.9g\n", localY, tempMesh->getPixy(1,2)/tempMesh->getP());
		
		fprintf(fBy, "%0.9g %0.9g\n", localY, tempMesh->getS(10));
	}
	
	for (int i=0;i < min(tempMesh->getYSize(),tempMesh->getXSize());i++) {
		if (!tempMesh->getActive(lPrintN,i,i)) break;
		double P = tempMesh->getP(lPrintN,i,i);
		double localR = sqrt( tempMesh->getX(lPrintN,i,i,1)*tempMesh->getX(lPrintN,i,i,1) 
							 + tempMesh->getX(lPrintN,i,i,2)*tempMesh->getX(lPrintN,i,i,2));
		
		fprintf(fTrr,"%0.9g %0.9g\n", localR, tempMesh->getTRR(lPrintN,i,i)/P);
		fprintf(fTpp,"%0.9g %0.9g\n", localR, tempMesh->getTPhiPhi(lPrintN,i,i)/P);
		fprintf(fTpr,"%0.9g %0.9g\n", localR, tempMesh->getTRPhi(lPrintN,i,i)/P);
	}
	
	/*
	if (!mPureBjorken)
		for (int i=0; i<tempMesh->getXSize(); i++) 
			for (int j=0;j<tempMesh->getYSize(); j++) {
				if (tempMesh->getActive(0,i,j)){
					tempMesh->selfUpdate(0,i,j);
					fprintf(fP,"%f %f %f %f\n",tempMesh->getX(1),tempMesh->getX(2),
							tempMesh->getS(1),tempMesh->getS(2));
				}
				if (tempMesh->getActive(4,i,j)){
					tempMesh->selfUpdate(4,i,j);
					fprintf(fS,"%f %f %f %f\n",tempMesh->getX(1),tempMesh->getX(2),
							tempMesh->getS(1),tempMesh->getS(2));
				}
				if (tempMesh->getActive(12,i,j)){
				tempMesh->selfUpdate(12,i,j);
					fprintf(fT,"%f %f %f %f\n",tempMesh->getX(1),tempMesh->getX(2),
							tempMesh->getS(1),tempMesh->getS(2));
				}
			}
	 */
	 
	closeSliceFiles();
}

void CHydro::printIntegrals() {

	double mTau = tempMesh->getTau();
	tempMesh->calcIntegrals();
	
	offMesh->calcLossIntegrals();
	
	double intS = tempMesh->integralS();
	double intE = tempMesh->integralE();
	
	if (mRK2) 
		if (mLogSinhT) {
			eLoss  += mDt * tanh(offMesh->getTau()) * offMesh->getELoss();
			sLoss  += mDt * tanh(offMesh->getTau()) * offMesh->getSLoss();
		} 
		else if (mLogT) {
			eLoss  += mDt * offMesh->getTau() * offMesh->getELoss();
			sLoss  += mDt * offMesh->getTau() * offMesh->getSLoss();
		}
		else {
			eLoss  += (tempMesh->getTau() - onMesh->getTau()) * offMesh->getELoss();
			sLoss  += (tempMesh->getTau() - onMesh->getTau()) * offMesh->getSLoss();
		}

	fprintf(fIE,"%0.9g %0.9g\n", mTau, (intE+eLoss)/intE0 - 1.);
	fprintf(fIS,"%0.9g %0.9g\n", mTau, (intS+sLoss)/intS0 - 1.);
	fprintf(fSLoss,"%0.6g %0.6g\n", mTau, sLoss/intS0);
	fprintf(fELoss,"%0.6g %0.6g\n", mTau, eLoss/intE0);
 
	fprintf(fIEP,"%0.6g %0.6g\n", mTau - mT0, tempMesh->integralEp());
	fprintf(fIEX,"%0.6g %0.6g\n", mTau - mT0, tempMesh->integralEx());
	fprintf(fIVR,"%0.6g %0.6g\n", mTau - mT0, tempMesh->averageVr());
}

void CHydro::printSpots() {
	
	double mTau = tempMesh->getTau();
	
	fprintf(fE0,"%0.9g %0.9g\n", mTau, tempMesh->getE(0,0,0));
	fprintf(fT0,"%0.6g %0.6g\n", mTau, tempMesh->getT(0,0,0));
	
	fprintf(fTxx0,"%0.6g %0.6g\n", mTau, tempMesh->getTxx(0,0,0));
		//	fprintf(fTzz0,"%0.6g %0.6g\n", mTau, tempMesh->getTzz(0,0,0));
	fprintf(fTzz0,"%0.9g %0.9g\n",mTau,tempMesh->getPixy(0,0,0,0,0));
	
	if (mPureBjorken) {
		/*
		fprintf(fUL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(0,tempMesh->getXSize()-1,0,1));
		fprintf(fEL,"%0.6g %0.6g\n", mTau,  tempMesh->getE(0,tempMesh->getXSize()-1,0));
		fprintf(fA1L,"%0.6g %0.6g\n", mTau, tempMesh->getS(0,tempMesh->getXSize()-1,0,5));
		fprintf(fA2L,"%0.6g %0.6g\n", mTau, tempMesh->getS(0,tempMesh->getXSize()-1,0,6));
		fprintf(fBL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(0,tempMesh->getXSize()-1,0,4));
		*/
		
			//		fprintf(fUL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(0,0,tempMesh->getYSize()-1,2));
			//fprintf(fEL,"%0.6g %0.6g\n", mTau,  tempMesh->getE(0,0,tempMesh->getYSize()-1));
		fprintf(fA1L,"%0.6g %0.6g\n", mTau, tempMesh->getS(0,0,tempMesh->getYSize()-1,5));
		fprintf(fA2L,"%0.6g %0.6g\n", mTau, tempMesh->getS(0,0,tempMesh->getYSize()-1,6));
			//fprintf(fBL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(0,0,tempMesh->getYSize()-1,10));
		
		fprintf(fEL,"%0.6g %0.6g\n", mTau,  tempMesh->getE(0,0,0));
		fprintf(fBL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(0,0,0,10));
	} 
	else {
		fprintf(fUL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(tempMesh->getNSize()-1,0,0,3));
		fprintf(fEL,"%0.6g %0.6g\n", mTau,  tempMesh->getE(tempMesh->getNSize()-1,0,0)/mE0);
		fprintf(fA1L,"%0.6g %0.6g\n", mTau, tempMesh->getS(tempMesh->getNSize()-1,0,0,5));
		fprintf(fA2L,"%0.6g %0.6g\n", mTau, tempMesh->getS(tempMesh->getNSize()-1,0,0,6));
		fprintf(fBL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(tempMesh->getNSize()-1,0,0,10));
	}

	fprintf(fFOSL,"%0.6g %0.6g %0.6g\n", tempMesh->getFOSX(), tempMesh->getFOSY(), mTau);
	fprintf(fFOSSST,"%0.6g %0.6g\n", mTau, tempMesh->getFOSSST());
	fprintf(fUL,"%0.6g %0.6g\n", mTau, tempMesh->getFOSV());
	fprintf(fTxxL,"%0.6g %0.6g\n", mTau, tempMesh->getPixy(0,mPrintX,mPrintY,1,1));
	fprintf(fTxxNSL,"%0.6g %0.6g\n", mTau, tempMesh->getPixyNS(0,mPrintX,mPrintY,1,1));
	fprintf(fTxxML,"%0.6g %0.6g\n", mTau, tempMesh->getPixy(0,mPrintX,mPrintY,1,1));
	fprintf(fDzz,"%0.6g %0.6g\n",mTau, mTau*tempMesh->getDULocal(0,0,0,3,3)-1.);
}

void CHydro::printMesh() {
	
#ifndef HYDRO_XDMF_OUTPUT
	string temp = mDataRoot;
	char num[2];
	sprintf(num,"%d",nPrintFull);
	temp.append("fFull");
	temp.append(num);
	temp.append(".txt");
	fFull = fopen(temp.c_str(),"w");
	nPrintFull++;
	
	int mN = tempMesh->getNSize();
	int mX = tempMesh->getXSize();
	int mY = tempMesh->getYSize();	
	fprintf(fFull,"%0.6g %d %d %d \n",tempMesh->getTau(),mN,mX,mY);
	if (!mPureBjorken)
	  for (int i=-1;i<=mN;i++) 
	    for (int j=-1;j<=mX;j++) 
	      for (int k=-1;k<=mY;k++) {
			  /*
			  if (!tempMesh->getActive(i,j,k)) {
				  fprintf(fFull,"0\n");
				  break;
			  }
			  fprintf(fFull,"1 ");
			   */
			  if (!tempMesh->getActive(i,j,k)) {
				  fprintf(fFull,"0 ");
				  for (int l=1;l<4;l++)  fprintf(fFull,"%0.6g ",tempMesh->getX(i,j,k,l));
				  for (int l=0;l<11;l++) fprintf(fFull,"0. ");
			  }
			  else {
				  fprintf(fFull,"1 ");
				  for (int l=1;l<4;l++)  fprintf(fFull,"%0.6g ",tempMesh->getX(i,j,k,l));
				  for (int l=0;l<11;l++) fprintf(fFull,"%0.9g ",tempMesh->getS(i,j,k,l));
			  }
			  fprintf(fFull,"\n");
		  }
	else 
	  for (int j=0;j<=mX;j++) 
	      for (int k=0;k<=mY;k++) {
			  if (!tempMesh->getActive(0,j,k)) {
				  fprintf(fFull,"0 \n");
				  break;
			  }
			  fprintf(fFull,"1 ");
			  for (int l=1;l<4;l++)  fprintf(fFull,"%0.6g ",tempMesh->getX(0,j,k,l));
			  for (int l=0;l<11;l++) fprintf(fFull,"%0.9g ",tempMesh->getS(0,j,k,l));
			  fprintf(fFull,"\n");
		  }
	
	fflush(fFull);
	fclose(fFull);
#else
	WriteXdmf();
		//	MoveH5();
#endif
}

void CHydro::printTQM() {
	double mZ = 0.;
	int mX = tempMesh->getXSize();
	int mY = tempMesh->getYSize();	
	for (int j=0;j<mX;j++) {
	  int k=0;
	  tempMesh->selfUpdate(0,j,k);
	  fprintf(fTQMX,"%0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g\n",
		  tempMesh->getTau(), tempMesh->getX(0,j,k,1), tempMesh->getX(0,j,k,2), sqrt(tempMesh->getX(0,j,k,1)*tempMesh->getX(0,j,k,1)+tempMesh->getX(0,j,k,2)*tempMesh->getX(0,j,k,2)), //1-4
		  tempMesh->getS(0,j,k,4), tempMesh->getS(0,j,k,1)/tempMesh->getS(0,j,k,0), tempMesh->getS(0,j,k,2)/tempMesh->getS(0,j,k,0),  //5-7
		  sqrt(tempMesh->getS(0,j,k,1)*tempMesh->getS(0,j,k,1) + tempMesh->getS(0,j,k,2)*tempMesh->getS(0,j,k,2))/tempMesh->getS(0,j,k,0), //8
		  tempMesh->getS(0,j,k,0), tempMesh->getDUDT(0,j,k,0),tempMesh->getTxy(0,j,k,0,0), tempMesh->getTxy(0,j,k,0,1), //9-12
		  tempMesh->getTxy(0,j,k,0,2), tempMesh->getP(0,j,k), tempMesh->getBV(0,j,k), //13-15
		  (tempMesh->getS(0,j,k,1)*tempMesh->getDUDT(0,j,k,0) + tempMesh->getS(0,j,k,2)*tempMesh->getDUDT(0,j,k,1) + tempMesh->getS(0,j,k,3)*tempMesh->getDUDT(0,j,k,2))/tempMesh->getS(0,j,k,0) 
		  + tempMesh->getDS(0,j,k,0,0) + tempMesh->getDS(0,j,k,1,1) + tempMesh->getDS(0,j,k,2,2), //16
		  tempMesh->getPixy(0,j,k,0,0), tempMesh->getPixy(0,j,k,0,1), tempMesh->getPixy(0,j,k,0,2), tempMesh->getPixy(0,j,k,3,3), //17-20
		  tempMesh->getPixy(0,j,k,1,1), tempMesh->getPixy(0,j,k,2,2), tempMesh->getPixy(0,j,k,1,2), mZ, mZ, mZ, //21-26
		  tempMesh->getPiZR(0,j,k), tempMesh->getPiZPhi(0,j,k), tempMesh->getPiRR(0,j,k), tempMesh->getPiPhiPhi(0,j,k), tempMesh->getPiRPhi(0,j,k), //27-31
		  tempMesh->getT(0,j,k));
	}
	fprintf(fTQMX,"\n");
	
	for (int k=0;k<mY;k++) {
	  int j=0;
	  tempMesh->selfUpdate(0,j,k);
	  fprintf(fTQMY,"%0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g\n",
		  tempMesh->getTau(), tempMesh->getX(0,j,k,1), tempMesh->getX(0,j,k,2), sqrt(tempMesh->getX(0,j,k,1)*tempMesh->getX(0,j,k,1)+tempMesh->getX(0,j,k,2)*tempMesh->getX(0,j,k,2)), //1-4
		  tempMesh->getS(0,j,k,4), tempMesh->getS(0,j,k,1)/tempMesh->getS(0,j,k,0), tempMesh->getS(0,j,k,2)/tempMesh->getS(0,j,k,0),  //5-7
		  sqrt(tempMesh->getS(0,j,k,1)*tempMesh->getS(0,j,k,1) + tempMesh->getS(0,j,k,2)*tempMesh->getS(0,j,k,2))/tempMesh->getS(0,j,k,0), //8
		  tempMesh->getS(0,j,k,0), tempMesh->getDUDT(0,j,k,0),tempMesh->getTxy(0,j,k,0,0), tempMesh->getTxy(0,j,k,0,1), //9-12
		  tempMesh->getTxy(0,j,k,0,2), tempMesh->getP(0,j,k), tempMesh->getBV(0,j,k), tempMesh->getDUDT(0,j,k,0) + tempMesh->getDS(0,j,k,0,0) + tempMesh->getDS(0,j,k,1,1) + tempMesh->getDS(0,j,k,2,2), //13-16
		  tempMesh->getPixy(0,j,k,0,0), tempMesh->getPixy(0,j,k,0,1), tempMesh->getPixy(0,j,k,0,2), tempMesh->getPixy(0,j,k,2,2), //17-20
		  tempMesh->getPixy(0,j,k,1,1), tempMesh->getPixy(0,j,k,2,2), tempMesh->getPixy(0,j,k,1,2), mZ, mZ, mZ, //21-26
		  tempMesh->getPiZR(0,j,k), tempMesh->getPiZPhi(0,j,k), tempMesh->getPiRR(0,j,k), tempMesh->getPiPhiPhi(0,j,k), tempMesh->getPiRPhi(0,j,k), //27-31
		  tempMesh->getT(0,j,k));
	}
	fprintf(fTQMY,"\n");
	
	for (int j=0;j<min(mX,mY);j++) {
	  int k=j;
	  tempMesh->selfUpdate(0,j,k);
	  fprintf(fTQMR,"%0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g\n",
		  tempMesh->getTau(), tempMesh->getX(0,j,k,1), tempMesh->getX(0,j,k,2), sqrt(tempMesh->getX(0,j,k,1)*tempMesh->getX(0,j,k,1)+tempMesh->getX(0,j,k,2)*tempMesh->getX(0,j,k,2)), //1-4
		  tempMesh->getS(0,j,k,4), tempMesh->getS(0,j,k,1)/tempMesh->getS(0,j,k,0), tempMesh->getS(0,j,k,2)/tempMesh->getS(0,j,k,0),  //5-7
		  sqrt(tempMesh->getS(0,j,k,1)*tempMesh->getS(0,j,k,1) + tempMesh->getS(0,j,k,2)*tempMesh->getS(0,j,k,2))/tempMesh->getS(0,j,k,0), //8
		  tempMesh->getS(0,j,k,0), tempMesh->getDUDT(0,j,k,0),tempMesh->getTxy(0,j,k,0,0), tempMesh->getTxy(0,j,k,0,1), //9-12
		  tempMesh->getTxy(0,j,k,0,2), tempMesh->getP(0,j,k), tempMesh->getBV(0,j,k), tempMesh->getDUDT(0,j,k,0) + tempMesh->getDS(0,j,k,0,0) + tempMesh->getDS(0,j,k,1,1) + tempMesh->getDS(0,j,k,2,2), //13-16
		  tempMesh->getPixy(0,j,k,0,0), tempMesh->getPixy(0,j,k,0,1), tempMesh->getPixy(0,j,k,0,2), tempMesh->getPixy(0,j,k,2,2), //17-20
		  tempMesh->getPixy(0,j,k,1,1), tempMesh->getPixy(0,j,k,2,2), tempMesh->getPixy(0,j,k,1,2), mZ, mZ, mZ, //21-26
		  tempMesh->getPiZR(0,j,k), tempMesh->getPiZPhi(0,j,k), tempMesh->getPiRR(0,j,k), tempMesh->getPiPhiPhi(0,j,k), tempMesh->getPiRPhi(0,j,k), //27-31
		  tempMesh->getT(0,j,k));
	}
	fprintf(fTQMR,"\n");
}

// FIXME usage of CMesh::getS(int, int, int, 0) is antiquated.
void CHydro::printOscarFull(int mT) {
	int mX = tempMesh->getXSize();
	int mY = tempMesh->getYSize();
	int mN = tempMesh->getNSize();
	
	double mSVRatio = parameter::getD(*pMap,"HYDRO_SVRATIO",0.0);
	double mBVRatio = parameter::getD(*pMap,"HYDRO_BVRATIO",0.0);
	
	if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false))
		for (int i=0;i<mX;i++)
			for (int j=0;j<mY;j++) {
				if (!tempMesh->getActive(0,i,j)) break;
				tempMesh->selfUpdate(0,i,j);
				fprintf(fOscarFull,"%d %d %d ",mT,i,j);
				fprintf(fOscarFull,"%f %f %f 1. ",tempMesh->getE(0,i,j),tempMesh->getP(0,i,j),tempMesh->getT(0,i,j));
				
				fprintf(fOscarFull,"%f %f ",tempMesh->getS(0,i,j,1)/tempMesh->getS(0,i,j,0), tempMesh->getS(0,i,j,2)/tempMesh->getS(0,i,j,0));
				
				if (mSVRatio > 0.) {
					double mAlpha = tempMesh->getISAlpha(0,i,j);
					for (int m=1;m<6;m++) 
						fprintf(fOscarFull,"%f ",mAlpha*tempMesh->getS(0,i,j,m+4));
				}
				if (mBVRatio > 0.) 
					fprintf(fOscarFull,"%f ",tempMesh->getS(0,i,j,10)*tempMesh->getISGamma(0,i,j));
				
				if (mSVRatio > 0.)  
					fprintf(fOscarFull,"%f %f",tempMesh->getSV(0,i,j),tempMesh->getTIS(0,i,j));
				if (mBVRatio > 0.)
					fprintf(fOscarFull,"%f %f",tempMesh->getBV(0,i,j),tempMesh->getTISB(0,i,j));
				
				fprintf(fOscarFull,"%f",tempMesh->getS(0,i,j));
				fprintf(fOscarFull,"\n");
			}
	else 
		for (int i=0;i<mX;i++)
			for (int j=0;j<mY;j++) 
				for (int k=0;k<mN;k++) {
					if (!tempMesh->getActive(k,i,j)) break;
					tempMesh->selfUpdate(k,i,j);
					fprintf(fOscarFull,"%d %d %d %d ",mT,i,j,k);
					fprintf(fOscarFull,"%g %g %g 1. ",tempMesh->getE(k,i,j),tempMesh->getP(k,i,j),tempMesh->getT(k,i,j));
					
					fprintf(fOscarFull,"%g %g %g ",tempMesh->getS(k,i,j,1)/tempMesh->getS(k,i,j,0),
							tempMesh->getS(k,i,j,2)/tempMesh->getS(k,i,j,0), 
							tempMesh->getS(k,i,j,3)/tempMesh->getS(k,i,j,0));
					
					if (mSVRatio > 0.) {
						double mAlpha = tempMesh->getISAlpha(k,i,j);
						for (int m=1;m<6;m++) 
							fprintf(fOscarFull,"%g ",mAlpha*tempMesh->getS(k,i,j,m+4));
					}
					if (mBVRatio > 0.) 
						fprintf(fOscarFull,"%g ",tempMesh->getS(k,i,j,10)*tempMesh->getISGamma(k,i,j));
					
					if (mSVRatio > 0.)  
						fprintf(fOscarFull,"%g %g ",tempMesh->getSV(k,i,j),tempMesh->getTIS(k,i,j));
					if (mBVRatio > 0.)
						fprintf(fOscarFull,"%g %g ",tempMesh->getBV(k,i,j),tempMesh->getTISB(k,i,j));
					
					fprintf(fOscarFull,"%g ",tempMesh->getS(k,i,j));
					fprintf(fOscarFull,"\n");
					fflush(fOscarFull);
				}
}

bool CHydro::printOscarHyper(int mT) {
	
		// it ix [iy iz] e p T R_qgp vx [vy y_L] [n(1)...n(C)] [mu(1)...mu(C)] 
		//dsig_t dsig_x [dsig_y dsig_eta] [Diss(1)...Diss(D)]  [Tr(1)...Tr(T)]
	double mFOS[XSIZE+YSIZE][2];
	double mFOSigma[XSIZE+YSIZE][2];
	double mFOVelo[XSIZE+YSIZE][2];
	double mFODiss[XSIZE+YSIZE][6];
	int fosSize;
	double temp = mFoTemp;
	double ed = mEos->getEGivenT(temp);
	
	fprintf(fOscarHyper,"time: %g\n",tempMesh->getTau());
	
	if (!mPureBjorken) {
		for (int k=0;k<tempMesh->getNSize();k++){
			if(!tempMesh->getFOS(k,mFOS,mFOSigma,mFOVelo,mFODiss,fosSize)){
				fflush(fOscarHyper);
				break;
			}
				
			fprintf(fOscarHyper,"eta: %g size: %d\n",tempMesh->getX(k,0,0,3),fosSize);
			for (int i=0;i<fosSize;i++) {
				fprintf(fOscarHyper,"%g %g ", mFOS[i][0], mFOS[i][1]);
				fprintf(fOscarHyper,"%g %g %g 1. ",ed,mEos->getP(ed),temp);
				for (int k=0;k<2;k++)
					fprintf(fOscarHyper,"%g ",mFOVelo[i][k]);
				fprintf(fOscarHyper,"0. %g %g ", mFOSigma[i][0], mFOSigma[i][1]);
				for (int k=0;k<6;k++)
					fprintf(fOscarHyper,"%g ",mFODiss[i][k]);
				fprintf(fOscarHyper,"\n");
			}
		}
	}
	else {
		if (!tempMesh->getFOS(mFOS,mFOSigma,mFOVelo,mFODiss,fosSize)){
			fflush(fOscarHyper);
			return false;
		}
		for (int i=0;i<fosSize;i++) {
			fprintf(fOscarHyper,"%g %g ", mFOS[i][0], mFOS[i][1]);
			fprintf(fOscarHyper,"%g %g %g 1. ",ed,mEos->getP(ed),temp);
			for (int k=0;k<2;k++)
				fprintf(fOscarHyper,"%g ",mFOVelo[i][k]);
			fprintf(fOscarHyper,"0. %g %g ", mFOSigma[i][0], mFOSigma[i][1]);
			for (int k=0;k<6;k++)
				fprintf(fOscarHyper,"%g ",mFODiss[i][k]);
			fprintf(fOscarHyper,"\n");
		}
	}
	fflush(fOscarHyper);
	return true;
}

bool CHydro::checkDirectory(string fName){
	size_t mSpot = fName.find_last_of('/');
	
	string fDir;
	for (int i=0;i<=mSpot;i++)
		fDir[i]=fName[i];
		
	struct stat st;
	if (stat(fDir.c_str(),&st) != 0)
		return false;
	else 
		return true;
}

bool CHydro::checkDirectoryExistence(string dirname){
	struct stat dirinfo;
	int intstat=stat(dirname.c_str(),&dirinfo);
	if(intstat==0) return true;
	else return false;
}

void CHydro::printDNs(CMesh* tempMesh, int mNSize, int mXSize, int mYSize) {
}

void CHydro::testFileOpen() {
	return;
	if (fEX != NULL) printf("problem with file fEX (%p)...\n",fEX);
	if (fEN != NULL) printf("problem with file fEN (%p)...\n",fEN);
	if (fE0 != NULL) printf("problem with file fE0 (%p)...\n",fE0);
	if (fIE != NULL) printf("problem with file fIE (%p)...\n",fIE);
	if (fA1x != NULL) printf("problem with file fA1X (%p)...\n",fA1x);
	if (fTxx != NULL) printf("problem with file fTxx (%p)...\n",fTxx);
	if (fTrr != NULL) printf("problem with file fTrr (%p)...\n",fTrr);
	if (fSLoss != NULL) printf("problem with file fSloss (%p)...\n",fSLoss);
	if (fUL != NULL) printf("problem with file fUL (%p)...\n",fUL);
	if (fTIS != NULL) printf("problem with file fTIS (%p)...\n",fTIS);
	if (fIVR != NULL) printf("problem with file fVR0 (%p)..\n",fIVR);
	if (fFOS != NULL) printf("problem with file fFOS (%p)...\n",fFOS);
	if (fTQMX != NULL) printf("problem with file fTQMX (%p)...\n",fTQMX);
	if (fTxxL != NULL) printf("problem with file fTxxL (%p)...\n",fTxxL);
	if (fFull != NULL) printf("problem with file fFull (%p)...\n",fFull);
	if (fPiDNDY != NULL) printf("problem with file fPiDNDYU (%p)...\n",fPiDNDY);
}

void CHydro::printVorticity() {
	string fn = "fVorticity.txt";
	fn = mDataRoot + fn;
	FILE* f = fopen( fn.c_str(),"w");

	CMesh* m = tempMesh;
	
	if (!mPureBjorken)
		for (int i=0;i<=tempMesh->getNSize();i++) 
			for (int j=0;j<=tempMesh->getXSize();j++)
				for (int k=0;k<=tempMesh->getYSize();k++)
					if (!m->getActive(i,j,k))
						break;
					else {
						m->update(i,j,k);
						fprintf(f,"%g %g %g %g %g %g\n",m->getX(i,j,k,1),m->getX(i,j,k,2),m->getX(i,j,k,3),
								m->getVorticity(i,j,k,1,2),m->getVorticity(i,j,k,1,3),m->getVorticity(i,j,k,2,3));
					}
	else 
		for (int j=0;j<=tempMesh->getXSize();j++)
			for (int k=0;k<=tempMesh->getYSize();k++) 
				if (m->getActive(0,j,k)) {
					m->update(0,j,k);
					fprintf(f,"%g %g %g %g %g %g\n",m->getX(0,j,k,1),m->getX(0,j,k,2),m->getE(0,j,k),
							m->getVorticity(0,j,k,1,2),m->getS(0,j,k,1),m->getS(0,j,k,2));
						//					fprintf(f,"%g %g 0. %g %g %g\n",m->getX(0,j,k,1),m->getX(0,j,k,2),
						//							m->getVorticity(0,j,k,1,2),m->getVorticity(0,j,k,1,3),m->getVorticity(0,j,k,2,3));
				}
				else break;
	
	fclose(f);
}

void CHydro::printActive() {
	FILE *fActive  = fopen("activeCells.txt","w");
	FILE *fActive2 = fopen("activeCells2.txt","w");
	
		// set loop ranges
	int lowX1=0, highX1=tempMesh->getXSize(),lowY1=0, highY1=tempMesh->getYSize();
	int lowX2=0, highX2=  onMesh->getXSize(),lowY2=0, highY2=  onMesh->getYSize();
	if (mOctant != 3){
		lowX1 = -highX1; lowY1 = -highY1;
		lowX2 = -highX2; lowY2 = -highY2;
	}

	for (int j=lowX1;j<=highX1;j++)
		for (int k=lowY1;k<=highY1;k++) 
			if (tempMesh->getActive(0,j,k))
				fprintf(fActive,"%g %g\n",tempMesh->getX(0,j,k,1),tempMesh->getX(0,j,k,2));
	
	for (int j=lowX2;j<=highX2;j++)
		for (int k=lowY2;k<=highY2;k++) 
			if (onMesh->getActive(0,j,k))
				fprintf(fActive2,"%g %g\n",onMesh->getX(0,j,k,1),onMesh->getX(0,j,k,2));
	
	fclose(fActive);
	fclose(fActive2);
}

#ifdef HYDRO_XDMF_OUTPUT
int CHydro::ReadDataH5(){
	
}

void CHydro::fillCellsH5(CCellH5& cH, CCell* p) {
		//	p->update();
	
	for (int i=0;i<4;i++)
		cH.x[i] = p->getX(i);
	for (int i=0;i<3;i++)
		cH.v[i] = p->getS(i+1);
	for (int i=0;i<5;i++)
		cH.a[i] = p->getA(i+1);
	cH.e = p->getE();
	cH.P = p->getPCalc();
	cH.T = p->getTCalc();
	cH.b = p->getB();
	cH.SV = p->getSVCalc();
	cH.tIS = p->getTISCalc();
	cH.BV = p->getBVCalc();
	cH.tISB = p->getTISBCalc();
	cH.S = p->getSCalc();
	
	cH.fQGP = 1.0;
}

int CHydro::WriteDataH5(){
	int size = 0;
	
	int mN = tempMesh->getNSize();
	int mX = tempMesh->getXSize();
	int mY = tempMesh->getYSize();
	
	if (!mPureBjorken)
		for (int i=0;i<mN;i++)
			for (int j=0;j<=mX;j++)
				for (int k=0;k<=mY;k++) 
					if (tempMesh->getActive(i,j,k))
						size++;	
					else break;
	
	CCellH5* cells = new CCellH5[size];
	hsize_t dim[] = {size};
	DataSpace space(1, dim);
	
	int lSize = 0;
	if (!mPureBjorken)
		for (int i=0;i<=mN;i++)
			for (int j=0;j<=mX;j++)
				for (int k=0;k<=mY;k++) 
					if (tempMesh->getActive(i,j,k)){
						fillCellsH5(cells[lSize],tempMesh->getCell(i,j,k));
						lSize++;
					}
					else break;
	
	
	H5File *h5outfile = new H5File("test.h5",H5F_ACC_TRUNC);
	CompType *ptype;
	
	ptype=new CompType(sizeof(CCellH5));
	
	hsize_t* xDim = new hsize_t[1];
	xDim[0] = 4;
	ArrayType hX(PredType::NATIVE_DOUBLE,1,xDim);
	ptype->insertMember("x", HOFFSET(CCellH5,x), hX);
	
	hsize_t* vDim = new hsize_t[1];
	vDim[0] = 3;
	ArrayType hV(PredType::NATIVE_DOUBLE,1,vDim);
	ptype->insertMember("v", HOFFSET(CCellH5,v), hV);
	
	hsize_t* aDim = new hsize_t[1];
	aDim[0] = 5;
	ArrayType hA(PredType::NATIVE_DOUBLE,1,aDim);
	ptype->insertMember("a", HOFFSET(CCellH5,a), hA);
	
	ptype->insertMember("e", HOFFSET(CCellH5,e), PredType::NATIVE_DOUBLE);
	ptype->insertMember("P", HOFFSET(CCellH5,P), PredType::NATIVE_DOUBLE);
	ptype->insertMember("T", HOFFSET(CCellH5,T), PredType::NATIVE_DOUBLE);
	ptype->insertMember("fQGP", HOFFSET(CCellH5,fQGP), PredType::NATIVE_DOUBLE);
	ptype->insertMember("b", HOFFSET(CCellH5,b), PredType::NATIVE_DOUBLE);
	ptype->insertMember("SV", HOFFSET(CCellH5,SV), PredType::NATIVE_DOUBLE);
	ptype->insertMember("tIS", HOFFSET(CCellH5,tIS), PredType::NATIVE_DOUBLE);
	ptype->insertMember("BV", HOFFSET(CCellH5,BV), PredType::NATIVE_DOUBLE);
	ptype->insertMember("tISB", HOFFSET(CCellH5,tISB), PredType::NATIVE_DOUBLE);
	ptype->insertMember("S", HOFFSET(CCellH5,S), PredType::NATIVE_DOUBLE); 
	
	DataSet * dataset;
	dataset = new DataSet(h5outfile->createDataSet("lilTest",*ptype,space));
	dataset->write(cells,*ptype);
	
	delete dataset;
	delete cells;
}

void CHydro::StartXdmf() {
	d = new XdmfDOM();
	root = new XdmfRoot();
	domain = new XdmfDomain();
	tGrid = new XdmfGrid();
	
	root->SetDOM(d);
	root->SetVersion(2.2); // Change the Version number because we can
	root->Build();
	
	root->Insert(domain);
	
	tGrid->SetName("Temporal Grid");
	tGrid->SetGridType(XDMF_GRID_COLLECTION);
	tGrid->SetCollectionType(XDMF_GRID_COLLECTION_TEMPORAL);
	domain->Insert(tGrid);
	
	tGrid->Build();
	
	sShape[0] = tempMesh->getNSizeOrig();
	sShape[1] = tempMesh->getXSizeOrig();
	sShape[2] = tempMesh->getYSizeOrig();
	
	for (int i=0;i<3;i++)
		vShape[i] = sShape[i];
	vShape[3] = 3;
	
		// symmetric tensor data shape
	for (int i=0;i<3;i++)
		tShape[i] = sShape[i];
	tShape[3] = 6;
}

void CHydro::WriteXdmf() {
	CMesh* m = tempMesh;
	
	char cTime[2];
	sprintf(cTime,"%d",nPrintFull);

	XdmfGrid grid;
	tGrid->Insert(&grid);
	
		// Topology
	XdmfTopology* topo = grid.GetTopology();
	topo->SetTopologyType(XDMF_3DCORECTMESH);
	topo->GetShapeDesc()->SetShape(3,sShape);
	
		// Geometry
	XdmfGeometry *geo = grid.GetGeometry();
	geo->SetGeometryType(XDMF_GEOMETRY_ORIGIN_DXDYDZ);
	geo->SetOrigin(0, 0, 0);
	geo->SetDxDyDz(parameter::getD(*pMap,"HYDRO_DX",0.1),
				   parameter::getD(*pMap,"HYDRO_DY",0.1),
				   parameter::getD(*pMap,"HYDRO_DN",0.1));
	
	XdmfTime* xTime = grid.GetTime();
	xTime->SetTimeType(XDMF_TIME_SINGLE);
	xTime->SetValue(m->getTau());
	grid.Insert(xTime);
	
	XdmfAttribute aE, aP, aT, aV, aTxy, aB, aSV, aSTIS, aBV, aBTIS, aS;
	XdmfArray * bTIS =  xdmfScalar(aBTIS,"Bulk Relaxation Time",sShape);
	XdmfArray * bv = xdmfScalar(aBV,"Bulk Viscosity",sShape);
	XdmfArray * e = xdmfScalar(aE,"Energy Density",sShape);
	XdmfArray * s = xdmfScalar(aS,"Entropy Density",sShape);
	XdmfArray * b = xdmfScalar(aB,"ISBulk",sShape);
	XdmfArray * p = xdmfScalar(aP,"Pressure",sShape);
	XdmfArray * sv = xdmfScalar(aSV,"Shear Viscosity",sShape);
	XdmfArray * sTIS = xdmfScalar(aSTIS,"Shear Relaxation Time",sShape);
	XdmfArray * t = xdmfScalar(aT,"Temperature",sShape);
	
	aV.SetName("Velocity");
	aV.SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_CELL);
	aV.SetAttributeType(XDMF_ATTRIBUTE_TYPE_VECTOR);
	
	XdmfArray *v = aV.GetValues();
	v->SetShape(4,vShape);
	
	aTxy.SetName("Shear Tensor");
	aTxy.SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_CELL);
	aTxy.SetAttributeType(XDMF_ATTRIBUTE_TYPE_TENSOR6);
	
	XdmfArray * txy = aTxy.GetValues();
	txy->SetShape(4,tShape);
	
	for (XdmfInt64 k=0;k<sShape[0];k++)
		for (XdmfInt64 j=0;j<sShape[1];j++)
			for (XdmfInt64 i=0;i<sShape[2];i++) {
				XdmfInt64 index = i+sShape[2]*(j+k*sShape[1]);
				m->selfUpdate(k,j,i);
				
				e->SetValue(index, m->getE());
				p->SetValue(index, m->getP());
				t->SetValue(index, m->getT());
				
				for (int n=0;n<3;n++)
					v->SetValue(3*index+n, m->getS(n+1));
				
				txy->SetValue(6*index,   m->getPixy(1,1));
				txy->SetValue(6*index+1, m->getPixy(1,2));
				txy->SetValue(6*index+2, m->getPixy(1,3));
				txy->SetValue(6*index+3, m->getPixy(2,2));
				txy->SetValue(6*index+4, m->getPixy(2,3));
				txy->SetValue(6*index+5, m->getPixy(3,3));
				
				b->SetValue(index, m->getB());
				sv->SetValue(index, m->getSV());
				sTIS->SetValue(index, m->getTIS());
				bv->SetValue(index, m->getBV());
				bTIS->SetValue(index, m->getTISB());
				s->SetValue(index, m->getS());
			}
	
	grid.Insert(&aE);
	grid.Insert(&aP);
	grid.Insert(&aT);
	grid.Insert(&aV);
	grid.Insert(&aTxy);
	grid.Insert(&aB);
	grid.Insert(&aSV);
	grid.Insert(&aSTIS);
	grid.Insert(&aBV);
	grid.Insert(&aBTIS);
	grid.Insert(&aS);
	
	d->Write(XdmfConstString("timeSeries.xmf"));
	
		// Update XML and Write Values to DataItems
	grid.Build();
}

void CHydro::StopXdmf() {
	/*
	string mFN = mDataRoot;
	mFN.append("timeSeries.xmf");
	XdmfConstString localFN(mFN.c_str());
	d->Write(localFN); // write to file
	*/
	
	d->Write(XdmfConstString("timeSeries.xmf")); 
	
	string comm("mv timeSeries.* ");
	comm.append(mDataRoot);
		//	printf("%s\n",comm.c_str());
	system(comm.c_str());
}

XdmfArray* CHydro::xdmfScalar(XdmfAttribute& mAtt, const char name[], XdmfInt64 sShape[3]) {
	mAtt.SetName(name);
	mAtt.SetAttributeCenter(XDMF_ATTRIBUTE_CENTER_CELL);
	mAtt.SetAttributeType(XDMF_ATTRIBUTE_TYPE_SCALAR);
	
	XdmfArray * mArray = mAtt.GetValues();
	mArray->SetShape(3,sShape);
	return mArray;
}

void CHydro::MoveH5() {
		// file location for .xmf
	string mX = mDataRoot;
	mX.append("t");
	
		//	char cTime[7];
		//	sprintf(cTime,"%0.6g",tempMesh->getTau());
	
	char cTime[2];
	sprintf(cTime,"%d",nPrintFull);
	nPrintFull++;
	mX.append(cTime);
	mX.append(".xmf");
	
	
		// file location for .h5
	string mH; // = mDataRoot;
	mH.append("t");
	mH.append(cTime);
	mH.append(".h5");
	
		//sed -i 's/Xdmf.h5/${mH}/' ${mX} > /tmp/temp.xmf
	string sedComm("sed \'s/Xdmf.h5/");
	sedComm.append(mH);
	sedComm.append("/\' ");
	sedComm.append(mX);
	sedComm.append(" > /tmp/temp.xmf");
	system(sedComm.c_str());
	
		//	mv /tmp/tempfile.tmp ${mX}
	string mvCommX("mv /tmp/temp.xmf ");
	mvCommX.append(mX);
	system(mvCommX.c_str());
	
		//	mv /tmp/tempfile.tmp ${mX}
	string mvCommH("mv Xdmf.h5 ");
	mvCommH.append(mDataRoot);
	mvCommH.append(mH);
	system(mvCommH.c_str());
	
}

void CHydro::ReadXdmf() {

	XdmfDOM *dom = new XdmfDOM();
	dom->Parse("junk.xmf");
	
	XdmfXmlNode ge = dom->FindElementByPath("/Xdmf/Domain/Grid");
	XdmfGrid *grid = new XdmfGrid();
	grid->SetDOM(dom);
	grid->SetElement(ge);
	grid->UpdateInformation();
	grid->Update();
	
		//	grid->AssignAttributeByName("Pressure");
		//	XdmfAttribute* p = grid->GetAssignedAttribute();
	for (int i=0;i<grid->GetNumberOfAttributes();i++){
	
		XdmfAttribute* p = grid->GetAttribute(i);
		cout << "Attribute " << i << " Name = " << p->Get("Name") << endl;
		p->Update();
		
		cout << p->Serialize() << endl;
		return;
		
		XdmfArray *arr = p->GetValues(); 
		
		cout << endl << "Emptying array with..." << endl;
		for (int i=0;i<5;i++)
			for (int j=0;j<5;j++)
				cout << arr->GetValues(5*(j+5.*i),5,1) << endl;
	}
	/*
	XdmfTopology *top = grid->GetTopology();
	top->DebugOn();
	
		//	XdmfArray *conn = top->GetConnectivity();
		//	cout << " Values = " << conn->GetValues() << endl;
	
	XdmfAttribute *p = dom->GetAttribute(ge, "Pressure");
	cout << dom->Serialize(ge);
	/*
	XdmfGeometry *geo = grid->GetGeometry();
	XdmfArray *points = geo->GetPoints();

	cout << "Geo Type = " << geo->GetGeometryTypeAsString() <<
	" # Points = " << geo->GetNumberOfPoints() << endl;
	cout << " Points(0,6) = " << points->GetValues(0,6) << endl;
	*/
}

#endif

string CHydro::mDataRoot, CHydro::mOldFileName;
bool CHydro::mDebug, CHydro::mIoIntegrals, CHydro::mIoSlices;
bool CHydro::mIoSpots, CHydro::mIoTechqm;
bool CHydro::mIoOscarFull, CHydro::mIoOscarHyper, CHydro::mIoFull;
bool CHydro::mIoSpectra, CHydro::mOldFile;
bool CHydro::mSawtooth, CHydro::mRK2, CHydro::mRK4;
bool CHydro::mLinT, CHydro::mLogT, CHydro::mLogSinhT;
bool CHydro::mBjorken, CHydro::mPureBjorken, CHydro::mHalving;
int CHydro::mOctant;
int CHydro::mTStep, CHydro::mXSize, CHydro::mYSize, CHydro::mNSize;
int CHydro::mIoSliceTStep, CHydro::mIoTechqmTStep, CHydro::mIoOscarTStep;
int CHydro::mIoSpotsTStep, CHydro::mIoFullTStep;
int CHydro::mPrintX, CHydro::mPrintY, CHydro::mPrintN;
double CHydro::mInitNS, CHydro::mInitNSTxxP, CHydro::mInitFlow, CHydro::mFoTemp;
double CHydro::mT0, CHydro::mDt, CHydro::mE0;
double CHydro::intS0, CHydro::intE0;
double CHydro::eLoss, CHydro::sLoss;