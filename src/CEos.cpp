/*
 *  CEos.cpp
 *  Created by Joshua Vredevoogd on 3/4/09.
 */

#include "CEos.h"

//constructor - read from EOS interpolation
CEosCalculator::CEosCalculator(parameterMap* pM) {
	lastAccess=10;

	mLatEos  = parameter::getB(*pM,"EQOFST_LATEOS",false);
	mFoTemp = parameter::getD(*pM,"HYDRO_FOTEMP",0.165);
	
	mSVRatio = parameter::getD(*pM,"HYDRO_SVRATIO",0.0);
	svSwitchTemp = parameter::getD(*pM,"EQOFST_SVSWITCHTEMP",0.0);
	
	bRealisticEtaS = parameter::getB(*pM,"EQOFST_REALISTIC_ETAS",false);
	mSvHighT = parameter::getD(*pM,"EQOFST_SV_HIGHT",mFoTemp);
	mSvHighTSlope = parameter::getD(*pM,"EQOFST_SV_HIGHT_SLOPE",2.0);
	
	mBVRatio = parameter::getD(*pM,"HYDRO_BVRATIO",0.0);
	mBVCent = parameter::getD(*pM,"EQOFST_BVCENT",mFoTemp);
	mBVWidth = parameter::getD(*pM,"EQOFST_BVWIDTH",0.015);
	
	if (svSwitchTemp > mSvHighT){
		printf("\n****** CEos requires that the region of high temperature be distinct ******\n");
		printf("****** from the region where the shear viscosity is proportional to energy density. Aborting. ******\n\n");
		exit(1);
	}
		
	
	if (mLatEos) {
		
		const int size=2400;
  
		temp = new double[size];
		ed   = new double[size]; 
		pr   = new double[size]; 
		sd   = new double[size];
		
		t2 = new double[size];
		e2 = new double[size];
		p2 = new double[size];
		s2 = new double[size];
		
		if (parameter::getS(*pM,"EQOFST_LATDATA","null") == string("null")){
			printf("\n** CEos requires parameterMap variable EQOFST_LATDATA if EQOFST_LATEOS == true **\nAborting...\n\n");
			exit(1);
		}
		
		string latFileName = parameter::getS(*pM,"EQOFST_LATDATA","none");
		if (!checkFile(latFileName)) {
			printf("\n** No file found at %s (EQOFST_LATDATA) **\nAborting...\n\n",latFileName.c_str());
			exit(1);
		}

		std::ifstream latFile(latFileName.c_str());
		
		double lowT  = parameter::getD(*pM,"HYDRO_FOTEMP",0.15);
		double highT = parameter::getD(*pM,"EQOFST_HIGH_MERGE_T",0.2);
		
		if (highT < lowT) {
			printf("Inconsistent values of HYDRO_FOTEMP and EQOFST_HIGH_CROSS_T.\n Aborting.\n");
			exit(1);
		}
		
		double lT[size], lE[size], lP[size], lS[size];
		double hT[size], hE[size], hP[size], hS[size];
		int lSize=0;
		int hSize=0;
		
		char test;
		latFile >> test;
		latFile.ignore(1000,'\n');
		
			// re-used counter
		int i;
		
			// check first line to determine file format
		if (test == 'e') { // EOS from Pasi
			double junk;
			for(i=0; !latFile.eof() && i<size; i++) {
				latFile >> lE[i] >> lP[i] >> lS[i] >> junk;
				lT[i] = (lE[i] + lP[i])/lS[i];
			}
			lSize = i-1;
			latFile.close();
		}
		else { // EOS from Ron/Fodor
			double junk;
			for(i=0; !latFile.eof() && i<size; i++) {
				latFile >> lT[i] >>  junk >> lE[i] >> lP[i] >> lS[i] >> junk >> junk;
				
				lE[i]   *= pow(lT[i],4.)/pow(JHBARC,3.);
				lP[i]   *= pow(lT[i],4.)/pow(JHBARC,3.);
				lS[i]   *= pow(lT[i],3.)/pow(JHBARC,3.);
			}	
			lSize = i-1;
		}
		
			// get hrg file name from parameterMap
		string hrgFileName = parameter::getS(*pM,"EQOFST_HRGDATA","none");		
		
			// not merging lattice with HRG data
		if (hrgFileName == string("none")){
				// move variables to permanent homes
			aSize = lSize;
			for (int j=0;j<aSize;j++) {
				temp[j] = lT[j]; 
				sd[j] = lS[j];
				ed[j] = lE[j];
				pr[j] = lP[j];
			} 
		}
		else if (!checkFile(hrgFileName)) {
			printf("Unable to find eos data at %s....\n Continuing with lattice data only....\n",
				   hrgFileName.c_str());
			
				// move variables to permanent homes
			aSize = lSize;
			for (int j=0;j<aSize;j++) {
				temp[j] = lT[j]; 
				sd[j] = lS[j];
				ed[j] = lE[j];
				pr[j] = lP[j];
			} 
		}
		else {
				// Open hrg eos file for merging
			std::ifstream hrgFile;
			hrgFile.open(hrgFileName.c_str());
			
				// skip first lines of hrg file
			for (int j=0;j<8;j++) 
				hrgFile.ignore(1000,'\n');
			
				// read hrg data
			for(i=0; !hrgFile.eof() && i<size; i++) {
				hrgFile >> hT[i] >> hS[i] >> hP[i] >> hE[i];
				hT[i] /= 1000.; hP[i] /= 1000.; hE[i] /= 1000.;
				if (hT[i] < 0.02) i--;
			}
			hSize = i-1;	
			
				// merge will fail if 
			if (hT[hSize] < highT) 
				if (hT[hSize] < lowT){
					printf("No hadron resonance gas data up to %g GeV (HYDRO_FOTEMP).\n",lowT);
					printf("Unclear how to proceed.  Aborting.\n");
				}
				else {
					printf("No hadron resonance gas data up to %g GeV (EQOFST_HIGH_CROSS_T).\n",highT);
					printf("Resetting to %g GeV, maximum value in %s (EQOFST_HRGDATA).\n",
						   hT[hSize],parameter::getS(*pM,"EQOFST_LATDATA","none").c_str());
					highT = hT[hSize];
				}
			
				// offset by tenth of step eliminating round-off error effects
			lowT  += (lT[1] - lT[0])/10.;
			highT += (lT[1] - lT[0])/10.;
			
				// find boundaries in each array
			int hlowTIndex=0, llowTIndex=0, highTIndex=0;
			while (hT[hlowTIndex] < lowT ) hlowTIndex++;
			while (lT[llowTIndex] < lowT ) llowTIndex++; 
			while (lT[highTIndex] < highT) highTIndex++;
			hlowTIndex--; llowTIndex--;
			
				// copy hrg data below TLow
			aSize=0;
			for (int j=0; j<=hlowTIndex; j++) {
				sd[j] = hS[j];
				temp[j] = hT[j];
				aSize++;
			}
			
				// merge between TLow and THigh
			for (int j=1;;j++) {
				temp[aSize] = lT[llowTIndex+j];
				if (temp[aSize] >= highT) break;
				double w = (tanh( tan((PI/2.)*( 2.*(highT-temp[aSize])/(highT-lowT) - 1.)))+1.)/2. ;
				sd[aSize] = w*hS[hlowTIndex+j] + (1.-w)*lS[llowTIndex+j];
				aSize++;
			}
			aSize--;

				// take lattice values above THigh
			for (int j=highTIndex-1;j<lSize;j++) {
				sd[aSize] = lS[j];
				temp[aSize] = lT[j];
				aSize++;
			}
			aSize--;
			
				// take hrg pressure or calculate by integral method
			for (int j=0;j<=aSize;j++) 
				if (j<= hlowTIndex)
					pr[j] = hP[j];
				else 
					pr[j] = pr[j-1] + 0.5*(temp[j]-temp[j-1])*(sd[j]+sd[j-1]);
			
				// energy density from therm identitiy
			for (int j=0;j<=aSize;j++) 
				ed[j] = temp[j]*sd[j] - pr[j];
		}		
				
			// add bump to eos for statistical runs
		if ( parameter::getB(*pM,"EQOFST_ADD_BUMP",false) ) {
			
				// temperature to add bump to s/T^3
				//			double T0     = parameter::getD(*pM,"EQOFST_BUMP_START",-1.);
			double T0 = parameter::getD(*pM,"HYDRO_FOTEMP",0.15);
				// width (in temperature [GeV]) of bump
			double width  = parameter::getD(*pM,"EQOFST_BUMP_WIDTH",-1.);
				// height of bump to s/T^3 
				// given as fraction of maximum possible change (for cs2 to stay in bounds)
			double height = parameter::getD(*pM,"EQOFST_BUMP_HEIGHT",-1.); 
			
				//
			bExpTail   = false;
			bGaussTail = false;
			bPowerTail = true; 
			
				// do nothing if any value is not set
			if (height == -1. || T0 == -1. || width == -1.) {
					// print an error message if they have one but not all
				if (height != -1. || T0 != -1. || width != -1.){
					printf("\n *** ParameterMap contains insufficient information to modify the EoS ***\n");
					printf(" *** Required variables: EQOFST_BUMP_START EQOFST_BUMP_WIDTH EQOFST_BUMP_HEIGHT ***\n");
					printf("Continuing without modification....\n\n");
				}
					// else just ignore the ADD_BUMP flag
				return;
			}
			else if (height > 1. || height < 0.) {
				printf("\n *** EQOFST_BUMP_HEIGHT out of bounds (must be between 0 and 1) ***\n");
				printf("Continuing without modification....\n\n");
			}
			
				// modify height
			else addBump(T0,width,height);
		}
		
		lastAccess = 1;
		
			// generate splines for each equation of state variable
		spline(ed,temp,t2,
			   getDeriv(ed[0],ed[1],ed[2],temp[0],temp[1],temp[2]),
			   getDeriv(ed[aSize-1],ed[aSize-2],ed[aSize-3],temp[aSize-1],temp[aSize-2],temp[aSize-3]));
		spline(ed,pr,p2,
			   getDeriv(ed[0],ed[1],ed[2],pr[0],pr[1],pr[2]),
			   getDeriv(ed[aSize-1],ed[aSize-2],ed[aSize-3],pr[aSize-1],pr[aSize-2],pr[aSize-3]));
		spline(ed,sd,s2,
			   getDeriv(ed[0],ed[1],ed[2],sd[0],sd[1],sd[2]),
			   getDeriv(ed[aSize-1],ed[aSize-2],ed[aSize-3],sd[aSize-1],sd[aSize-2],sd[aSize-3]));
		
			//		eAtT = 0.;
			//		printEos("preTest.dat");
		
	}
	
		// so we don't have to find these each time
	eAtT = getEGivenT(svSwitchTemp);
	sAtT = getSGivenT(svSwitchTemp);
	
	printEos("mEos.txt");
}

CEosCalculator::~CEosCalculator() {
		// only if we generated the arrays
	if (mLatEos){
		delete [] temp;
		delete [] ed;
		delete [] pr;
		delete [] sd;
		delete [] t2;
		delete [] e2;
		delete [] p2;
		delete [] s2;
	}
}

CEosCalculator::CEosCalculator() {
	lastAccess = 1;
}

CEosCalculator::CEosCalculator(CEosCalculator& p) {
	lastAccess = p.lastAccess;
}

CEosCalculator::CEosCalculator(CEosCalculator* p) {
	lastAccess = p->lastAccess;
}

double CEosCalculator::getCs2(double e){
#ifdef FLAT_EOS
	return (1./3.);
#endif
	
	if (!mLatEos) 
		return (1./3.);
	else {
		double value = splintDeriv(pr,p2,e);
		if (value == 0.)
				//			return pr[0]/ed[0];
			return getCs2(ed[0]+1E-11);
		return value;
	}
}

double CEosCalculator::getP(double e) {
#ifdef FLAT_EOS
	return (e/3.);
#endif
	
	if (!mLatEos) 
		return (e/3.);
	else {
		double value = splint(pr,p2,e);
		if (value == 0.)
			return (pr[0]/ed[0])*e;
		return value;
	}
}

double CEosCalculator::getS(double e)  {
#ifdef FLAT_EOS
	return (e + getP(e))/getT(e);
#endif
	
	if (!mLatEos) 
		return (e + getP(e))/getT(e);
	else {
		double value = splint(sd,s2,e);
		if (value == 0.)
			return sd[0] * pow(e/ed[0],0.75);
		return value;
	}
}

double CEosCalculator::getTIS(double e) {
#ifdef FLAT_EOS
	return 2.;
#endif
	
	if (mSVRatio == 0.)
		return 1.;	
	if (getT(e) < svSwitchTemp || bRealisticEtaS)
		return 3.*getSV(e)/(e+getP(e));
	else 
		return 3.*mSVRatio*JHBARC/getT(e);
}

double CEosCalculator::getTISB(double e) {
#ifdef FLAT_EOS
	return 2.;
#endif
	
	if (mBVRatio == 0.)
		return 1.;
	
		//return getTIS(e);
	return 1.0;

	if (getT(e) < svSwitchTemp || bRealisticEtaS)
		return 3.*getBV(e)/(e+getP(e));
	else 
		return (3.*mBVRatio*JHBARC)/getT(e);
}

double CEosCalculator::getSV(double e) {
#ifdef FLAT_EOS
	return mSVRatio*e;
#endif
	
	double localT = getT(e);
	
	if (localT < svSwitchTemp)
		return JHBARC*mSVRatio*(sAtT/eAtT)*e;
	else if (bRealisticEtaS) {
		if (localT < mSvHighT) {
			if (svSwitchTemp == 0.)
				return JHBARC*getS(e)*(mSVRatio + (0.6-mSVRatio)/(1.+exp(-(0.1-localT)/0.02)));
			else 
				return JHBARC*getS(e)*mSVRatio;
		}
		else if(svSwitchTemp == 0.)
			return JHBARC*getS(e)*(mSVRatio + (0.6-mSVRatio)/(1.+exp(-(0.1-localT)/0.02)) + mSvHighTSlope*log(localT/mSvHighT));
		else 
			return JHBARC*getS(e)*(mSVRatio + mSvHighTSlope*log(localT/mSvHighT));
	}
	else 
		return JHBARC*getS(e)*mSVRatio;
}

double CEosCalculator::getBV(double e) {
#ifdef FLAT_EOS
	return mBVRatio*e;
#endif
		
		//	return 2.*getSV(e)*(1./3. - getCs2(e));
	
	return JHBARC * mBVRatio * getS(e) * exp( -0.5*pow((getT(e)-mBVCent)/mBVWidth,2));
	/*
	if (getT(e) < svSwitchTemp)
		return JHBARC*mBVRatio*(sAtT/eAtT)*e;
	else 
		return JHBARC*mBVRatio*getS(e);
	 */
}

double CEosCalculator::getT(double e) {
#ifdef FLAT_EOS
	return 0.5179033247*pow(e,0.25)*pow(JHBARC,0.75);
#endif
	
	// e = (16 + 21/2 * NDOF)* PI^2/30 * T^4/HBARC^3 for NDOF=2.5
	if (!mLatEos) return 0.5179033247*pow(e,0.25)*pow(JHBARC,0.75);  
	else {
		double value = splint(temp,t2,e);
		if (value == 0.)
			return temp[0] * pow(e/ed[0],0.25);
		return value;
	}
}

double CEosCalculator::getSigmaA(double e) {
	return getISAlpha(e)/sqrt(getS(e));
}

double CEosCalculator::getSigmaB(double e) {
	return getISGamma(e)/sqrt(getS(e));
}

double CEosCalculator::getISAlpha(double e){
	if (mSVRatio == 0.) 
		return 1.;
	if (ISSCALE == 'j') 
		return e+getP(e);
	else if (ISSCALE == 'e') 
		return e;
	else if (ISSCALE == 's') 
		return getS(e);
	else 
		return 1.;
}

double CEosCalculator::getISBeta(double e) {
//	return getSV(e)/(getISAlpha(e)*getTIS(e));
	return getSV(e)/getTIS(e);
}

double CEosCalculator::getISGamma(double e){
	return getBV(e)/getTISB(e);
}

double CEosCalculator::getDISAlphaDE(double e) {
	if      (ISSCALE == 'j') return 1.+getCs2(e);
	else if (ISSCALE == 'e') return 1.;
	else if (ISSCALE == 's') return (1./(4.*e));
	else if (ISSCALE == 'c') return 0.;
	else{
	    printf("CEos::getDISAlphaDE, ISSCALE not equl to j, e, s or c\n");
	    return ( getISAlpha(1.001*e) - getISAlpha(0.999*e))/(0.002*e);
	}
}

double CEosCalculator::getDISGammaDE(double e) {
	if      (ISSCALE == 'j') return 1.+getCs2(e);
	else if (ISSCALE == 'e') return 1.;
	else if (ISSCALE == 's') return (1./(4.*e));
	else if (ISSCALE == 'c') return 0.;
	else{
	    printf("CEos::getDISAlphaDE, ISSCALE not equl to j, e, s or c\n");
	    return ( getISAlpha(1.001*e) - getISAlpha(0.999*e))/(0.002*e);
	}
}

double CEosCalculator::getAMax(double e){
	return 2.*getP(e)/ROOT3;
}

double CEosCalculator::getBMax(double e){
	return getP(e);
}

double CEosCalculator::getDAMaxDE(double e){
	return 2.*getCs2(e)/ROOT3;
}

double CEosCalculator::getDBMaxDE(double e){
	return getCs2(e);
}

double CEosCalculator::getTA(double e) {
  if (mLatEos) 
	  return e-3.*getP(e);
  else 
	  return 0.;
}

double CEosCalculator::getEGivenT(double T) {
  if (!mLatEos) return 1811.507976 * pow(T,4);
  else {
    int i=1;
	for (;i<aSize;i++) if (temp[i] >= T) break;
	if (T==temp[i]) return ed[i];
	return ((temp[i] - T)*ed[i-1] + (T-temp[i-1])*ed[i])/(temp[i] - temp[i-1]);
  }
}

double CEosCalculator::getSGivenT(double T) {
	if (!mLatEos) 
		return (4./3.)*getEGivenT(T)/T;
	else {
		int i=1;
		for (;i<aSize;i++) if (temp[i] >= T) break;
		if (T==temp[i]) return sd[i];
		return ((temp[i] - T)*sd[i-1] + (T-temp[i-1])*sd[i])/(temp[i] - temp[i-1]);
	}
}

double CEosCalculator::getEGivenS(double mS) {
	if (!mLatEos) 
		return 1.;
	else {
		if (mS < sd[0])
			return ed[0]*pow(mS/sd[0],4./3.);
		
		int i=1;
		for (;i<aSize;i++) if (sd[i] >= mS) break;
		if (mS==sd[i]) return ed[i];
		return ((sd[i] - mS)*ed[i-1] + (mS-sd[i-1])*ed[i])/(sd[i] - sd[i-1]);
	}
}

void CEosCalculator::addBump(double mT0, double mW, double fH) {
	int mPow = 2;
	float T0 = float(mT0);
	float width = float(mW);
	
	float minH, maxH;	
	getMaxAmp(T0,width,minH,maxH);
	float height = (1.-fH)*minH + fH*maxH;
	
	float a, amp;
	if (bGaussTail) {
		a = float(mPow)/(width*width);
		amp = height * pow( a/float(mPow), mPow) * exp(mPow);
	}
	else if (bExpTail) {
		a = ((float)mPow+1.)/width;
		amp = height * pow( a/(float(mPow) + 1.), mPow+1) * exp( float(mPow) + 1.);
	}
	else if (bPowerTail) {
		a = width;
		amp = height * pow(2./width, mPow);
	}
	
	for (int j=0;j<aSize;j++) {
		float x = temp[j] - T0;
		if (temp[j] < T0) continue;
		else {
			
				// modification is to s/T^3
			if (bGaussTail)
				sd[j] = sd[j]/pow(temp[j],3) + amp*pow(x,2*mPow)*exp(-a*x*x);
			else if (bExpTail)
				sd[j] = sd[j]/pow(temp[j],3) + amp*pow(x,mPow+1)*exp(-a*x);
			else if (bPowerTail)
				sd[j] = sd[j]/pow(temp[j],3) + amp*pow(x/(1. + (x/width)*(x/width)), mPow);
			
				// decreasing s/T^3 == cs2 < 0 || 1/3 < cs2
			if (sd[j] < sd[j-1]) {
				cout << "\n ***** Speed of Sound out of bounds! ***** Aborting! ******* \n\n";
				exit(1);
			}
			
				// s/T^3 => s
			sd[j] *= pow(temp[j],3);
			
			pr[j] = pr[j-1] + 0.5 * (sd[j] + sd[j-1]) * (temp[j] - temp[j-1]);
		}
		
		ed[j] = temp[j]*sd[j] - pr[j];
	}
	
}

void CEosCalculator::getMaxAmp(double T0, double width, double &minH, double &maxH ) {
	float mMinH;
	float mMaxH;
	getMaxAmp( float(T0), float(width), mMinH, mMaxH);
	minH = double(mMinH);
	maxH = double(mMaxH);
}

void CEosCalculator::getMaxAmp(float T0, float mWidth, float &minH, float &maxH) {
	float dH = 10.;
	
		// find max height
	for (float height=dH; ; height+=dH) 
		if (!bounds(T0,mWidth,height)) break;
		else maxH = height;
		// find min height 
	for (float height=-dH; ; height-=dH) 
		if (!bounds(T0,mWidth,height)) break;
		else minH = height;
}

bool CEosCalculator::bounds (float T0, float width, float height) {
	int mPow=2;
	
	float a, amp;
	
	if (bExpTail) {
		a = ((float)mPow+1.)/width;
		amp = height * pow( a/(float(mPow) + 1.), mPow+1) * exp( float(mPow) + 1.);
	}
	else if (bGaussTail) {
		a = float(mPow)/(width*width);
		amp = height * pow( a/float(mPow), mPow) * exp(mPow);
	}
	else if (bPowerTail) {
		a = width;
		amp = height * pow(2./width, mPow);
	}
	else {
		printf("\n ****** Unknown Option! *******\n\n");
		exit(1);
	}
	
	float mST[aSize];
	
	for (int j=0;j<aSize;j++) {
		float x = temp[j] - T0;
		if (x < 0.)
			mST[j] = sd[j]/pow(temp[j],3);
		else {
			if (bGaussTail)
				mST[j] = sd[j]/pow(temp[j],3) + amp*pow(x,2*mPow)*exp(-a*x*x);
			else if (bExpTail)
				mST[j] = sd[j]/pow(temp[j],3) + amp*pow(x,mPow+1)*exp(-a*x);
			else if (bPowerTail)
				mST[j] = sd[j]/pow(temp[j],3) + amp*pow(x/(1. + (x/width)*(x/width)), mPow);
			if ( mST[j] < mST[j-1])
				return false;
		}
	}
	return true;
}

	// fills lastAccess with the array index just lower than x
	// returns -1 for below ed[0], returns aSize-1 for above ed[aSize-1]
void CEosCalculator::search(double e) {
		//	printf("cec::search(%g)...\n",e);
	
	if (lastAccess==-1) lastAccess++;
	
	if (e > ed[lastAccess]) {
		for (;lastAccess<aSize-1;lastAccess++)
			if (ed[lastAccess+1] > e) 
				break;
	}
	else {
		for (;lastAccess>=0;lastAccess--)
			if (ed[lastAccess] < e) 
				break;
	}
	
	if (lastAccess == 0 && ed[lastAccess] > e)
		lastAccess--;
}

	// compute the second derivatives via tridiagonalization
	// called once during construction
	// required for splint() interpolation 
	// Algorithm from Numerical Recipes for C++ Press et al. (Section 3.3)
	// Inputs are (1,2) coordinates of the data points to be fit
	// (3) second derivative at each data point
	// (4,5) first derivative at the first and last data point
void CEosCalculator::spline(double *x, double *y, double *y2, double yp1, double ypA){
	double p,sig;
	double u[aSize-1];
	
	y2[0]=-0.5;
	u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	
	for (int i=1;i<aSize-1;i++) {
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sig*y2[i-1]+2.0;
		y2[i] = (sig-1.0)/p;
		u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
	}
	
	double qn = 0.5; 
	double un = (3.0/(x[aSize-1]-x[aSize-2]))*(ypA - (y[aSize-1]-y[aSize-2])/(x[aSize-1]-x[aSize-2]));
	
	y2[aSize-1] = (un-qn*u[aSize-2])/(qn*y2[aSize-2]+1.0);
	
		// backsubstitution
	for (int k=aSize-2;k>=0;k--)
		y2[k] = y2[k]*y2[k+1] + u[k];
}

	// compute the interpolated value using spline (from NR)
	// Algorithm from Numerical Recipes for C++, Press et al. (Section 3.3)
double CEosCalculator::splint(double *y, double *y2, double e){
	search(e);
	
	if (lastAccess==-1 || lastAccess==aSize-1) {
			//		printf("splint() called out of bounds....\n");
		return 0.;
	}
	
	double h,a,b;
	h = ed[lastAccess+1] - ed[lastAccess];
	a = (ed[lastAccess+1]-e)/h;
	b = (e - ed[lastAccess])/h;
	
	return a*y[lastAccess]+b*y[lastAccess+1] 
	+ (a*(a*a-1.)*y2[lastAccess] + b*(b*b-1.)*y2[lastAccess+1])*(h*h)/6.0;
}

	// instead of value, returns derivative (dydx) at x=e
	// used for interpolating c_s^2
	// Algorithm from Numerical Recipes for C++, Press et al. (Section 3.3)
double CEosCalculator::splintDeriv(double *y, double *y2, double e){
	search(e);
	
	if (lastAccess==-1 || lastAccess==aSize-1) {
			//		printf("splintDeriv() called out of bounds....\n");
		return 0.;
	}
	
	double h,a,b;
	h = ed[lastAccess+1] - ed[lastAccess];
	a = (ed[lastAccess+1]-e)/h;
	b = (e - ed[lastAccess])/h;
	
	return (y[lastAccess+1]-y[lastAccess])/h 
	- (3.*a*a-1.)*h*y2[lastAccess]/6.0 + (3.*b*b-1.)*h*y2[lastAccess+1]/6.0;
}

void CEosCalculator::printEos(const char * mFileName) {
	FILE *mF = fopen(mFileName,"w");
	printEos(mF);
	fclose(mF);
}

void CEosCalculator::printEosEnergyUnits(const char * mFileName) {
	FILE *mF = fopen(mFileName,"w");
	printEosEnergyUnits(mF);
	fclose(mF);
}

void CEosCalculator::printEos(FILE* mF) {
	if (mLatEos)
		for (int i=0;i<aSize;i++)
			fprintf(mF,"%g %g %g %g %g %g %g\n",temp[i],sd[i],ed[i],pr[i],getCs2(ed[i]),getSV(ed[i]),getBV(ed[i]));
	else {
		double deltaE = 1E-4;
		for (double de=deltaE;de < 1.; de+=deltaE)
			fprintf(mF,"%g %g %g %g %g %g %g\n",getT(de),getS(de),de,getP(de),getCs2(de),getSV(de),getBV(de));
				//fprintf(mF,"%g %g %g\n",de,getP(de),getCs2(de));
	} 
}

void CEosCalculator::printEosEnergyUnits(FILE* mF) {
	double mHBCubed = pow(JHBARC,3);
	if (mLatEos)
		for (int i=0;i<aSize;i++) {
			double mT4 = pow(temp[i],4);
			fprintf(mF,"%g %g %g %g %g %g\n",temp[i],sd[i]*(temp[i]*mHBCubed/mT4),ed[i]*(mHBCubed/mT4),
					pr[i]*(mHBCubed/mT4),getCs2(ed[i]),getSV(ed[i])*(temp[i]*mHBCubed/mT4));
		}
	else 
		for (double de=0.001;de < 1.; de+=0.001) {
			double mT4 = pow(getT(de),4);
			fprintf(mF,"%g %g %g %g %g\n",getT(de),getS(de)*(getT(de)*mHBCubed/mT4),de*(mHBCubed/mT4),
					getP(de)*(mHBCubed/mT4),getCs2(de));
		} 
}

	// calculates derivative dydx(x=x1) using quadratic ansatz
double CEosCalculator::getDeriv(double x1, double x2, double x3, double y1, double y2, double y3){
	double a = ((y3-y1)*(x2-x1) - (y2-y1)*(x3-x1))/((x3*x3-x1*x1)*(x2-x1) - (x2*x2-x1*x1)*(x3-x1));
	double b = ((y2-y1)*(x3*x3-x1*x1) - (y3-y1)*(x2*x2-x1*x1))/((x2-x1)*(x3*x3-x1*x1) - (x3-x1)*(x2*x2-x1*x1));
	return 2*a*x1+b;
}

	// verifies a file's existence
	// swiped wholesale from: http://www.techbytes.ca/techbyte103.html
bool CEosCalculator::checkFile(string fName) {
		//return true;
	
	struct stat stFileInfo; 
		// Attempt to get the file attributes 
		// if we succeed, the file exists
	if(stat(fName.c_str(),&stFileInfo) == 0) 
		return true; 
	else 
		return false; 
}

CEos::CEos() {
	mEosC = new CEosCalculator();
}

CEos::CEos(parameterMap* p){
	mEosC = new CEosCalculator(p);
	
		// these variables don't necessarily get set
	aMax = 0.;
	bMax = 0.;
	dAMaxDE = 0.;
	dBMaxDE = 0.;
	
	mISMax = parameter::getB(*p,"HYDRO_ISMAX",false);
}

	// calculator destructor automagically called
CEos::~CEos(){
	delete mEosC;
}

CEos::CEos(CEos& e) {
	cS2 = e.cS2;
	P = e.P;
	S = e.S;
	tIS = e.tIS;
	tISB = e.tISB;
	SV = e.SV;
	BV = e.BV;
	T = e.T;
	
	sigmaA = e.sigmaA; 
	sigmaB = e.sigmaB;
	alphaIS = e.alphaIS;
	gammaIS = e.gammaIS;
	betaIS = e.betaIS;
	aMax = e.aMax;
	bMax = e.bMax;

	dAlphaISDE = e.dAlphaISDE;
	dGammaISDE = e.dGammaISDE;
	dAMaxDE = e.dAMaxDE;
	dBMaxDE = e.dBMaxDE;
	
	mEosC = e.mEosC;
}

CEos::CEos(CEos* p) {
	cS2 = p->cS2;
	P = p->P;
	S = p->S;
	tIS = p->tIS;
	tISB = p->tISB;
	SV = p->SV;
	BV = p->BV;
	T = p->T;
	
	sigmaA = p->sigmaA; 
	sigmaB = p->sigmaB;
	alphaIS = p->alphaIS;
	gammaIS = p->gammaIS;
	betaIS = p->betaIS;
	aMax = p->aMax;
	bMax = p->bMax;
	
	dAlphaISDE = p->dAlphaISDE;
	dGammaISDE = p->dGammaISDE;
	dAMaxDE = p->dAMaxDE;
	dBMaxDE = p->dBMaxDE;
	
		// make a new CEosCalculator but remember where we were
	mEosC = new CEosCalculator(p->mEosC);
}

	// Use CEosCalculator to fill CEos variables
void CEos::fill (double e) {
	P = getP(e);
	cS2	= getCs2(e);
	S = getS(e);
	T = getT(e);
	SV	= getSV(e);
	BV	= getBV(e);
	
	tIS	= getTIS(e);
	tISB  = getTISB(e);
	
	alphaIS= getISAlpha(e);
	gammaIS= getISGamma(e);
	
	dAlphaISDE = getDISAlphaDE(e);
	dGammaISDE = getDISGammaDE(e);
	
	sigmaA = alphaIS/sqrt(S);
	sigmaB = gammaIS/sqrt(S);
	
//	betaIS = SV/(alphaIS*tIS);
	betaIS = SV/tIS;
	
		// the variables for max
	if (mISMax) {
		aMax = getAMax(e);
		bMax = getBMax(e);
		dAMaxDE = getDAMaxDE(e);
		dBMaxDE = getDBMaxDE(e);
	}	
}

	// for cleanliness cells ask CEos to adjust
	// if the shear viscosity is position dependent
	// NOTE :: changes here need to migrate to:
	// CCell::getSVCalc(double) and CCell::getTISCalc(double)
void CEos::trimSV(double r) {
	SV  /= (1. + exp( (r - 10.)/0.6));
	tIS /= (1. + exp( (r - 10.)/0.6));
}

	// static variables for the calculator
bool CEosCalculator::mLatEos;
int CEosCalculator::aSize;
double CEosCalculator::mSVRatio, CEosCalculator::mBVRatio, CEosCalculator::eAtT;
double CEosCalculator::sAtT, CEosCalculator::svSwitchTemp, CEosCalculator::mFoTemp;
double CEosCalculator::mSvHighT, CEosCalculator::mSvHighTSlope;
double CEosCalculator::mBVCent, CEosCalculator::mBVWidth;
bool CEosCalculator::bExpTail, CEosCalculator::bPowerTail;
bool CEosCalculator::bGaussTail, CEosCalculator::bRealisticEtaS;

	// equation of state variables
double *CEosCalculator::temp, *CEosCalculator::ed, *CEosCalculator::pr, *CEosCalculator::sd;
	// second derivatives for splines
double *CEosCalculator::t2, *CEosCalculator::e2, *CEosCalculator::p2, *CEosCalculator::s2;

	// static variables for eos
bool CEos::mISMax;