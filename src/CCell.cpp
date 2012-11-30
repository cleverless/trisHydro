/*  CCell.cpp
 *  isHydro
 *  Created by Joshua Vredevoogd on 2/11/09.
 */

#include "CCell.h"

CCellHelper::CCellHelper(){
	/* noop */
}

CCellHelper::CCellHelper(CCellHelper& h){
	/* noop */
}

CCellHelper::CCellHelper(CCellHelper* p){
	/* noop */
}

CCellHelper::~CCellHelper(){
	/* noop */
}

CCell::CCell(parameterMap* pM){
	helper = new CCellHelper();
	pMap = pM;
	paramFill();
	
	helper->g[0] = 1.;
	helper->g[1] = -1.;
	helper->g[2] = -1.;
	helper->g[3] = -1.;
		
}

CCell::CCell(double t0,double n0, double x0, double y0) { 
	if (mLinT) x[0] = t0;
	else if (mLogT) x[0] = log(t0/mT0);
	else if (mLogSinhT) x[0] = log( sinh(t0));
	
	x[3] = n0;
	x[1] = x0;
	x[2] = y0;
	
		//	printf("CC:CC() : %g %g %g %g\n",x[0],x[1],x[2],x[3]);
	
	for (int i=0;i<11;i++) 
		s[i] = 0.;
	
	active = true;
}

CCell::CCell(CCell* cell) {
	for (int i=0;i<4;i++) 
		x[i] = cell->x[i];
	for (int i=0;i<11;i++) 
		s[i] = cell->s[i];
	for (int i=0;i<3;i++)
		dUdT[i] = cell->dUdT[i];
	
	active = cell->active;
}  

CCell::~CCell() {
		//noop
}

double CCell::getTau() {
	if (mLinT) return x[0];
	else if (mLogT) return mT0*exp(x[0]);
	else if (mLogSinhT) return asinh(exp(x[0]));
	else{
		printf("CCell::getTau, didn't satisfy conditions\n");
		return 0.0;
	}
}

void CCell::paramFill() {
	
	mDebug  = parameter::getB(*pMap,"HYDRO_DEBUG",false);
	mSVTrim = parameter::getB(*pMap,"HYDRO_SVTRIM",false);
	mSVTrimInit = parameter::getB(*pMap,"HYDRO_SVTRIMINIT",false);
	mISMax = parameter::getB(*pMap,"HYDRO_IS_MAX",false);
	
	mViscNS = parameter::getB(*pMap,"HYDRO_VISCNS",false);
	mPureBjorken = parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false);
	mBjorken = parameter::getB(*pMap,"HYDRO_BJORKEN",true);
	mLinT = parameter::getB(*pMap,"HYDRO_LINT",false);
	mLogT = parameter::getB(*pMap,"HYDRO_LOGT",false);
	mLogSinhT = parameter::getB(*pMap,"HYDRO_LOGSINHT",true);
	mISVort = parameter::getB(*pMap,"HYDRO_IS_VORT",false);
	mISMax = parameter::getB(*pMap,"HYDRO_IS_MAX",false);
	mLaxFried = parameter::getB(*pMap,"HYDRO_LAX_FRIED",false);
	mSlopeLimit = parameter::getB(*pMap,"HYDRO_SLOPE_LIMIT",false);
	
	mT0 = parameter::getD(*pMap,"HYDRO_T0",1.0);
	mSVRatio = parameter::getD(*pMap,"HYDRO_SVRATIO",0.0);
	mBVRatio = parameter::getD(*pMap,"HYDRO_BVRATIO",0.0);
	
	mInitFlow = parameter::getD(*pMap,"HYDRO_INIT_FLOW",0.0);
	mInitNS = parameter::getD(*pMap,"HYDRO_INIT_NS",0.0);
	mInitNSTxxP = parameter::getD(*pMap,"HYDRO_INIT_NS_TXXP",1.0);
	mSlopeLimitTheta = parameter::getD(*pMap,"HYDRO_SLOPE_LIMIT_THETA",1.1);
	
	mOctant = parameter::getI(*pMap,"HYDRO_OCTANT",3);
	
	dx[1] = parameter::getD(*pMap,"HYDRO_DX",0.1);
	dx[2] = parameter::getD(*pMap,"HYDRO_DY",0.1);
	dx[3] = parameter::getD(*pMap,"HYDRO_DN",0.1);
}

void CCell::update() {
	
	fillEosVar();
	
	helper->g[0] =  1.;
	helper->g[1] = -1.;
	helper->g[2] = -1.;
	
	if (mBjorken)
		helper->g[3] = -getTau()*getTau();
	else 
		helper->g[3] = -1.;
	
	if (mViscNS) 
		initNS();
	
	helper->gamma = sqrt(1.+s[1]*s[1]+s[2]*s[2]-helper->g[3]*s[3]*s[3]);
	helper->y = sqrt(s[5]*s[5] + s[6]*s[6]);
	helper->a = getAMax() * tanh( getISAlpha() * helper->y / getAMax());
	
	fillPixy();
	
	calcDeriv();
	fillM();
}

	// version of update not relient on any derivative information
void CCell::selfUpdate() {
	fillEosVar();

	if (mBjorken)
		helper->g[3] = -getTau()*getTau();

	if (mViscNS) 
		initNS();
	
	helper->gamma = sqrt(1.+s[1]*s[1]+s[2]*s[2]-helper->g[3]*s[3]*s[3]);
	helper->y = sqrt(s[5]*s[5] + s[6]*s[6]);
	helper->a = getAMax() * tanh( getISAlpha() * helper->y / getAMax());
	
	fillPixy();
}

	// for calculating slope limited derivatives
double CCell::minmod( double d1, double d2, double d3) {
	d2 *= mSlopeLimitTheta;
	d3 *= mSlopeLimitTheta;
	
	if (d1 > 0.) {
		if (d2 < 0. || d3 < 0.) return 0.;
		else return min(d1, d2, d3);
	}
	else {
		if (d2 > 0. || d3 > 0.) return 0.;
		else return max(d1, d2, d3);
	}
}

double CCell::max(double a, double b, double c){
	 if ( b > a) {
		 if (b > c) return b;
		 else return c;
	 }
	else {
		if (a > c) return a;
		else return c;
	}
}

double CCell::min(double d1, double d2, double d3){
	if ( d2 < d1) {
		if (d2 < d3) return d2;
		else return d3;
	}
	else {
		if (d1 < d3) return d1;
		else return d3;
	}
}

double CCell::minmod( double d1, double d2) {
	d2 *= mSlopeLimitTheta;

	if (d1 > 0.) {
		if (d2 < 0.) return 0.;
		else return min(d1, d2);
	}
	else {
		if (d2 > 0.) return 0.;
		else return max(d1, d2);
	}
}

double CCell::max(double d1, double d2){
	if (d2 > d1) return d2;
	else return d1;
}

double CCell::min(double d1, double d2){
	if (d2 < d1) return d2;
	else return d1;
}

void CCell::calcDeriv() {
	if (!active) {
		cout << "\ncalcDeriv() called in dead Cell! (" << getTau() 
		<< " " << x[1] << " " << x[2] << " " << x[3] <<")\n";
		selfPrint();
		fflush(stdout);
	}
	
	fillDS();
	
	for (int i=0;i<3;i++) 
		helper->dTISS[i] = getDISAlphaDE() * getE() * helper->dS[0][i];

	fillDU();
	fillDPixy();
}

void CCell::fillDS(){
	for (int j=0;j<11;j++)
		for (int i=0;i<3;i++)
			if (neighbors[i][1]->getActive() && neighbors[i][0]->getActive()) {
				if (ISSCALE == 'c' && j>3 ) {
					helper->dS[j][i] = getE()*(  neighbors[i][1]->s[j]/neighbors[i][1]->getE() 
											   - neighbors[i][0]->s[j]/neighbors[i][0]->getE() ) / (2.*dx[i+1]);
					helper->dS[j][i] += s[j]*helper->dS[0][i];
				}
				else 
					helper->dS[j][i] = (neighbors[i][1]->s[j] - neighbors[i][0]->s[j])/ (2.*dx[i+1]);
			}
			else if (!neighbors[i][1]->getActive()) {
				
				if (!neighbors[i][0]->getActive()){
					printf("\ntrouble in the old neighborhood!\n");
					selfPrint();
				}
				
				if (ISSCALE == 'c' && j>3 ) {
						//helper->dS[j][i] = (3.*s[j] + getE() * ( neighbors[i][0]->neighbors[i][0]->s[j]/neighbors[i][0]->neighbors[i][0]->getE() 
						//									- 4.*neighbors[i][0]->s[j]/neighbors[i][0]->getE())) / (2.*dx[i+1]);
					
						//					helper->dS[j][i] = (s[j] - getE() * neighbors[i][0]->s[j]/neighbors[i][0]->getE())/dx[i+1];
					
					helper->dS[j][i] = minmod( (3.*s[j] + getE() * ( neighbors[i][0]->neighbors[i][0]->s[j]/neighbors[i][0]->neighbors[i][0]->getE() 
																	- 4.*neighbors[i][0]->s[j]/neighbors[i][0]->getE())) / (2.*dx[i+1]),
											  (s[j] - getE() * neighbors[i][0]->s[j]/neighbors[i][0]->getE())/dx[i+1]);
					helper->dS[j][i] += s[j]*helper->dS[0][i];
				}
				else
						//helper->dS[j][i] = (3.*s[j] - 4.*neighbors[i][0]->s[j] 
						//				+ neighbors[i][0]->neighbors[i][0]->s[j]) / (2.*dx[i+1]);
					
						//helper->dS[j][i] = (s[j] - neighbors[i][0]->s[j])/dx[i+1];
					
					helper->dS[j][i] = minmod((3.*s[j] - 4.*neighbors[i][0]->s[j] 
											   + neighbors[i][0]->neighbors[i][0]->s[j]) / (2.*dx[i+1]),
											  (s[j] - neighbors[i][0]->s[j])/dx[i+1]);
			}
			else {
				if (ISSCALE == 'c' && j>3) {
						//helper->dS[j][i] = (-3.*s[j] + getE() * (-neighbors[i][1]->neighbors[i][1]->s[j]/neighbors[i][1]->neighbors[i][1]->getE() 
						//									+ 4.*neighbors[i][1]->s[j]/neighbors[i][1]->getE())) / (2.*dx[i+1]);

						//helper->dS[j][i] = (getE() * neighbors[i][1]->s[j]/neighbors[i][1]->getE() - s[j])/dx[i+1];
					
					helper->dS[j][i] = minmod ((-3.*s[j] + getE() * (-neighbors[i][1]->neighbors[i][1]->s[j]/neighbors[i][1]->neighbors[i][1]->getE() 
																	 + 4.*neighbors[i][1]->s[j]/neighbors[i][1]->getE())) / (2.*dx[i+1]),
											  (getE() * neighbors[i][1]->s[j]/neighbors[i][1]->getE() - s[j])/dx[i+1]);
					
					helper->dS[j][i] += s[j]*helper->dS[0][i];
				}
				else 
						//helper->dS[j][i] = (-3.*s[j] + 4.*neighbors[i][1]->s[j] 
						//				- neighbors[i][1]->neighbors[i][1]->s[j]) / (2.*dx[i+1]);

						//helper->dS[j][i] = (neighbors[i][1]->s[j] - s[j])/dx[i+1];
					
					helper->dS[j][i] = minmod ( (-3.*s[j] + 4.*neighbors[i][1]->s[j] 
												 - neighbors[i][1]->neighbors[i][1]->s[j]) / (2.*dx[i+1]),
											   (neighbors[i][1]->s[j] - s[j])/dx[i+1]);
			}
	
	if (mPureBjorken)
		for (int j=0;j<11;j++)
			helper->dS[j][2] = 0.;
}

void CCell::fillPixy(){
	helper->Pixy[1][1] = s[4];
	helper->Pixy[2][2] = s[5];
	helper->Pixy[3][3] = s[6];
	helper->Pixy[1][2] = s[7];
	helper->Pixy[2][1] = helper->Pixy[1][2];
	helper->Pixy[1][3] = s[8];
	helper->Pixy[3][1] = helper->Pixy[1][3];
	helper->Pixy[2][3] = s[9];
	helper->Pixy[3][2] = helper->Pixy[2][3];
	
	for (int i=1;i<4;i++){
		helper->Pixy[0][i] = (s[1]*helper->Pixy[1][i] + s[2]*helper->Pixy[2][i] - helper->g[3]*s[3]*helper->Pixy[3][i])/helper->gamma;
		helper->Pixy[i][0] = helper->Pixy[0][i];
	}
	
		//	helper->Pixy[0][0] = helper->Pixy[1][1] + helper->Pixy[2][2] - helper->g[3]*helper->Pixy[3][3];
	helper->Pixy[0][0] = (s[1]*helper->Pixy[0][1]+s[2]*helper->Pixy[0][2]-helper->g[3]*s[3]*helper->Pixy[0][3])/helper->gamma;
	
	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++)
			helper->Pixy[i][j] *= getISAlpha();
	
	helper->BPi = s[10]*getISAlpha();
}

void CCell::fillDU(){
		// copy and reorient
	for (int i=0;i<3;i++)
		for (int j=1;j<4;j++)
			helper->dU[i+1][j] = helper->dS[j][i];
	
		// initialize connections of time derivatives
	helper->dU[0][1] = 0.;
	helper->dU[0][2] = 0.;
	helper->dU[0][3] = 0.;
	
		// add affine connections
	if (mBjorken) {
		helper->dU[0][3]  = s[3]/getTau();
			//		helper->dU[3][0]  = getTau()*s[3];  captured in next loop
		helper->dU[3][3] += helper->gamma/getTau();
		helper->dU[0][0]  = getTau()*s[3]*s[3]/helper->gamma;
	}
	
	for (int i=1;i<4;i++){
		helper->dU[i][0] = 0.;
		for (int j=1;j<4;j++)
			helper->dU[i][0] -= helper->g[j]*s[j]*helper->dU[i][j]/helper->gamma;
	}
}

void CCell::fillDPixy(){
	
		// not all of these will be initialized elsewhere
	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++)
			helper->DPixy[0][i][j] = 0.;
	
		// derivatives of spatial elements
	for (int i=0;i<3;i++){
		helper->DPixy[i+1][1][1] = helper->dS[4][i]*getISAlpha() + helper->dTISS[i]*s[4];
		helper->DPixy[i+1][2][2] = helper->dS[5][i]*getISAlpha() + helper->dTISS[i]*s[5];
		helper->DPixy[i+1][3][3] = helper->dS[6][i]*getISAlpha() + helper->dTISS[i]*s[6];
		helper->DPixy[i+1][1][2] = helper->dS[7][i]*getISAlpha() + helper->dTISS[i]*s[7];
		helper->DPixy[i+1][1][3] = helper->dS[8][i]*getISAlpha() + helper->dTISS[i]*s[8];
		helper->DPixy[i+1][2][3] = helper->dS[9][i]*getISAlpha() + helper->dTISS[i]*s[9];
	}
	
		// add affine connections to spatial elements
	if (mBjorken) {
		for (int i=1;i<4;i++) {
			helper->DPixy[0][i][3] += helper->Pixy[i][3]/getTau();
			helper->DPixy[3][i][3] += helper->Pixy[0][i]/getTau();
		}
		
		helper->DPixy[0][3][3] += helper->Pixy[3][3]/getTau();
		helper->DPixy[3][3][3] += helper->Pixy[0][3]/getTau();
	}
	
		//enforce symmetry on spatial elements
	for (int i=0;i<4;i++)
		for (int j=1;j<4;j++)
			for (int k=j+1;k<4;k++)
				helper->DPixy[i][k][j] = helper->DPixy[i][j][k];

		// calculate temporal elements
	for (int i=0;i<4;i++){
		for (int j=1;j<4;j++) {	
			helper->DPixy[i][0][j] = 0.;
			for (int k=0;k<4;k++) 
				helper->DPixy[i][0][j] -= helper->g[k]*helper->Pixy[k][j]*helper->dU[i][k]/helper->gamma;
			for (int k=1;k<4;k++) 
				helper->DPixy[i][0][j] -= helper->g[k]*s[k]*helper->DPixy[i][k][j]/helper->gamma;
		}
	}
	
		// enforce symmetry on temporal elements
	for (int i=0;i<4;i++)
		for (int j=1;j<4;j++) 
			helper->DPixy[i][j][0] = helper->DPixy[i][0][j];
			
		// calc D_i Pi_{00} elements
	for (int i=0;i<4;i++)
		helper->DPixy[i][0][0] = helper->DPixy[i][1][1]+helper->DPixy[i][2][2]-helper->g[3]*helper->DPixy[i][3][3];

}

void CCell::fillM() {

	for (int i=0;i<12;i++) 
		helper->M[0][i] = getECoeff(i);
	
	for (int i=1;i<4;i++) 
		for (int j=0;j<12;j++) 
			helper->M[i][j] = getPCoeff(i,j);

		for (int i=4;i<11;i++) 
			for (int j=0;j<12;j++)
				helper->M[i][j] = getPiCoeff(i,j);

}

double CCell::getECoeff(int j) {
	double value = 0.;
	if (j==0) 
		value = getE()*helper->gamma;
	else if (j<4) 
		value = -helper->g[j]*( (s[j]/helper->gamma)*(getE()+getP()+helper->BPi-helper->Pixy[0][0]) 
							   + helper->Pixy[0][j]);
	else if (j==11) {
		for (int m=0;m<4;m++) {
			value -= (getE()+getP()+helper->BPi)*helper->dU[m][m];
			for (int n=0;n<4;n++)
				value += helper->g[m]*helper->Pixy[m][n]*helper->dU[n][m];
		}
		for (int m=1;m<4;m++)
			value -= s[m]*getE()*helper->dS[0][m-1];
		
		//BULKTAG
	}
		
	return value;
}

double CCell::getPCoeff(int i, int j) {
	double value = 0.;
	
	if (j==0) {
		value += s[i]*helper->gamma*getCs2();
		value += getDISAlphaDE()*helper->Pixy[i][0]/getISAlpha();
		value *= getE();
		value += s[i]*helper->gamma*getDISAlphaDE()*helper->BPi/getISAlpha();
	}
	else if (j<4) {
		if (i==j)
			value += helper->gamma*(getE()+getP()+helper->BPi);
		value += s[i]*helper->g[j]*(helper->Pixy[0][j] - (helper->Pixy[0][0]*s[j]/helper->gamma));
		value += (helper->g[j]/helper->gamma)*((helper->Pixy[i][0]*s[j]/helper->gamma) - helper->Pixy[i][j]);
	}
	else if (j==4){ //Pi_xx
		if (i==1) 
			value += getISAlpha()*s[1]/helper->gamma;
	}
	else if (j==5){ //Pi_yy
		if (i==2)
			value += getISAlpha()*s[2]/helper->gamma;
	}
	else if (j==6){ //Pi_zz
		if (i==3)
			value -= getISAlpha()*helper->g[3]*s[3]/helper->gamma;
	}
	else if (j==7){ //Pi_xy
		if (i==1)
			value += getISAlpha()*s[2]/helper->gamma;
		else if (i==2)
			value += getISAlpha()*s[1]/helper->gamma;
	}
	else if (j==8){ //Pi_xz
		if (i==1)
			value -= getISAlpha()*helper->g[3]*s[3]/helper->gamma;
		else if (i==3)
			value += getISAlpha()*s[1]/helper->gamma;
	}
	else if (j==9){ //Pi_yz
		if (i==2)
			value -= getISAlpha()*helper->g[3]*s[3]/helper->gamma;
		else if (i==3)
			value += getISAlpha()*s[2]/helper->gamma;
	}
	else if (j==10){
		value = -s[i]*helper->gamma;
	}
	else if (j==11){
		value -= (getE()+getP()+helper->BPi)*helper->gamma*helper->dU[0][i];
		value += getCs2()*getE()*helper->dS[0][i-1]/helper->g[i];
		value += (helper->dS[10][i-1]*getISAlpha() + helper->BPi*getDISAlphaDE()*helper->dS[0][i-1])/helper->g[i];
		
		//BULKTAG
		
		for (int m=1;m<4;m++){
			value -= (getE()+getP()+helper->BPi)*s[m]*helper->dU[m][i];
			value -= s[i]*s[m]*getCs2()*getE()*helper->dS[0][m-1];
			value -= s[i]*s[m]*(helper->dS[10][i-1]*getISAlpha() + helper->BPi*getDISAlphaDE()*helper->dS[0][i-1]);
		}

		for (int m=0;m<4;m++){
			value -= helper->DPixy[m][i][m];
			for (int n=0;n<4;n++){
				value -= s[i]*helper->g[n]*helper->Pixy[m][n]*helper->dU[m][n];
			}
		}
	}
	return value;
}

	// shear components
double CCell::getPiCoeff(int i, int j) {
	double value = 0.;
	
	if (i==10) return getPiCoeff(j);
	
	if (mSVRatio==0.) {
		if (i==j) return 1.;
		else return 0.;
	}
	
	if (j>3 && j<11) {
		// no equations mix time derivatives of Pi's.
	 	if (i!=j)
			return 0.;
		else 
			return helper->gamma*getISAlpha();
	}
	
	if (j==0) 
//		return -helper->gamma * s[i] * getDISAlphaDE();
		return 0.;

	// local indices to clean up code
	int mu, nu; 
	
	if (i==4){
		mu = 1;
		nu = 1;
	} 
	else if (i==5){
		mu = 2;
		nu = 2;
	} 
	else if (i==6){
		mu = 3;
		nu = 3;
	} 
	else if (i==7){
		mu = 1;
		nu = 2;
	} 
	else if (i==8){
		mu = 1;
		nu = 3;
	}
	else if (i==9){
		mu = 2;
		nu = 3;
	}
		
	if (j<4) {
		if (j==mu)
			value += helper->gamma * getISBeta() * s[nu];
		if (j==nu)
			value += helper->gamma * getISBeta() * s[mu];
			
		value += helper->gamma * helper->g[j] *( s[mu]*helper->Pixy[nu][j] + s[nu]*helper->Pixy[mu][j] );
		value -= (s[mu]*helper->Pixy[nu][0]+s[nu]*helper->Pixy[mu][0]) * helper->g[j]*s[j];
		
		if (ISRESCALE)
			value -= (4./3.)*helper->Pixy[mu][nu]*helper->g[j]*s[j]/helper->gamma;
		
		value += ((2.*getISBeta())/(3.*helper->gamma)) * helper->g[j] * s[j] * s[nu] * s[mu];
		
		if (mu == nu)
			value -= ((2.*getISBeta())/(3.*helper->gamma)) * helper->g[j] * s[j]/helper->g[mu];

		return value;
	}
	else { //(j==10)
		double temp;

		value -= helper->Pixy[mu][nu]/getTIS();
		value -= helper->gamma * helper->DPixy[0][mu][nu];
		
		for (int m=1;m<4;m++) {
			value -= s[m]*helper->DPixy[m][mu][nu];
			value += s[m]*helper->Pixy[mu][nu]/getISAlpha() * helper->dTISS[m-1];
		}
		
		for (int m=0;m<4;m++){
			value -= helper->g[m]*(s[mu]*helper->Pixy[nu][m] + s[nu]*helper->Pixy[mu][m])*helper->gamma*helper->dU[0][m];
			for (int n=1;n<4;n++) 
				value -= helper->g[m]*(s[mu]*helper->Pixy[nu][m]+s[nu]*helper->Pixy[mu][m])*s[n]*helper->dU[n][m];
		}
		
		if (ISRESCALE) {
			temp = (4./3.)*helper->Pixy[mu][nu];
			for (int m=0;m<4;m++)
				value -= temp*helper->dU[m][m];
		}
		
		value += getISBeta() * ((helper->dU[mu][nu]/helper->g[mu]) + (helper->dU[nu][mu]/helper->g[nu]));
		value -= getISBeta() * helper->gamma*(s[mu]*helper->dU[0][nu] + s[nu]*helper->dU[0][mu]);
		
		for (int m=1;m<4;m++) 
			value -= getISBeta() * s[m] * (s[mu]*helper->dU[m][nu] + s[nu]*helper->dU[m][mu]);
		
		if (mu==nu)
			temp = getISBeta() * (2./3.) * ((1./helper->g[mu])-s[mu]*s[mu]);
		else 
			temp = -getISBeta() * (2./3.)* s[mu]*s[nu];
		
		for (int m=0;m<4;m++)
			value -= temp * helper->dU[m][m];
		return value;
	}
}

	// bulk component
double CCell::getPiCoeff(int j) {
	
		//BULKTAG - 0.5 D (beta_0/T) missing
	if (j==0)
		return 0.;
	else if (j<3){
		if (ISRESCALE)
			return -((4./3.)*helper->BPi + 0.5*helper->BPi + getISGamma())*helper->g[j]*s[j]/helper->gamma;
		else 
			return -s[j]*helper->g[j]*(0.5*helper->BPi+getISGamma())/helper->gamma;
	}
	if (j<10)
		return 0.;
	else if (j==10)
		return helper->gamma*getISAlpha();
	else {
		double value = 0.;
		value -= helper->BPi/getTISB();
		for (int m=1;m<4;m++) 
			value -= s[m] * getISAlpha() * helper->dS[10][m-1];
		
		if (ISRESCALE) {
			double temp = (4./3.)*helper->BPi + 0.5*helper->BPi + getISGamma();
			for (int m=0;m<4;m++)
				value -= temp*helper->dU[m][m];
		}
		else 
			for (int m=0;m<4;m++) 
				value -= (0.5*helper->BPi+getISGamma())*helper->dU[m][m];
		return value;
	} 
}

void CCell::fillEosVar(){ 
	eos->fill(getE());
	
	if (mSVTrim)
		eos->trimSV(sqrt(x[1]*x[1]+x[2]*x[2]));
}

void CCell::initNS() {
	if (mInitNS != 0.) {
			// do what would normally happen in an CCell::update();
			// except ignore derivatives of shear tensor 
		
		fillEosVar();
		
		if (mBjorken)
			helper->g[3] = -getTau()*getTau();
		
		helper->gamma = sqrt(1.+s[1]*s[1]+s[2]*s[2]-helper->g[3]*s[3]*s[3]);
		helper->y = sqrt(s[5]*s[5] + s[6]*s[6] + s[5]*s[6]);
		helper->a = getAMax() * tanh( getISAlpha() * helper->y / getAMax());
		
		helper->BPi=s[10]*getISAlpha();
		
		fillPixy();
		
		for (int j=0;j<4;j++)
			for (int i=0;i<3;i++)
				if (neighbors[i][1]->getActive()) {
					if (ISSCALE == 'c' && j>3) {
						helper->dS[j][i] = getE()*(neighbors[i][1]->s[j]/neighbors[i][1]->getE() - neighbors[i][0]->s[j]/neighbors[i][0]->getE())/(2.*dx[i+1]);
						helper->dS[j][i] += s[j]*helper->dS[0][i];
					}
					else 
						helper->dS[j][i] = (neighbors[i][1]->s[j] - neighbors[i][0]->s[j])/ (2.*dx[i+1]);
				}
				else {
					if (ISSCALE == 'c' && j>3) {
						helper->dS[j][i] = (3.*s[j] + getE() * ( neighbors[i][0]->neighbors[i][0]->s[j]/neighbors[i][0]->neighbors[i][0]->getE() 
																- 4.*neighbors[i][0]->s[j]/neighbors[i][0]->getE())) / (2.*dx[i+1]);
						helper->dS[j][i] += s[j]*helper->dS[0][i];
					}
					else
						helper->dS[j][i] = (3.*s[j] - 4.*neighbors[i][0]->s[j] 
											+ neighbors[i][0]->neighbors[i][0]->s[j]) / (2.*dx[i+1]);
				}
		
		for (int i=0;i<3;i++) 
			helper->dTISS[i] = getDISAlphaDE() * getE() * helper->dS[0][i];
		
		
		if (mPureBjorken)
			for (int j=0;j<4;j++)
				helper->dS[j][2] = 0.;
		
		for (int j=4;j<11;j++)
			for (int i=0;i<3;i++) {
				helper->dS[j][i] = 0.;
			}
		
		fillDU();
		fillDPixy();
		fillM();
		
		if (fabs(x[1] - 5.) < .01 && fabs(x[2]) < .01 && fabs(x[3]) < 0.01 && false) {
			print();
			printM();
		}
		
		if (!mSVTrimInit) {
			for (int i=4;i<11;i++) 
					//			s[i] = (mInitNS*getTIS()/getISAlpha())*(helper->M[i][10] - helper->M[i][1]*dUdT[0] - helper->M[i][2]*dUdT[1]);			
				s[i] = (mInitNS*getTIS()/getISAlpha())*helper->M[i][11];
		}
		else {		
			for (int i=4;i<11;i++) 
				s[i] = (mInitNS*getTIS()/getISAlpha())*(helper->M[i][11]/(1.+exp((sqrt(x[1]*x[1]+x[2]*x[2]) - 10.)/0.6))
														- helper->M[i][1]*dUdT[1] - helper->M[i][2]*dUdT[2]);
		}	
		
		if (fabs(x[1]) < .01 && fabs(x[2]) < .01 && false) {
			fillPixy();
			print();
			printM();
		}
	}
	else { // mInitNSTxxP != 0.
		s[6] = -mInitNSTxxP*getTxyCalc(3,3)/getISAlphaCalc();
		s[4] = -0.5*s[6]*getTau()*getTau();
		s[5] = s[4];
		s[7] = 0.;
		s[8] = 0.;
		s[9] = 0.;
		s[10]= 0.;
	}
}

void CCell::addInitialFlow(){
	fillDS();
	fillEosVar();
	
	s[1] = -mInitFlow * 0.5 * helper->dS[0][0]/(1+getP()/getE()) * getTau();
	dUdT[0] = s[1]/getTau();

	s[2] = -mInitFlow * 0.5 * helper->dS[0][1]/(1+getP()/getE()) * getTau();
	dUdT[1] = s[2]/getTau();
}

	// turn this cell off; if impossible return false
	// FIXME :assumes octant=true
CCell* CCell::deactivate() {
	active = false;
	CCell* returnCell;
	int nCells=0;
	
	if ( mOctant == 3) {
		for (int i=0;i<3;i++) {
			if ( fabs(neighbors[i][0]->x[i+1]) < 0.0001 ) {
				returnCell = neighbors[i][0]->deactivate();
				if (returnCell == NULL)
					returnCell = neighbors[i][0];
				else if (returnCell == neighbors[i][1])
					nCells++;
				nCells++;
			}
		}
	}
	else {
		for (int i=0;i<3;i++){
			if (neighbors[i][0]->getActive())
				if (!neighbors[i][0]->neighbors[i][0]->getActive()){
					returnCell = neighbors[i][0]->deactivate();
					if (returnCell == NULL)
						returnCell = neighbors[i][0];
					else if (returnCell == neighbors[i][0])
						nCells++;
					nCells++;
				}
			if (neighbors[i][1]->getActive())
				if (!neighbors[i][1]->neighbors[i][1]->getActive()){
					returnCell = neighbors[i][1]->deactivate();
					if (returnCell == NULL)
						returnCell = neighbors[i][1];
					else if (returnCell == neighbors[i][1])
						nCells++;
					nCells++;
				}
		}
	}
	
	if (nCells == 0)
		return NULL;
	else if (nCells == 1)
		return returnCell;
	else 
		return this;
}

double CCell::getExpansionRate(){
	double value = 0.;
	for (int i=1;i<4;i++){
		value += helper->dU[i][i];
		value -= helper->g[i]*s[i]*dUdT[i-1]/helper->gamma;
	}
	return value;
}

double CCell::getTxy(int x, int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknown case of CCell::getTxy(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif

		//	double value = getPixy(x,y);
	double value = helper->Pixy[x][y];
	
	if (x==0)
		if (y==0)
			value += helper->gamma*helper->gamma*(getE()+getP()+helper->BPi);
		else
			value += helper->gamma*s[y]*(getE()+getP()+helper->BPi);
	else
		if (y==0)
			value += helper->gamma*s[x]*(getE()+getP()+helper->BPi);
		else 
			value += s[x]*s[y]*(getE()+getP()+helper->BPi);
		
	if (x==y)
		value -= (getP()+helper->BPi)/helper->g[x];
		
	return value;
}

double CCell::getPixy(int x, int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknown case of CCell::getPixy(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	
	return helper->Pixy[x][y];
	
	/*
	if (x==y){
			// diagonals
		if (mISMax) {
			return 0.; //FIXME
			if (helper->y==0.) return 0.;
			if (x==1) 
				return (helper->a/helper->y)*(s[5] + s[6]/ROOT3);
			else if (x==2) 
				return (helper->a/helper->y)*(-s[5] + s[6]/ROOT3);
			else 
				return - (helper->a/helper->y) * (2./ROOT3) * s[6];
		} 
		else {
			if (x==1) 
				return s[4]*getISAlpha();
			else if (x==2) 
				return s[5]*getISAlpha();
			else if (x==3)
				return s[6]*getISAlpha();
			else 
				return -(helper->g[1]*s[4]+helper->g[2]*s[5]+helper->g[3]*s[6])*getISAlpha();
		}
	}

	if (x>y){
		int temp = x;
		x=y;
		y=temp;
	}
	
	if (x==0){
		if (y==1) 
			return -(helper->g[1]*s[1]*s[4]+ helper->g[2]*s[2]*s[7] +helper->g[3]*s[3]*s[8])*getISAlpha()/helper->gamma;
		if (y==2) 
			return -(helper->g[1]*s[1]*s[7]+ helper->g[2]*s[2]*s[5] +helper->g[3]*s[3]*s[9])*getISAlpha()/helper->gamma;
		if (y==3) 
			return -(helper->g[1]*s[1]*s[8]+ helper->g[2]*s[2]*s[9] +helper->g[3]*s[3]*s[6])*getISAlpha()/helper->gamma;
	}
	
		// off-diagonals - not subject to maximums
	if (x==1) {
		if (y==2) 
			return s[7]*getISAlpha();
		else 
			return s[8]*getISAlpha();
	} 
	else 
		return s[9]*getISAlpha();
	*/
}

double CCell::getPixyBulk(int x, int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknown case of CCell::getPixyBulk(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	double value = helper->Pixy[x][y];
	value -= s[x]*s[y] * getB();
	if (x==y)
		value += helper->g[x]*getB();
	
	return value;
}

double CCell::getTxyLocal(int x, int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknown case of CCell::getTxyLocal(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	
	if (x==0 && y==0) return getE();
	
	if (x==0 || y==0) return 0.;
	
	double value = getPixyLocal(x,y);

	if (x==y) {
		if (x==3 && mBjorken) value += (getP()+helper->BPi)/pow(getTau(),2);
		else value += getP()+helper->BPi;
	}
	
	return value;
}

double CCell::getPixyLocal(int x, int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknown case of CCell::getPixyLocal(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	
	if (x==0 || y==0) return 0.;
	
	double value = helper->Pixy[x][y];
	value -= s[x]/(1+helper->gamma) * helper->Pixy[y][0];
	value -= s[y]/(1+helper->gamma) * helper->Pixy[x][0];
	value += s[x]*s[y]/pow(1+helper->gamma,2) * helper->Pixy[0][0];
	
	return value;
}

double CCell::getTxyCalc(int x, int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknown case of CCell::getTxy(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	
	double value = getPixyCalc(x,y);
		
	if (x==y)
		if (x==0) 
			value -= getPCalc() + getBPiCalc();
		else if (x<3 || !mBjorken)
			value += getPCalc() + getBPiCalc();
		else 
			value += (getPCalc() + getBPiCalc())/pow(getTau(),2);
		
	if (x>y){
		int temp = x;
		x=y;
		y=temp;
	}
	
	if (x==0){
		if (y==0) {
			if (mBjorken)
				value += (1.+s[1]*s[1]+s[2]*s[2]+getTau()*getTau()*s[3]*s[3])*(getE()+getPCalc()+getBPiCalc());
			else 
				value += (1.+s[1]*s[1]+s[2]*s[2]+s[3]*s[3])*(getE()+getPCalc()+getBPiCalc());
		}
		else{
			if (mBjorken)
				value += sqrt(1.+s[1]*s[1]+s[2]*s[2]+getTau()*getTau()*s[3]*s[3])*s[y]*(getE()+getPCalc()+getBPiCalc());
			else 
				value += sqrt(1.+s[1]*s[1]+s[2]*s[2]+s[3]*s[3])*s[y]*(getE()+getPCalc()+getBPiCalc());
		}
			
	}
	else
		value += s[x]*s[y]*(getE()+getPCalc()+getBPiCalc());
	
	
	return value;
}

double CCell::getPixyCalc(int x, int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknown case of CCell::getPixyCalc(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	
	if (x==y){
			// diagonals
		if (mISMax) {
			return 0.;  //FIXME
			double mY = sqrt(s[5]*s[5]+s[6]*s[6]);
			double mA = getAMaxCalc()*tanh(getISAlphaCalc()*mY/getAMaxCalc());

			if (mY==0.) return getPCalc();

			if (x==1) 
				return getPCalc() + (mA/mY)*(s[5] + s[6]/ROOT3);
			else if (x==2) 
				return getPCalc() + (mA/mY)*(-s[5] + s[6]/ROOT3);
			else 
				return getPCalc() - (mA/mY) * (2./ROOT3) * s[6];
		} 
		else {
			if (x==1) 
				return s[4]*getISAlphaCalc();
			else if (x==2) 
				return s[5]*getISAlphaCalc();
			else if (x==3)
				return s[6]*getISAlphaCalc();
			else {
				if (mBjorken)
					return (s[4]+s[5]+getTau()*getTau()*s[6])*getISAlphaCalc();
				else 
					return (s[4]+s[5]+s[6])*getISAlphaCalc();
			} 
		}
	}

	if (x>y){
		int temp = x;
		x=y;
		y=temp;
	}
	
	if (x==0){
		if (mBjorken){
			if (y==1) 
				return (s[1]*s[4]+ s[2]*s[7] + getTau()*getTau()*s[3]*s[8])*getISAlphaCalc()/sqrt(1.+s[1]*s[1]+s[2]*s[2]+getTau()*getTau()*s[3]*s[3]);
			if (y==2) 
				return (s[1]*s[7]+ s[2]*s[5] + getTau()*getTau()*s[3]*s[9])*getISAlphaCalc()/sqrt(1.+s[1]*s[1]+s[2]*s[2]+getTau()*getTau()*s[3]*s[3]);
			if (y==3) 
				return (s[1]*s[8]+ s[2]*s[9] + getTau()*getTau()*s[3]*s[6])*getISAlphaCalc()/sqrt(1.+s[1]*s[1]+s[2]*s[2]+getTau()*getTau()*s[3]*s[3]);
		}
		else {
			if (y==1) 
				return (s[1]*s[4]+ s[2]*s[7] + s[3]*s[8])*getISAlphaCalc()/sqrt(1.+s[1]*s[1]+s[2]*s[2]+s[3]*s[3]);
			if (y==2) 
				return (s[1]*s[7]+ s[2]*s[5] + s[3]*s[9])*getISAlphaCalc()/sqrt(1.+s[1]*s[1]+s[2]*s[2]+s[3]*s[3]);
			if (y==3) 
				return (s[1]*s[8]+ s[2]*s[9] + s[3]*s[6])*getISAlphaCalc()/sqrt(1.+s[1]*s[1]+s[2]*s[2]+s[3]*s[3]);
		}
	}
	
	
		// off-diagonals - not subject to maximums
	if (x==1) {
		if (y==2) 
			return s[7]*getISAlphaCalc();
		else 
			return s[8]*getISAlphaCalc();
	} 
	else 
		return s[9]*getISAlphaCalc();
}

double CCell::getTxyLocalCalc(int x, int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknwon case of CCell::getTxyMesh(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	
	if (x==0 && y==0) return getE();
	
	if (x==0 || y==0) return 0.;
	
	double value = getPixyLocalCalc(x,y);
	
	if (x==y) {
		if (x==3 && mBjorken) value += (getPCalc()+getBPiCalc())/pow(getTau(),2);
		else value += getPCalc() + getBPiCalc();
	}
	
	return value;	
}

double CCell::getPixyLocalCalc(int x, int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknwon case of CCell::getPixyLocalCalc(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	
	if (x==0 || y==0) return 0.;
	
	double value = getPixyCalc(x,y);
	double mG;
	if (mBjorken)
		mG = sqrt(1 + s[1]*s[1] + s[2]*s[2] + pow(getTau()*s[3],2));
	else 
		mG = sqrt(1 + s[1]*s[1] + s[2]*s[2] + s[3]*s[3]);
	
	value -= s[x]/(1+mG) * getPixyCalc(y,0);
	value -= s[y]/(1+mG) * getPixyCalc(x,0);
	value += s[x]*s[y]/pow(1+mG,2) * getPixyCalc(0,0);
	
	return value;
}

double CCell::getPixyNS(int x,int y) {
#ifdef DEBUG
	if (x<1 || x>3 || y<1 || y>3) {
		std::cout << "Unknwon case of CCell::getPixyNS(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	
	double value = 0.;
	value -= helper->gamma*s[x]*dUdT[y] + helper->gamma*s[y]*dUdT[x];
	/*
	for (int i=0;i<3;i++)
		if (x==y)
			value += (2./3.)*((1./helper->g[x]) - s[x]*s[y])*helper->g[i+1]*s[i+1]*dUdT[i]/helper->gamma;
		else 
			value -= (2./3.)*s[x]*s[y]*
*/
	
	return -getSV()*value;
}

double CCell::getPixyNSLocal(int x,int y) {
#ifdef DEBUG
	if (x<0 || x>3 || y<0 || y>3) {
		std::cout << "Unknwon case of CCell::getPixyNSMesh(" << x << ',' << y << "):" << std::endl;
		return 0.;
	}
#endif
	
	return 0.;
}

double CCell::getTZR() {
	return 0.;
}

double CCell::getTZPhi() {
	return 0.;
}

double CCell::getTRR() {
	return 0.;
}

double CCell::getTRPhi() {
	return 0.;
}

double CCell::getTPhiPhi() {
	return 0.;
}

double CCell::getPiZR() {
	return 0.;
}

double CCell::getPiZPhi() {
	return 0.;
}

double CCell::getPiRR() {
	return 1.;
}

double CCell::getPiRPhi() {
	return 0.;
}

double CCell::getPiPhiPhi() {
	return 0.;
}

double CCell::getVorticity(int i, int j) {
#ifdef DEBUG
	if (j<1 || j>3 || i<1 || i>3) {
		std::cout << "Unknown case of CCell::getVorticity(" << i << ',' << j << ',' << k << "):" << std::endl;
		return 0.;
	}
#endif
	
	double value = helper->dU[i][j]/helper->g[i] - helper->dU[j][i]/helper->g[j];
	value += s[j] * helper->gamma * dUdT[i-1] - s[i] * helper->gamma * dUdT[j-1];
	
	for (int k=1;k<4;k++) {
		value += s[j] * s[k] * helper->dU[k][i] - s[i] * s[k] * helper->dU[k][j];
	}
	
	return 0.5 * value;
}

void CCell::setU(double x0, double y0, double z0) {
	s[1] = x0;
	s[2] = y0;
	s[3] = z0;
}

void CCell::setS(int i, double v) {
#ifdef DEBUG
	if (i<0 || i>10) {
		std::cout << "Unknown case of setS("<< i << ',' << v << ")!" << std::endl;
		return;
	}
#endif
/*
	if (i==0) 
		s[i]=log(v);
	else 
*/
	s[i] = v;
}

	// Integrates this cell forward in time
	// NOTE: update() must be called for all cells before forward
void CCell::forward(CCell* mCell) {
		// fill the mixed partials and M_ij
	update();

	// NOTE:  commented lines are the intent of a particular row operation
	// these elements may be non-zero in memory, but would be zero
	// upon completion of the implied row operation.
	
			// remove tau_ij dep from cons eqs                                                                             
	for (int k=1;k<4;k++) 
		for (int j=4;j<11;j++) {
			for (int i=0; i<j;i++) 
				helper->M[k][i] -= (helper->M[k][j]/helper->M[j][j])*helper->M[j][i];
			helper->M[k][11] -= (helper->M[k][j]/helper->M[j][j])*helper->M[j][11];
			helper->M[k][j] = 0.;
		}
		
		// finish LT form
	for (int k=3;k>0;k--)
		for (int j=k-1;j>=0;j--) {
			helper->M[j][11] -= (helper->M[j][k]/helper->M[k][k])*helper->M[k][11];
			for (int i=0;i<k;i++)  
				helper->M[j][i] -= (helper->M[j][k]/helper->M[k][k])*helper->M[k][i];
			helper->M[j][k] = 0.;
		}
		
			// diagonalize M_ii, calculate change to vector side (M_j_11)                                 
	for (int j=1;j<4;j++) {
		helper->M[j][11] -= (helper->M[j][0]/helper->M[0][0])*helper->M[0][11];	
		helper->M[j][0] = 0.;
	}
		
	for (int i=1;i<4;i++)
		for (int j=i+1;j<12;j++) {
			helper->M[j][11] -= (helper->M[j][i]/helper->M[i][i])*helper->M[i][11];
			helper->M[j][i] = 0.;
		}

	if (x[1] == 0. && x[2] == 0. && x[3] == 0. && false){
			/*		printM();
		printf("%0.9g %0.9g %0.9g %0.9g\n\n",s[0],
			   dx[0] * helper->M[0][11]/helper->M[0][0],s[0] + dx[0] * helper->M[0][11]/helper->M[0][0],
			   exp(s[0] + dx[0] * helper->M[0][11]/helper->M[0][0]));
			 */
		printf("%0.9g %0.9g %0.9g %0.9g \n",helper->M[0][11]/helper->M[0][0],helper->M[1][11]/helper->M[1][1],
			   helper->M[4][11]/helper->M[4][4],helper->M[6][11]/helper->M[6][6]);
	}
	
	
	
	if (mLinT) {
		for (int i=0;i<11;i++) 
			mCell->setS(i, s[i] + dx[0] * helper->M[i][11]/helper->M[i][i]);
		mCell->setTau( x[0] + dx[0]);
	}
	else if (mLogT) {
		for (int i=0;i<11;i++) 
			mCell->setS(i, s[i] +  mT0*exp(x[0]) * dx[0] * helper->M[i][11]/helper->M[i][i]);	  
		mCell->setTau(mT0*exp(x[0] + dx[0]));
	}
	else if (mLogSinhT) {
		for (int i=0;i<11;i++) {
			if (mLaxFried) {
				double averageS = 0.;
				double avDenom = 0.;
				
				if (!mPureBjorken){
					for (int m=0;m<3;m++) 
						if (neighbors[m][1]->getActive()) {
							for (int n=0;n<2;n++) 
								averageS += neighbors[m][n]->getS(i);
							avDenom += 2.;
						}			
				}
				else 
					for (int m=0;m<2;m++) 
						if (neighbors[m][1]->getActive()) {
							for (int n=0;n<2;n++) 
								averageS += neighbors[m][n]->getS(i);
							avDenom += 2.;
						}
				
				if (avDenom == 0.) {
					averageS = s[i];
					avDenom = 1.;
				}
				
					//				printf("(%g %g) - %g/%g = %g [%g]\n",x[1],x[2],averageS,avDenom,averageS/avDenom,s[i]);
				
				mCell->setS(i, averageS/avDenom + tanh(getTau())*dx[0]*helper->M[i][10]/helper->M[i][i]);
			} 
			else mCell->setS(i, s[i] + tanh(getTau())*dx[0]*helper->M[i][11]/helper->M[i][i]);
			
			
		}
		mCell->setTau(asinh(exp(x[0] + dx[0])));
	}
	
	
//	printf("\n%g %g %g %g %g\n\n",helper->M[0][10]/helper->M[0][0],helper->M[1][10]/helper->M[1][1],
//			helper->M[4][10]/helper->M[4][4],helper->M[5][10]/helper->M[5][5],helper->M[6][10]/helper->M[6][6]);
	
	for (int i=1;i<4;i++) 
		mCell->dUdT[i-1] = (helper->M[i][11]/helper->M[i][i]);
}

	// Integrates this cell forward in time
void CCell::forward(CCell* onCell, CCell* offCell) {
		// fill the mixed partials and M_ij
	offCell->update();
	
		//	CCellHelper *mH = helper;
		//	helper = offCell->getHelper();
	
	// NOTE:  commented lines are the intent of a particular row operation
	// these elements may be non-zero in memory, but would be zero
	// upon completion of the implied row operation.
	
			// remove tau_ij dep from cons eqs                                                                             
	for (int k=1;k<4;k++) 
		for (int j=4;j<11;j++) {
			for (int i=0; i<j;i++) 
				helper->M[k][i] -= (helper->M[k][j]/helper->M[j][j])*helper->M[j][i];
			helper->M[k][11] -= (helper->M[k][j]/helper->M[j][j])*helper->M[j][11];
//			helper->M[k][j] = 0.;
		}
		
		// finish LT form
	for (int k=3;k>0;k--)
		for (int j=k-1;j>=0;j--) {
			helper->M[j][11] -= (helper->M[j][k]/helper->M[k][k])*helper->M[k][11];
			for (int i=0;i<k;i++)  
				helper->M[j][i] -= (helper->M[j][k]/helper->M[k][k])*helper->M[k][i];
//			helper->M[j][k] = 0.;
		}
		
			// diagonalize M_ii, calculate change to vector side (M_j_11)                                 
	for (int j=1;j<4;j++) {
		helper->M[j][11] -= (helper->M[j][0]/helper->M[0][0])*helper->M[0][11];	
//		helper->M[j][0] = 0.;
	}
		
	for (int i=1;i<4;i++)
		for (int j=i+1;j<12;j++) {
			helper->M[j][11] -= (helper->M[j][i]/helper->M[i][i])*helper->M[i][11];
//			helper->M[j][i] = 0.;
		}
	
		// note - off-diagonal elements HAVE NOT been set to zero, but HAVE been taken in to account
	if (mLinT) {
		for (int i=0;i<11;i++) 
			onCell->setS(i, s[i] + dx[0]*helper->M[i][11]/helper->M[i][i]);
		onCell->setTau( x[0] + dx[0]);
	}
	else if (mLogT) {
		for (int i=0;i<11;i++) 
			onCell->setS(i, s[i] + mT0*exp(offCell->x[0]) * dx[0] * helper->M[i][11]/helper->M[i][i]);
		onCell->setTau( mT0 * exp( x[0] + dx[0]));
	} 
	else if (mLogSinhT) {
		for (int i=0;i<11;i++)
			onCell->setS(i, s[i] + dx[0] * tanh(offCell->getTau()) * helper->M[i][11]/helper->M[i][i]);
		onCell->setTau(asinh(exp(x[0] + dx[0])));
	}	
	
	for (int i=1;i<4;i++) 
		onCell->dUdT[i-1] = helper->M[i][11]/helper->M[i][i];
	
		//	helper = mH;
}

void CCell::forward(CCell* k0, CCell* k1, CCell* k2, CCell* k3){
		// since the matrix for RK4 has contributions from 4 different
		// time steps, we'll compile it here from the usual update()
		// of the four cells to be used.
	double mM[10][11];
	
	k0->update();
	for (int i=0;i<10;i++) {
		for (int j=0;j<10;j++) mM[i][j] = (helper->M[i][j]/6.);
		if (mLogT)
			if (i!=3) mM[i][10] = helper->M[i][10]*exp(k0->x[0])/6.;
			else mM[i][10] = helper->M[i][10]*exp(-k0->s[4] + k0->x[0])/6.;
			else if (mLogSinhT)
				if (i!=3) mM[i][10] = helper->M[i][10]*tanh(k0->getTau())/6.;
				else mM[i][10] = helper->M[i][10]*exp(-k0->s[4])*tanh(k0->getTau())/6.;
				else if (mLinT)
					if (i!=3) mM[i][10] = helper->M[i][10]/6.;
					else mM[i][10] = helper->M[i][10]*exp(-k0->s[4])/6.;
	} 
	
	k3->update();
	for (int i=0;i<10;i++) {
		for (int j=0;j<10;j++) mM[i][j] += (helper->M[i][j]/6.);
		if (mLogT)
			if (i!=3) mM[i][10] += (helper->M[i][10]*exp(k3->x[0])/6.);
			else mM[i][10] += (helper->M[i][10]*exp(-k3->s[4] + k3->x[0])/6.);
			else if (mLogSinhT)
				if (i!=3) mM[i][10] += (helper->M[i][10]*tanh(k3->getTau())/6.);
				else mM[i][10] += (helper->M[i][10]*exp(-k3->s[4])*tanh(k3->getTau())/6.);
				else if (mLinT)
					if (i!=3) mM[i][10] = helper->M[i][10]/6.;
					else mM[i][10] = helper->M[i][10]*exp(-k3->s[4])/6.;
	}
	
	k1->update();
	for (int i=0;i<10;i++) {
		for (int j=0;j<10;j++) mM[i][j] += helper->M[i][j]/3.;
		if (mLogT)
			if (i!=3)  mM[i][10] += helper->M[i][10]*exp(k1->x[0])/3.;
			else mM[i][10] += helper->M[i][10]*exp(-k1->s[4] + k1->x[0])/3.;
			else if (mLogSinhT)
				if (i!=3)  mM[i][10] += helper->M[i][10]*tanh(k1->getTau())/3.;
				else mM[i][10] += helper->M[i][10]*exp(-k1->s[4])*tanh(k1->getTau())/3.;
				else if (mLinT)
					if (i!=3) mM[i][10] = helper->M[i][10]/6.;
					else mM[i][10] = helper->M[i][10]*exp(-k1->s[4])/6.;
	}
	
	k2->update();
	for (int i=0;i<10;i++) {
		for (int j=0;j<10;j++) mM[i][j] += helper->M[i][j]/3.;
		if (mLogT)
			if (i!=3) mM[i][10] += helper->M[i][10]*exp( k2->x[0])/3.;
			else mM[i][10] += helper->M[i][10]*exp(-k2->s[4] + k2->x[0])/3.;
			else if (mLogSinhT)
				if (i!=3) mM[i][10] += helper->M[i][10]*tanh(k2->getTau())/3.;
				else mM[i][10] += helper->M[i][10]*exp(-k2->s[4])*tanh(k2->getTau())/3.;
				else if (mLinT)
					if (i!=3) mM[i][10] = helper->M[i][10]/6.;
					else mM[i][10] = helper->M[i][10]*exp(-k2->s[4])/3.;
	}
	
		// remove tau_ij dep from cons eqs                                                                             
	for (int k=0;k<4;k++) {
		for (int j=4;j<10;j++) {
			for (int i=0; i<j;i++) mM[k][i] -= (mM[k][j]/mM[j][j])*mM[j][i];
			mM[k][10] -= (mM[k][j]/mM[j][j])*mM[j][10];
		}
	}
		// finish LT form                                                                                              
	for (int k=3;k>0;k--)
		for (int j=k-1;j>=0;j--) {
			for (int i=0;i<k;i++)  mM[j][i] -= (mM[j][k]/mM[k][k])*mM[k][i];
			mM[j][10] -= (mM[j][k]/mM[k][k])*mM[k][10];
			mM[j][k] = 0.;
		}
		// diagonal M_ii is now as it will be, calculate change to vector side (M_j_10)                                
	for (int i=0;i<=3;i++)
		for (int j=i+1;j<11;j++) 
			mM[j][10] -= (mM[j][i]/mM[i][i])*mM[i][10];
	
	double mV[3];
	for (int i=0;i<3;i++) 
		if (mLogT) mV[i] = mT0 * dx[0] * mM[i][10]/mM[i][i];
		else      mV[i] = dx[0] * mM[i][10]/mM[i][i];
	
		//Now move it forward
	double uv = s[1]*mV[0] + s[2]*mV[1] + s[3]*mV[2];
	for (int i=0;i<3;i++)
		s[i+1] =  mV[i] + s[i+1]*sqrt(mV[0]*mV[0]+mV[1]*mV[1]+mV[2]*mV[2]) + (uv/(1.+s[0])) * s[i+1];
	
		//Now move it forward
	if (mLogT) 
		for (int i=4;i<10;i++) 
			s[i+1] += mT0 * dx[0] * mM[i][10]/mM[i][i];
	else if (mLogSinhT || mLinT)
		for (int i=4;i<10;i++)
			s[i+1] += (dx[0]) * mM[i][10]/mM[i][i];
	
	s[0] = sqrt(1. + s[1]*s[1] + s[2]*s[2] + s[3]*s[3]);
	x[0] += dx[0];
}

void CCell::copy(CCell* mC) {
	for (int i=0;i<4;i++)
		x[i] = mC->x[i];
	
	for (int i=0;i<11;i++)
		s[i] = mC->s[i];
	
	active = mC->active;
}

void CCell::copyST(CCell* mC) {
	x[0] = mC->x[0];
	
	for (int i=0;i<11;i++)
		s[i] = mC->s[i];
	
	active = mC->active;
}

void CCell::flip(int i) {
	s[i] = -s[i];
	
	if (i==1) {
		s[7] = -s[7];
		s[8] = -s[8];
	}
	else if (i==2) {
		s[7] = -s[7];
		s[9] = -s[9];
	} 
	else {
		s[8] = -s[8];
		s[9] = -s[9];
	}
}

void CCell::print() { 
	if (active)
		std::cout << "Active ";
	else 
		std::cout << "Dead ";
	
	std::cout << "Cell (" << getTau() << ','
	<< x[1] << ',' << x[2] << ',' << x[3] << "):" << std::endl;
	
	std::cout << "dx (" << dx[0] << ','
	<< dx[1] << ',' << dx[2] << ',' << dx[3] << "):" << std::endl;
	
	std::cout << "g (" << helper->g[0] << ','
	<< helper->g[1] << ',' << helper->g[2] << ',' << helper->g[3] << "):" << std::endl;
	
	for (int i=0;i<11;i++){
		if (i!=0) std::cout << s[i] << '\t';
		else std::cout << getE() << '\t';
		if (i==3) std::cout << std::endl << '\t';
	}
	std::cout << std::endl;
	
	std::cout << std::endl;
	std::cout << "gamma = " << helper->gamma << std::endl;
	std::cout << "sv = " << getSV() << std::endl;
	std::cout << "bv = " << getBV() << std::endl;
	std::cout << "TIS = " << getTIS() << std::endl;
	std::cout << "ISAlpha = " << getISAlpha() << std::endl;
	std::cout << "ISBeta = " << getISBeta() << std::endl;
	std::cout << "ISGamma = " << getISGamma() << std::endl;
	std::cout << "pres = " << getP() << std::endl;
	std::cout << "T = " << getT() << std::endl;
	std::cout << "S = " << getS() << std::endl;
	std::cout << "dTISS:  ";
	for (int i=0;i<3;i++) std::cout << helper->dTISS[i] << "  ";
	
	std::cout << std::endl << "dUdT "<< dUdT[0] << " " << dUdT[1] << " " << dUdT[2] << std::endl << std::endl;
	
	std::cout << std::endl << "expRate:  " << getExpansionRate() << std::endl;
	
	std::cout << std::endl << "getDISAlphaDE():  " << getDISAlphaDE() << std::endl;
	std::cout << "cs2 = " << getCs2() << std::endl;
	
	for (int i=0;i<3;i++) {
		for (int j=0;j<11;j++) 
			if (fabs(helper->dS[j][i]) > 0.) 
				std::cout << "dS[" << j << "][" << i << "] = " << helper->dS[j][i] << std::endl;
		std::cout << std::endl;
	}
	
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++) 
			if (fabs(helper->dU[j][i]) > 0.) 
				std::cout << "dU[" << j << "][" << i << "] = " << helper->dU[j][i] << std::endl;
		std::cout << std::endl;
	}
	
	if (mBjorken) {
		std::cout << "Txy:" << std::endl;
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) 
				if (i!=3 && j!=3)
					std::cout << getTxy(i,j) << '\t';
				else if ((i==3 && !(j==3)) || (j==3 && !(i==3)))
					std::cout << getTau() * getTxy(i,j) << '\t';
				else 
					std::cout << getTau() * getTau() * getTxy(i,j) << '\t';
			std::cout << std::endl;
		}
	}
	else {
		std::cout << "Txy:" << std::endl;
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) 
				std::cout << getTxy(i,j) << '\t';
			std::cout << std::endl;
		}
	}
	
	if (mSVRatio == 0.0) return;
	
	if (mBjorken) {
		std::cout << "Pixy:" << std::endl;
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) 
				if (i!=3 && j!=3)
					std::cout << helper->Pixy[i][j] << '\t';
				else if ((i==3 && !(j==3)) || (j==3 && !(i==3)))
					std::cout << getTau() * helper->Pixy[i][j] << '\t';
				else 
					std::cout << getTau() * getTau() * helper->Pixy[i][j] << '\t';
			std::cout << std::endl;
		}
	}
	else {
		std::cout << "Pixy:" << std::endl;
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) 
				std::cout << helper->Pixy[i][j] << '\t';
			std::cout << std::endl;
		}
	}
	
	std::cout << std::endl << "DPixy:" << std::endl;
	for (int i=0;i<4;i++){
		for (int j=0;j<4;j++){
			for (int k=0;k<4;k++)
				printf("%g ",helper->DPixy[i][j][k]);
			printf("\n");
		}
		printf("\n");
	}
	printf("\n Neighbors active?? : ");

	for (int i=0;i<3;i++)
		for (int j=0;j<2;j++)
			printf("%d ",neighbors[i][j]->getActive());
	printf("\n\n");
	
	return;
}

void CCell::selfPrint() {
	if (active)
		std::cout << "Active ";
	else 
		std::cout << "Dead ";
	
	std::cout << "Cell (" << getTau() << ','
	<< x[1] << ',' << x[2] << ',' << x[3] << "):" << std::endl;
	
	std::cout << getE() << '\t';
	for (int i=1;i<11;i++)
		std::cout << s[i] << '\t';
	std::cout << std::endl;
	
	std::cout << "Neighbors active?? : ";
	
	for (int i=0;i<3;i++)
		for (int j=0;j<2;j++)
			std:: cout << neighbors[i][j]->getActive() << " ";
	std::cout << std::endl << std::endl;
	
	return;
	
	std::cout << "gamma = " << helper->gamma << std::endl;
	
	if (mBjorken) {
		std::cout << "Pixy:" << std::endl;
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) 
				if (i!=3 && j!=3)
					std::cout << helper->Pixy[i][j] << '\t';
				else if ((i==3 && !(j==3)) || (j==3 && !(i==3)))
					std::cout << getTau() * helper->Pixy[i][j] << '\t';
				else 
					std::cout << getTau() * getTau() * helper->Pixy[i][j] << '\t';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	else {
		std::cout << "Pixy:" << std::endl;
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) 
				std::cout << helper->Pixy[i][j] << '\t';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	
}

void CCell::verify() {
	printf("CCell::verify() - %g %g %g %g\n",getTau(),x[1],x[2],x[3]);
	double uPi[4];
	for (int i=0;i<4;i++) 
		uPi[i] = helper->gamma * helper->Pixy[0][i];
	
	for (int i=1;i<4;i++)
		for (int j=0;j<4;j++)
			uPi[j] += helper->g[i] * s[i] * helper->Pixy[i][j];
	
	printf("uPi[] = %g %g %g %g\n",uPi[0],uPi[1],uPi[2],uPi[3]);
	printf("Trace = %g\n",helper->Pixy[0][0]-helper->Pixy[1][1]-helper->Pixy[2][2]-getTau()*getTau()*helper->Pixy[3][3]);
}

	// prints the matrix inverted in forward()                                                                       
	// should only be called for debugging purposes...                                                               
void CCell::printM(){

	std::cout << "Matrix for Cell(" << getTau() << ',' << x[1] << ',' << x[2] << ',' << x[3] << ")" << std::endl;
	if (mSVRatio != 0. || mBVRatio != 0.)
		for (int i=0;i<11;i++) {
			for (int j=0;j<12;j++) {
				if (j==0)
					std::cout << helper->M[i][j]/getE() << '\t';
				else std::cout << helper->M[i][j] << '\t';
			} 
			std::cout << std::endl << std::endl;
		}
	else
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) std::cout << helper->M[i][j] << '\t';
			std::cout << helper->M[i][11] << std::endl << std::endl;
		}
}

void CCell::printMm(double mM[11][12]) {
	std::cout << "mMatrix for Cell(" << getTau() << ',' << x[1] << ',' << x[2] << ',' << x[3] << ")" << std::endl;
	
	if (mSVRatio != 0. || mBVRatio != 0.) 
		for (int i=0;i<11;i++) {
			for (int j=0;j<12;j++)
				std::cout << mM[i][j] << '\t' << '\t';
			std::cout << std::endl << std::endl;
		}
	else 
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) std::cout << mM[i][j] << '\t';
			std::cout << mM[i][11] << std::endl << std::endl;
		}
}

void CCell::getM(double mM[11][12]) {
	for (int i=0;i<11;i++)
		for (int j=0;j<12;j++)
			mM[i][j] = helper->M[i][j];
}

bool operator!= (CCell &c1, CCell &c2){
	return c1.x[1]!=c2.x[1] || c1.x[2]!=c2.x[2] || c1.x[3]!=c2.x[3];
}

	//statics
#ifndef HYDRO_BOOST_THREADS
CEos*  CCell::eos = NULL;
CCellHelper* CCell::helper = NULL;
#endif
bool CCell::mDebug, CCell::mSVTrim, CCell::mSVTrimInit;
bool CCell::mViscNS, CCell::mPureBjorken, CCell::mBjorken;
bool CCell::mLinT, CCell::mLogT, CCell::mLogSinhT;
bool CCell::mISVort, CCell::mISMax, CCell::mLaxFried, CCell::mSlopeLimit;
double CCell::mT0, CCell::mSVRatio, CCell::mBVRatio;
double CCell::mISAMax, CCell::mISBMax, CCell::mInitFlow;
double CCell::mInitNS, CCell::mInitNSTxxP, CCell::mSlopeLimitTheta;
double CCell::dx[4];
int CCell::mOctant;
parameterMap* CCell::pMap;
