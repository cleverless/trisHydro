/*  Cmesh.cpp
 *  isHydro3
 *  Created by Joshua Vredevoogd on 2/11/09.
 */

#include "CMesh.h"

	// Assumes that this constructor is only called once 
	// or that it is totally cleaned up between calls.
CMesh::CMesh(parameterMap* pM) {
	pMap = pM;
	fillParams();
	
	CCell* dummyCell = new CCell(pMap);
	
	if (mInitNS != 0.0 && mInitNSTxxP != 0.0) {
		std::cout << "\nCMesh::CMesh(parameterMap* pM)....\nConflicting instructions for shear tensor initialization...\n";
		std::cout << "Use either HYDRO_INIT_NS or HYDRO_INIT_NS_TXXP but not both!\n  Aborting.\n\n";
		fflush(stdout);
		exit(1);
	}
	
#ifdef HYDRO_BOOST_THREADS
	for (int i=0;i<NTHREADS;i++){
		eosVector[i] = new CEos();
		helperVector[i] = new CCellHelper();
	}
#endif
	
	wnNormalize();
	
	localE4 = new double***[2];
	for (int i1=0; i1 < 2; i1++) {
		localE4[i1] = new double**[2];
		for (int i2=0; i2 < 2; i2++) {
			localE4[i1][i2] = new double*[2];
			for (int i3=0; i3 < 2; i3++) {
				localE4[i1][i2][i3] = new double[2];
			}
		}
	}
	
	localE3 = new double**[2];
	for (int i1=0; i1 < 2; i1++) {
		localE3[i1] = new double*[2];
		for (int i2=0; i2 < 2; i2++) 
			localE3[i1][i2] = new double[2];
	}
	
	dxCorn = new double[4];
	
	delete dummyCell;
}

CMesh::~CMesh() {
	activeCells.clear();
}

void CMesh::fillParams() {
	mE0 = parameter::getD(*pMap,"HYDRO_E0",1.0);
	
	mDataRoot = parameter::getS(*pMap,"HYDRO_OUTPUT_DATAROOT","./");
	
	mOctant = parameter::getI(*pMap,"HYDRO_OCTANT",3);
	mPureBjorken = parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false);
	mBjorken = parameter::getB(*pMap,"HYDRO_BJORKEN",true);
	mPrintMs = parameter::getB(*pMap,"HYDRO_PRINTMS",false);
	mSVTrimInit = parameter::getB(*pMap,"HYDRO_SVTRIMINIT",false);
	mSVTrim = parameter::getB(*pMap,"HYDRO_SVTRIM",false);
	mFOTemp = parameter::getD(*pMap,"HYDRO_FOTEMP",0.165);
	mDeadT = parameter::getD(*pMap,"HYDRO_DEADT",0.025);
	
	mNSize = parameter::getI(*pMap,"HYDRO_NSIZE",20);
	mXSize = parameter::getI(*pMap,"HYDRO_XSIZE",60);
	mYSize = parameter::getI(*pMap,"HYDRO_YSIZE",60);
	mNSizeOrig = mNSize;
	mXSizeOrig = mXSize;
	mYSizeOrig = mYSize;
	
	mPrintN = parameter::getI(*pMap,"HYDRO_PRINTN",0);
	mPrintX = parameter::getI(*pMap,"HYDRO_PRINTX",0);
	mPrintY = parameter::getI(*pMap,"HYDRO_PRINTY",0);
	
	mT0 = parameter::getD(*pMap,"HYDRO_T0",1.0);
	mDn = parameter::getD(*pMap,"HYDRO_DN",0.1);
	mDx = parameter::getD(*pMap,"HYDRO_DX",0.1);
	mDy = parameter::getD(*pMap,"HYDRO_DY",0.1);
	
	wnRho0 = parameter::getD(*pMap,"GLAUBER_Rho0",0.17);
	wnRAu  = parameter::getD(*pMap,"GLAUBER_RAu",6.37);
	wnXi   = parameter::getD(*pMap,"GLAUBER_Xi",0.54);
	wnSigma= parameter::getD(*pMap,"GLAUBER_Sigma",4.2);
	wnA    = parameter::getD(*pMap,"GLAUBER_A",197.);
	wnB    = parameter::getD(*pMap,"GLAUBER_B",0.);
	wnKTau = parameter::getD(*pMap,"GLAUBER_K_TAU",4.3);
	mRn = parameter::getD(*pMap,"HYDRO_RN",1.4);
	
	bBGK   = parameter::getB(*pMap,"GLAUBER_BGK",false);
	if (bBGK && (mPureBjorken || mOctant!=0)){
		std::cout << "BGK only makes sense for longitudinal hydro." << std::endl;
		exit(1);
	}
	
	bGSat  = parameter::getB(*pMap,"GLAUBER_SAT",false);
	gSatN  = parameter::getD(*pMap,"GLAUBER_SAT_N",2.);
	mWnBinRatio = parameter::getD(*pMap,"GLAUBER_WNBIN_RATIO",1.0);
	
	fracPPE = parameter::getD(*pMap,"GLAUBER_PP_ENERGY_FRAC",1.0);
	sigmaAtt = parameter::getD(*pMap,"GLAUBER_SIGMA_ATT",wnSigma);
	wnAttRatio = parameter::getD(*pMap,"GLAUBER_WN_ATT_RATIO",1.0);
	
	collRootS = parameter::getD(*pMap,"GLAUBER_ROOT_S",200.);
	
	mInitFlow = parameter::getD(*pMap,"HYDRO_INIT_FLOW",0.0);
	mInitNS = parameter::getD(*pMap,"HYDRO_INIT_NS",0.0);
	mInitNSTxxP = parameter::getD(*pMap,"HYDRO_INIT_NS_TXXP",0.0);
	
	icNiceA = parameter::getD(*pMap, "IC_NICE_A", 1.4/ROOT2);
	icNiceB = parameter::getD(*pMap, "IC_NICE_B", 2.);
}

double CMesh::wnT(double x,double y) {
	wnTx=x;
	wnTy=y; 
	meshMemFunc mf = &CMesh::wnTIntegrand; 
	return 2.*qromb(mf,0.,15.);
}

void CMesh::wnNormalize() {
	meshMemFunc mf = &CMesh::wnRhoIntegrand; 
	wnRho0 *= wnA/(4.*JPI*qromb(mf,0.,15.));
}

double CMesh::wnE(double x, double y) {
	double value1 = wnT(x+wnB/2.,y);
	double value2 = wnT(x-wnB/2.,y);
	
	if (bGSat) {
		double mV1 = gSatN*value1*value2/(value1 + value2); 
		double mV2 = (value1 + value2)/gSatN; 

			// 2.613 GeV is average transverse energy per rapidity
			// in 200 GeV pp collisions
		return 2.613 * fracPPE * (wnSigma/sigmaAtt) *
		( 2. * (1. - wnAttRatio) * mV1 * (1.-exp(-sigmaAtt*mV2))
		 + wnAttRatio * (  value1*(1.-exp(-sigmaAtt*value2)) 
						 + value2*(1.-exp(-sigmaAtt*value1))));
	}
	
	return wnKTau*0.5*mWnBinRatio*(  value1 * ( 1. - pow(( 1. - (wnSigma/wnA)*value2),wnA))
								   + value2 * ( 1. - pow(( 1. - (wnSigma/wnA)*value1),wnA)))
	     + wnKTau *(1. - mWnBinRatio) * wnSigma * value1 * value2;
}

double CMesh::qromb(meshMemFunc mf, double a, double b) {
	const int JMAX=20, JMAXP=JMAX+1, K=5;
	const double EPS=1E-10;
	
	double ss, dss;
	std::vector<double> s, h, s_t, h_t;
	int i,j;
	
	h.push_back(1.0);
	for (j=1;j<=JMAX;j++){
		s.push_back(trapzd(mf,a,b,j));
		
		if (j>=K){
			for (i=0;i<K;i++){
				h_t.push_back(h[j-K+i]);
				s_t.push_back(s[j-K+i]);
			}
			polint(h_t,s_t,0.0,ss,dss);
			if (fabs(dss) <= EPS*fabs(ss)) {
				return ss;
			}
			
			h_t.clear();
			s_t.clear();
		}
		h.push_back(0.25*h[j-1]);
	} 
	
	fprintf(stderr,"Too many steps in routine qromb\n");
	exit(1);
	return 0.0;
}

double CMesh::trapzd(meshMemFunc mf, const double a, const double b, const int n){
	double x,tnm,sum,del;
	static double s;
	int it,j;
	
	if (n==1){
		return (s=0.5*(b-a)*((this->*(mf))(a) + (this->*(mf))(b)));
	}
	else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=0;j<it;j++,x+=del) 
			sum += (this->*(mf))(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

void CMesh::polint(std::vector<double> &xa, std::vector<double> &ya, const double x, double &y, double &dy){
	int i,m,ns=0;
	double den,dif,dift,ho,hp,w;
	
	int n=xa.size();
	std::vector<double> c(n,0.), d(n,0.);
	dif = fabs(x-xa[0]);
	for (i=0;i<n;i++){
		if ((dift=fabs(x-xa[i])) < dif){
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	y=ya[ns--];
	for (m=1;m<n;m++){
		for (i=0;i<n-m;i++){
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ((den=ho-hp) == 0.0) 
				fprintf(stderr,"Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
	}
}

	// assigns the cells initial conditions
void CMesh::initialCondition(CCell* mCell) {
		// if this should be non-zero due to initial flow (et cetera),
		// it will be taken care of later.
	mCell->setU( 0., 0., 0.);
	
	for (int l=4;l<11;l++) 
		mCell->setS(l,0.);
	
	mCell->setE( wnE(mCell->getX(), mCell->getY())/mCell->getTau());
	
}

	// assigns the cells initial conditions
void CMesh::initialCondition(CCell* mCell, CCell* zCell) {
	
		// if this should be non-zero due to initial flow (et cetera),
		// it will be taken care of later.
	mCell->setU(0., 0., 0.);
	
	for (int l=4;l<11;l++) 
		mCell->setS(l,0.);
	
	mCell->setE( zCell->getE() * exp(-.5*pow( (mCell->getEta()/mRn),2)));
	
}

void CMesh::addInitialFlow() {
	flipEdge();
	list<CCell*>::iterator it, end;
	for ( it=activeCells.begin(), end=activeCells.end(); it != end; ++it){
#ifdef HYDRO_BOOST_THREADS
		(*it)->setHelper(helperVector[0]);
#endif
		(*it)->addInitialFlow();
	}
	flipEdge();		
}

void CMesh::initNS() {	
	list<CCell*>::iterator it, end;
	for ( it=activeCells.begin(), end=activeCells.end(); it != end; ++it)
		(*it)->initNS();
	flipEdge();	
}

void CMesh::setTau(double mT) {
	list<CCell*>::iterator it, end;
	for ( it=activeCells.begin(), end=activeCells.end(); it != end; ++it)
		(*it)->setTau(mT);
}

double CMesh::getF(double beta, double dF/* =0.*/){
	double mF;
	if (mFBose){
		return 1./(exp(beta) - 1.) * (1. + dF);
	}
	else if (mFFerm)
		return 1./(exp(beta) + 1.) * (1. + dF);
	else 
		return exp(-beta) * (1. + dF);
}

double CMesh::getDULocal(int eta, int x, int y, int u, int v){
	if (v==1) 
		return (getS(eta,x+1,y,u) - getS(eta,x-1,y,u))/(getX(eta,x+1,y,v) - getX(eta,x-1,y,v));
	if (v==2) 
		return (getS(eta,x,y+1,u) - getS(eta,x,y-1,u))/(getX(eta,x,y+1,v) - getX(eta,x,y-1,v));
	if (v==3) 
		return getU0(eta,x,y)/getTau() 
		+ (getS(eta+1,x,y,u) - getS(eta-1,x,y,u))/(getX(eta+1,x,y,v) - getX(eta-1,x,y,v));
	
	return 0.;
}

bool CMesh::anyActiveCells() {
	list<CCell*>::iterator it, end;
	for ( it=activeCells.begin(), end=activeCells.end(); it != end; ++it) {
#ifdef HYDRO_BOOST_THREADS
		(*it)->setEos(eosVector[0]);
#endif
		if ((*it)->getTCalc() > mFOTemp)
			return true;
	}
	
	return false;
}

void CMesh::copyActive(CMesh* mMesh) {
	mMesh->setNSize(mNSize);
	mMesh->setXSize(mXSize);
	mMesh->setYSize(mYSize);
	
	list<CCell*>::iterator it1, end1, it2, end2;
	it1  = activeCells.begin();
	end1 = activeCells.end();
	
	it2  = mMesh->activeCells.begin();
	end2 = mMesh->activeCells.end();
	
	while (it1 != end1 && it2 != end2 ) {
		if (*(*it1) != *(*it2)){
			mMesh->activeCells.remove(*it2);
			(*it2)->deactivate();
			it2++;
		}
		else {
			it1++;
			it2++;
		}
	}
	
	while (it2 != end2){
		mMesh->activeCells.remove(*it2);
		(*it2)->deactivate();
		it2++;
	}
		
}

void CMesh::copyActive(CMesh* mM1, CMesh* mM2) {
	mM1->setNSize(mNSize);
	mM1->setXSize(mXSize);
	mM1->setYSize(mYSize);
	
	mM2->setNSize(mNSize);
	mM2->setXSize(mXSize);
	mM2->setYSize(mYSize);
	
	list<CCell*>::iterator it1, end1, it2, end2, it3, end3;
	it1  = activeCells.begin();
	end1 = activeCells.end();
	
	it2  = mM1->activeCells.begin();
	end2 = mM1->activeCells.end();
	
	it3  = mM2->activeCells.begin();
	end3 = mM2->activeCells.end();
	
	while (it1 != end1) {
		
		if (*(*it1) != *(*it2)){
			(*it2)->setActive(false);
			(*it3)->setActive(false);
			
#ifdef HYDRO_BOOST_THREADS
			mM1->removeCell(it2++);
			mM2->removeCell(it3++);
#else 
			mM1->activeCells.erase(it2++);
			mM2->activeCells.erase(it3++);
#endif
			
		}
		else {
			++it1;
			++it2;
			++it3;
		}
	}
	
	while (it2 != end2){
		(*it2)->setActive(false);
		mM1->activeCells.erase(it2++);
		
		(*it3)->setActive(false);
		mM2->activeCells.erase(it3++);
	}
	
		//	assert(activeCells.size() == mM1->activeCells.size());
		//	assert(activeCells.size() == mM2->activeCells.size());
}

void CMesh::cleanActiveCells() {
	list<CCell*>::iterator it = activeCells.begin(), end = activeCells.end();
	for (;it != end; ++it) 
		if ( !(*it)->getActive()) {
#ifdef HYDRO_BOOST_THREADS
			removeCell(it);
#else
			activeCells.erase(it);
#endif
		}
	
}

bool CMesh::detectCrash(){
	list<CCell*>::iterator it, end;
	for ( it=activeCells.begin(), end=activeCells.end(); it != end; ++it)
		for (int l=0;l<11;l++) 
			if ( (*it)->getS(l) != (*it)->getS(l)) {
				printf("\n\ncrash cell (t=%0.6g) : n=%g x=%g y=%g \n",
					   getTau(),(*it)->getX(3),(*it)->getX(1),(*it)->getX(2));
				(*it)->print();		
				return false;
			}
	
	return true;
}

	// calculates how this cell should change
	// and puts the answer into the provided mesh, same cell
void CMesh::forward(CMesh* mMesh, double mDt) {
	
		// setting dt anywhere changes it everywhere
	getCell(0,0,0)->setDt(mDt);
	sLoss = 0.; eLoss = 0.;
	
	if (mPrintMs) {
		std::cout << std::endl << "pre-half..." << std::endl;
		CCell *mC = getCell(mPrintN,mPrintX,mPrintY);
#ifdef HYDRO_BOOST_THREADS
		mC->setEos(eosVector[0]);
		mC->setHelper(helperVector[0]);
#endif
		mC->update();
		mC->print();
		mC->printM();
	}
	
	/*
	for (int i=0;i<NTHREADS-1;i++) {
		std::cout << std::endl << i << "  0 **********" << std::endl;
		(*marks[i][0])->selfPrint();
		std::cout << std::endl <<  i << "  1 ***********" << std::endl;
		(*marks[i][1])->selfPrint(); fflush(stdout);
		std::cout << std::endl;
	}
	(*marks[NTHREADS-1][0])->selfPrint();
	
	std::cout << activeCells.size() << " " << mMesh->activeCells.size() << std::endl; fflush(stdout);
	*/
	 
#ifndef HYDRO_BOOST_THREADS
	std::list<CCell*>::iterator it1, end1, it2, end2;
	for (       it1=activeCells.begin(),        end1=activeCells.end(),
		 it2=mMesh->activeCells.begin(), end2=mMesh->activeCells.end(); 
		 it1 != end1; ++it1, ++it2) {
		(*it1)->forward(*it2);
	}
#else 
	boost::thread_group g;
	for (int i=0;i<NTHREADS;i++)
			g.add_thread( new boost::thread(&CMesh::tForward, this, mMesh, i));
	g.join_all();
	
#endif
	mMesh->flipEdge();
}

	// calculates how this cell should change
	// and puts the answer into the provided mesh, same cell
void CMesh::forward(CMesh* onMesh, CMesh* offMesh, double mDt) {
		// setting dt anywhere changes it everywhere
	
	getCell(0,0,0)->setDt(mDt);
	sLoss = 0.; eLoss = 0.;
	
	if (mPrintMs) {
		std::cout << std::endl << "pre-whole ..." << std::endl;
		CCell *mC = offMesh->getCell(mPrintN,mPrintX,mPrintY);
#ifdef HYDRO_BOOST_THREADS
		mC->setEos(eosVector[0]);
		mC->setHelper(helperVector[0]);
#endif
		mC->update();
		mC->print();
		mC->printM();
	}
	
#ifndef HYDRO_BOOST_THREADS
	list<CCell*>::iterator it1, end1, it2, end2, it3, end3;
	for (         it1=activeCells.begin(),          end1=activeCells.end(),
		  it2=onMesh->activeCells.begin(),  end2=onMesh->activeCells.end(),
		 it3=offMesh->activeCells.begin(), end3=offMesh->activeCells.end(); 
		 it1 != end1; ++it1, ++it2, ++it3) 
		(*it1)->forward(*it2,*it3);
#else 
	boost::thread_group g;
	for (int i=0;i<NTHREADS;i++)
		g.add_thread( new boost::thread(&CMesh::tForward2, this, onMesh, offMesh, i));
	g.join_all();
#endif
	
	if (mPrintMs) {
		std::cout << "post-whole(t=" << onMesh->getTau() << ')' << std::endl;
		CCell *mC = onMesh->getCell(mPrintN,mPrintX,mPrintY);
#ifdef HYDRO_BOOST_THREADS
		mC->setEos(eosVector[0]);
		mC->setHelper(helperVector[0]);
#endif
		mC->update();
		mC->print();
		mC->printM();
	}
	
	onMesh->flipEdge();
}

void CMesh::forward(CMesh* k0, CMesh* k1, CMesh* k2, CMesh* k3, double mDt) {
	std::cout << "CMesh::forward(mRK4) just placeholder..." << std::endl;
	exit(1);
}

#ifdef HYDRO_BOOST_THREADS
void CMesh::tForward(CMesh* mMesh, int iThread){
	std::list<CCell*>::iterator it1;
	std::list<CCell*>::iterator it2;
	
	for ( it1 = marks[iThread][0], it2 = mMesh->marks[iThread][0];
		 it1 != marks[iThread][1]; ++it1, ++it2) {
		
		(*it1)->setEos(eosVector[iThread]);
		(*it1)->setHelper(helperVector[iThread]);
		(*it2)->setEos(eosVector[iThread]);
		(*it2)->setHelper(helperVector[iThread]);
		
		(*it1)->forward(*it2);
	}
}

void CMesh::tForward2(CMesh* mM1, CMesh* mM2, int iThread){
	std::list<CCell*>::iterator it1;
	std::list<CCell*>::iterator it2;
	std::list<CCell*>::iterator it3;
	
	for ( it1 = marks[iThread][0], it2 = mM1->marks[iThread][0], it3 = mM2->marks[iThread][0];
		 it1 != marks[iThread][1]; ++it1, ++it2, ++it3) {
		
		(*it1)->setEos(eosVector[iThread]);
		(*it1)->setHelper(helperVector[iThread]);
		(*it2)->setEos(eosVector[iThread]);
		(*it2)->setHelper(helperVector[iThread]);
		(*it3)->setEos(eosVector[iThread]);
		(*it3)->setHelper(helperVector[iThread]);
		
		(*it1)->forward(*it2,*it3);
	}
}

void CMesh::removeCell(list<CCell*>::iterator it) {
	bool ifMatch = false;
	for (int i=0;i<NTHREADS && !ifMatch;i++)
		if (it == marks[i][0]){
			ifMatch = true;
			list<CCell*>::iterator nMark = activeCells.erase(it);
			marks[i][0] = nMark;
			if (i>0)
				marks[i-1][1] = nMark;
		}
	
	if (it == marks[NTHREADS-1][1]){
		ifMatch = true;
		marks[NTHREADS-1][1] = activeCells.end();
	}
	
	if (!ifMatch)
		activeCells.erase(it);
}

void CMesh::setEos(CEos* mEos) {
	list<CCell*>::iterator it, end;
	for ( it=activeCells.begin(), end=activeCells.end(); it != end; ++it)
		(*it)->setEos(mEos);
}

void CMesh::setHelper(CCellHelper* mHelper) {
	list<CCell*>::iterator it, end;
	for ( it=activeCells.begin(), end=activeCells.end(); it != end; ++it)
		(*it)->setHelper(mHelper);
}
#endif 

	//quadratic extrapolation from three evenly spaced cells to their neighbor
	//used to smoothen grid edges
void CMesh::smoothEdge(CCell* c1,CCell* c2, CCell* c3, CCell* cOut){
	
		// quadratic fits to velocity
	for (int i=1;i<4;i++)	 
		if (c3->getS(i) == 0.) cOut->setS(i,0.);
		else cOut->setS(i, c1->getS(i) - 3.*(c2->getS(i) - c3->getS(i)));
		//  	  else cOut->setS(i, - c2->getS(i) + 2.*c3->getS(i));
	
		// quadratic fits to scaled momentum anisotropies (a_i/alpha)
	if (!ISRESCALE)
		for (int i=4;i<11;i++) 
			cOut->setS(i, c1->getS(i) - 3.*c2->getS(i) + 3.*c3->getS(i));
	else 
		for (int i=4;i<11;i++) {
			if (c2->getS(i) != 0.)
				cOut->setS(i, c1->getS(i) * pow(c3->getS(i)/c2->getS(i),3));
			else 
				cOut->setS(i,0.);
		}
		//	cOut->setS(i, - c2->getS(i) + 2.*c3->getS(i));
	
		//logarithmic smoothing for e
	
		// forces exponential tail
		//	cOut->setE( c3->getE() * (c3->getE()/c2->getE()));
	
		// allows exponential and gaussian components -- *******seems unstable to tail curling*******
	cOut->setE(c1->getE()*pow( (c3->getE()/c2->getE()) ,3));
	
}

void CMesh::calcIntegrals(){
	intS  = 0.; intE  = 0.;	intEp = 0.;	intEx = 0.;	intVr = 0.;
	
	double intEpD = 0., intExD = 0., intVrD = 0.;
	
	list<CCell*>::iterator it, end;
	for ( it=activeCells.begin(), end=activeCells.end(); it != end; ++it) {
		
		double fact=1.;
		for (int i=0;i<3;i++)
			for (int j=0;j<2;j++)
				if ( !(*it)->getNeighbor(i,j)->getActive()) 
					fact *= 0.5;
		
		(*it)->selfUpdate();
		
		intS   += fact * (*it)->getU0() * (*it)->getS();
		intE   += fact * (*it)->getTxy(0,0);
		intEp  += fact * ((*it)->getTxy(1,1) - (*it)->getTxy(2,2));
		intEpD += fact * ((*it)->getTxy(1,1) + (*it)->getTxy(2,2));
		intEx  += fact * ( pow((*it)->getX(2),2) - pow((*it)->getX(1),2)) * (*it)->getE();
		intExD += fact * ( pow((*it)->getX(2),2) + pow((*it)->getX(1),2)) * (*it)->getE();
		intVr  += fact * (*it)->getE() * (*it)->getGammaTrans() * (*it)->getVr();
		intVrD += fact * (*it)->getE() * (*it)->getGammaTrans();
	}
	
	if (mBjorken) {
		intS *= getTau()*mDx*mDy*mDn;
		intE *= getTau()*mDx*mDy*mDn;
	}
	else {
		intS *= mDx*mDy*mDn;
		intE *= mDx*mDy*mDn;
	}
	
	intEp /= intEpD;
	intEx /= intExD;
	intVr /= intVrD;
}

void CMesh::calcLossIntegrals(){
	
	double tempE[3]; tempE[0]=0.; tempE[1]=0.; tempE[2]=0.;
	double tempS[3]; tempS[0]=0.; tempS[1]=0.; tempS[2]=0.;
	
	double mTau = getTau();
	double temp = 0.;
	
	eLoss = 0.; sLoss = 0.;
	
	list<CCell*>::iterator it, end;
	for ( it=activeCells.begin(), end=activeCells.end(); it != end; ++it) {
		
		double fact=1.;
		for (int i=0;i<3;i++)
			for (int j=0;j<2;j++)
				if ( !(*it)->getNeighbor(i,j)->getActive()) 
					fact *= 0.5;
					
		temp += fact * (*it)->getTxyCalc(3,3);
		
		if (fact < 1.0) {
			for (int i=0;i<3;i++)
				for (int j=0;j<2;j++)
					if ( !(*it)->getNeighbor(i,j)->getActive()) {
						tempE[i] += fact * (*it)->getTxy(0,i);
						tempS[i] += fact * (*it)->getS() * (*it)->getS(i);
					}
		}
	}
	
	if (mBjorken)
		eLoss = mTau*mTau*mDx*mDy*mDn*temp;
	else 
		eLoss = mDx*mDy*mDn*temp;
	
	if (mBjorken){
		eLoss += mTau*(mDx*mDy*tempE[2] + mDx*mDn*tempE[1] + mDy*mDn*tempE[0]);
		sLoss += mTau*(mDx*mDy*tempS[2] + mDx*mDn*tempS[1] + mDy*mDn*tempS[0]);
	}
	else {
		eLoss += mDx*mDy*tempE[2] + mDx*mDn*tempE[1] + mDy*mDn*tempE[0];
		sLoss += mDx*mDy*tempS[2] + mDx*mDn*tempS[1] + mDy*mDn*tempS[0];
	}
	
}

double CMesh::getFOSX() {
	for (int i=mXSize; i>0; i--) {
		double mT = getT(0,i-1,0);
		if ( mT > mFOTemp) {
			double v1 =  getX(0,i-1,0,1);
			double v2 =  getX(0,i,0,1);
			return v1 + (v2-v1) * ((mFOTemp - mT) /(getT(0,i,0) - mT));
		}
	}
	return 0.;
}

double CMesh::getFOSY() {
	for (int i=mYSize; i>0; i--) {
		double mT = getT(0,0,i-1);
		if ( mT > mFOTemp) {
			double v1 =  getX(0,0,i-1,2);
			double v2 =  getX(0,0,i,2);
			return v1 + (v2-v1) * ((mFOTemp - mT) /(getT(0,0,i) - mT));
		}
	}
	return 0.;
}

double CMesh::getFOSSST() {
	for (int i=mXSize; i> 0; i--) {
		double mT = getT(0,i-1,0);
		if ( mT > mFOTemp) {
			double v1 =  ( getPixy(0,i-1,0,1,1) + getPixy(0,i-1,0,2,2))	/( getE(0,i-1,0) + getP(0,i-1,0));
			double v2 =  ( getPixy(0,i,0,1,1) + getPixy(0,i,0,2,2))	/( getE(0,i,0) + getP(0,i,0));
			return (1./JHBARC) * (v1 + (v2- v1) * ((mFOTemp - mT) /(getT(0,i,0) - mT)));
		}
	}
	return 0.;
}

double CMesh::getFOSV() {
	for (int i=mXSize; i>0; i--) {
		double mT = getT(0,i-1,0);
		if ( mT > mFOTemp) {
			double v1 =  getS(0,i-1,0,1);
			double v2 =  getS(0,i,0,1);
			return v1 + (v2-v1) * ((mFOTemp - mT) /(getT(0,i,0) - mT));
		}
	}
	return 0.;
}

bool CMesh::containsSurf(double &fosE, double localE []){
	int ieLength;
	
	if (mPureBjorken)
		ieLength = 8;
	else 
		ieLength = 16;
	
	if (localE[0] > fosE){
		for (int i=1;i<ieLength;i++)
			if (localE[i] < fosE) return true;
		return false;
	} 
	else 
		for (int i=1;i<ieLength;i++)
			if (localE[i] > fosE) return true;
	return false;
}

	//bool CMesh::containsSurf(double &fosE, double localE[2][2][2][2]){
bool CMesh::containsSurf(double &fosE, double**** localE){
	
	if (localE[0][0][0][0] > fosE){
		for (int i=0;i<2;i++)
			for (int j=0;j<2;j++)
				for (int k=0;k<2;k++)
					for (int l=0;l<2;l++)
						if (localE[i][j][k][l] < fosE) return true;
		return false;
	} 
	else 
		for (int i=0;i<2;i++)
			for (int j=0;j<2;j++)
				for (int k=0;k<2;k++)
					for (int l=0;l<2;l++)
						if (localE[i][j][k][l] > fosE) return true;
	return false;
}

void CMesh::cubeInterp(double mS[11], double midX[3], CCell* cube[8], double mDt) {
	double mT = midX[0]/mDt;
	double mX = midX[1]/mDx;
	double mY = midX[2]/mDy;
	
	for (int i=0;i<11;i++){
		mS[i]  = cube[0]->getS(i)*(1.-mT)*(1.-mX)*(1.-mY);
		mS[i] += cube[1]->getS(i)*(1.-mT)*mX*(1.-mY);
		mS[i] += cube[2]->getS(i)*(1.-mT)*mX*mY;
		mS[i] += cube[3]->getS(i)*(1.-mT)*(1.-mX)*mY;
		mS[i] += cube[4]->getS(i)*mT*(1.-mX)*(1.-mY);
		mS[i] += cube[5]->getS(i)*mT*mX*(1.-mY);
		mS[i] += cube[6]->getS(i)*mT*mX*mY;
		mS[i] += cube[7]->getS(i)*mT*(1.-mX)*mY;
	}
}

void CMesh::hcubeInterp(double mS[11], double foPoint[4], CCell* cube[2][2][2][2], double mDt){
	for (int i=0;i<11;i++) {
		mS[i] = 0.;
		for (int j=0;j<2;j++)
			for (int k=0;k<2;k++)
				for (int l=0;l<2;l++)
					for (int m=0;m<2;m++){
						double delta = (1.- fabs(foPoint[0]-cube[j][k][l][m]->getTau())/mDt);
						delta *= (1.- fabs(foPoint[1]-cube[j][k][l][m]->getX())/mDx);
						delta *= (1.- fabs(foPoint[2]-cube[j][k][l][m]->getY())/mDy);
						delta *= (1.- fabs(foPoint[3]-cube[j][k][l][m]->getEta())/mDn);
						
						mS[i] += cube[j][k][l][m]->getS(i)*delta;
					}
	}
}

int CMesh::cubical2p1(double e0, double eCubic[8], double surfvec[3], 
					  int* nsurfptr, double Vmid[3], double dt, double dx, 
					  double dy, int* Nerrptr){
    
		double ibit[12][12] = {{0,1,2,1,1,1,0,0,2,0,0,0},
			{1,0,1,2,0,1,1,0,0,2,0,0},
			{2,1,0,1,0,0,1,1,0,0,2,0},
			{1,2,1,0,1,0,0,1,0,0,0,2},
			{1,0,0,1,0,2,0,2,1,0,0,1},
			{1,1,0,0,2,0,2,0,1,1,0,0},
			{0,1,1,0,0,2,0,2,0,1,1,0},
			{0,0,1,1,2,0,2,0,0,0,1,1},
			{2,0,0,0,1,1,0,0,0,1,2,1},
			{0,2,0,0,0,1,1,0,1,0,1,2},
			{0,0,2,0,0,0,1,1,2,1,0,1},
			{0,0,0,2,1,0,0,1,1,2,1,0}
		};
		
		double vl0 = 0.0;        //to determine the positive direction of the normal vector of the surface
		double vl1 = 0.0;
		double vl2 = 0.0;
		double vh0 = 0.0;
		double vh1 = 0.0;
		double vh2 = 0.0;
		double vd0 = 0.0;
		double vd1 = 0.0;
		double vd2 = 0.0;
		double elsum = 0.0;
		double ehsum = 0.0;
		
		int nsurf = 0;   //store the number of the edge that is cut by the surface
		int Nerr = *Nerrptr;
		int iedge[12];           //the number of the edge
		double cuts[12][3];      //store the intersection along each edge fo the cubic
		int Isid[6];             //used for sorting the cuts
		double surf[12][3];
		double ad[3],bd[3];      //for calculating normal vector
		
			//   initialized the array
		for(int i=0;i<12;i++)
		{
			for(int j=0;j<3;j++)
			{
				cuts[i][j] = 0.0;
				surf[i][j] = 0.0;
			}
			iedge[i] = 0;
		}
		for(int i=0;i<6;i++)
			Isid[i] = 0;
		for(int i=0;i<3;i++)
		{
			ad[i]=bd[i]=0.0;
			surfvec[i]=0.0;
		}
		
			//   determine the intersection on each edge of the cubic
			//   there are totally 12 edges on the cublic
		
		double ek = eCubic[0];
		double el = eCubic[1];
		double dek = (e0-ek);
		double adek = fabs(dek);
		if(dek*(el-e0) >= 0.0)
		{
			iedge[nsurf] = 0;
			cuts[nsurf][0] = 0.0;
			cuts[nsurf][1] = (e0-ek)/(el-ek)*dx;
			cuts[nsurf][2] = 0.0;
			nsurf = nsurf + 1;
		}
		if(dek > 0.0)
			elsum += adek;
		else
			ehsum += adek;
		
		ek = eCubic[1];
		el = eCubic[2];
		dek = (e0-ek);
		adek = fabs(dek);
		if(dek*(el-e0) >= 0.0)
		{
			iedge[nsurf] = 1;
			cuts[nsurf][0] = 0.0;
			cuts[nsurf][1] = dx;
			cuts[nsurf][2] = (e0-ek)/(el-ek)*dy;
			nsurf = nsurf + 1;
		}
		if(dek > 0.0)
		{
			elsum += adek;
			vl1 += adek;
		}
		else
		{
			ehsum += adek;
			vh1 += adek; 
		}
		
		ek = eCubic[3];
		el = eCubic[2];
		dek = (el-e0);
		adek = fabs(dek);
		if(dek*(e0-ek) >= 0.0)
		{
			iedge[nsurf] = 2;
			cuts[nsurf][0] = 0.0;
			cuts[nsurf][1] = (e0-ek)/(el-ek)*dx;
			cuts[nsurf][2] = dy;
			nsurf = nsurf + 1;
		}
		if(dek < 0.0)
		{
			elsum += adek;
			vl1 += adek;
			vl2 += adek;
		}
		else
		{
			ehsum += adek;
			vh1 += adek; 
			vh2 += adek;
		}
		
		ek = eCubic[0];
		el = eCubic[3];
		dek = (el-e0);
		adek = fabs(dek);
		if(dek*(e0-ek) >= 0.0)
		{
			iedge[nsurf] = 3;
			cuts[nsurf][0] = 0.0;
			cuts[nsurf][1] = 0.0;
			cuts[nsurf][2] = (e0-ek)/(el-ek)*dy;
			nsurf = nsurf + 1;
		}
		if(dek < 0.0)
		{
			elsum += adek;
			vl2 += adek;
		}
		else
		{
			ehsum += adek;
			vh2 += adek;
		}
		
		ek = eCubic[0];
		el = eCubic[4];
		dek = (el-e0);
		adek = fabs(dek);
		if(dek*(e0-ek) >= 0.0)
		{
			iedge[nsurf] = 4;
			cuts[nsurf][0] = (e0-ek)/(el-ek)*dt;
			cuts[nsurf][1] = 0.0;
			cuts[nsurf][2] = 0.0;
			nsurf = nsurf + 1;
		}
		if(dek < 0.0)
		{
			elsum += adek;
			vl0 += adek;
		}
		else
		{
			ehsum += adek;
			vh0 += adek; 
		}
		
		ek = eCubic[1];
		el = eCubic[5];
		dek = (el-e0);
		adek = fabs(dek);
		if(dek*(e0-ek) >= 0.0)
		{
			iedge[nsurf] = 5;
			cuts[nsurf][0] = (e0-ek)/(el-ek)*dt;
			cuts[nsurf][1] = dx;
			cuts[nsurf][2] = 0.0;
			nsurf = nsurf + 1;
		}
		if(dek < 0.0)
		{
			elsum += adek;
			vl0 += adek;
			vl1 += adek;
		}
		else
		{
			ehsum += adek;
			vh0 += adek; 
			vh1 += adek;
		}
		
		ek = eCubic[2];
		el = eCubic[6];
		dek = (el-e0);
		adek = fabs(dek);
		if(dek*(e0-ek) >= 0.0)
		{
			iedge[nsurf] = 6;
			cuts[nsurf][0] = (e0-ek)/(el-ek)*dt;
			cuts[nsurf][1] = dx;
			cuts[nsurf][2] = dy;
			nsurf = nsurf + 1;
		}
		if(dek < 0.0)
		{
			elsum += adek;
			vl0 += adek;
			vl1 += adek;
			vl2 += adek;
		}
		else
		{
			ehsum += adek;
			vh0 += adek; 
			vh1 += adek;
			vh2 += adek;
		}
		
		ek = eCubic[3];
		el = eCubic[7];
		dek = (el-e0);
		adek = fabs(dek);
		if(dek*(e0-ek) >= 0.0)
		{
			iedge[nsurf] = 7;
			cuts[nsurf][0] = (e0-ek)/(el-ek)*dt;
			cuts[nsurf][1] = 0.0;
			cuts[nsurf][2] = dy;
			nsurf = nsurf + 1;
		}
		if(dek < 0.0)
		{
			elsum += adek;
			vl0 += adek;
			vl2 += adek;
		}
		else
		{
			ehsum += adek;
			vh0 += adek; 
			vh2 += adek;
		}
		
		ek = eCubic[4];
		el = eCubic[5];
		if((e0-ek)*(el-e0) >= 0.0)
		{
			iedge[nsurf] = 8;
			cuts[nsurf][0] = dt;
			cuts[nsurf][1] = (e0-ek)/(el-ek)*dx;
			cuts[nsurf][2] = 0.0;
			nsurf = nsurf + 1;
		}
		
		ek = eCubic[5];
		el = eCubic[6];
		if((e0-ek)*(el-e0) >= 0.0)
		{
			iedge[nsurf] = 9;
			cuts[nsurf][0] = dt;
			cuts[nsurf][1] = dx;
			cuts[nsurf][2] = (e0-ek)/(el-ek)*dy;
			nsurf = nsurf + 1;
		}
		
		ek = eCubic[7];
		el = eCubic[6];
		if((e0-ek)*(el-e0) >= 0.0)
		{
			iedge[nsurf] = 10;
			cuts[nsurf][0] = dt;
			cuts[nsurf][1] = (e0-ek)/(el-ek)*dx;
			cuts[nsurf][2] = dy;
			nsurf = nsurf + 1;
		}
		
		ek = eCubic[4];
		el = eCubic[7];
		if((e0-ek)*(el-e0) >= 0.0)
		{
			iedge[nsurf] = 11;
			cuts[nsurf][0] = dt;
			cuts[nsurf][1] = 0.0;
			cuts[nsurf][2] = (e0-ek)/(el-ek)*dy;
			nsurf = nsurf + 1;
		}
		
		/*****************************************************************
		 ** VDN:S CALCULATED BELOW ARE THE COMPONENTS OF THE DIFFERENCE **
		 ** OF WEIGHTED AVERAGES OF POINTS WITH E<E0 AND E>E0.          **
		 ** THIS VECTOR IS USED TO DEFINE THE POSITIVE SIDE OF SURFACE  **
		 ** ELEMENTS IN THIS CUBE TO BE TOWARDS LOWER ENERGY.           **
		 *****************************************************************/
		
		if(elsum != 0.0)
		{
			vl0 = vl0*dt/elsum;
			vl1 = vl1*dx/elsum;
			vl2 = vl2*dy/elsum;
		}
		
		if(ehsum != 0.0)
		{
			vh0 = vh0*dt/ehsum;
			vh1 = vh1*dx/ehsum;
			vh2 = vh2*dy/ehsum;
		}
		
		vd0 = vl0 - vh0;
		vd1 = vl1 - vh1;
		vd2 = vl2 - vh2;
		
		/*********************************************************
		 *** SORT THE INTERSECTION POINTS INTO A CIRCULAR ORDER **
		 *** USING A BITCHART IN MATRIX 'IBIT'.                 **
		 *********************************************************/
		
		for(int i=0; i<(nsurf-2); i++)
		{
			int Ie = iedge[i];
			int Is = 0;
			for(int j=i+1;j<nsurf;j++)
			{
				int Je = iedge[j];
				if(ibit[Ie][Je] != 0)
				{
					Isid[Is]=j;
					Is++;
				}
			}
			int Jmin = 0;
			if(Is != 1)
			{
				if(Is == 0)
				{
					Nerr++;
					cout<<"error" << endl;
					surfvec[0] = 0.0;
					surfvec[1] = 0.0;
					surfvec[2] = 0.0;
					*Nerrptr = Nerr;
					*nsurfptr = nsurf;
					return(0);
				}
				int minpts = 100000;
				for(int j=0;j<Is;j++)
				{
					int Je = iedge[Isid[j]];
					int Ipts = 0;
					for(int k=0;k<i-1;k++)
					{
						int Ke = iedge[k];
						if(ibit[Ke][Je] > 1)
							Ipts += 10;
						else
							Ipts += 10*ibit[Ke][Je];
					}
					Ipts = Ipts + ibit[Je][Ie] -1;
					if(Ipts < minpts)
					{
						Jmin = Isid[j];
						minpts = Ipts;
					}
				}
			}
			else
			{
				Jmin = Isid[0];
			}
			for(int k=0;k<3;k++)
			{   
				double apu = cuts[Jmin][k];
				cuts[Jmin][k] = cuts[i+1][k];
				cuts[i+1][k] = apu;
			}
			Ie = iedge[Jmin];
			iedge[Jmin] = iedge[i+1];
			iedge[i+1] = Ie;
		}
		
		int Nse = iedge[nsurf-1];
		int Nsm = iedge[nsurf-2];
		int I1 = iedge[0];
		
		if(ibit[Nse][Nsm] == 0 || ibit[Nse][I1] ==0)
		{
			Nerr++;
			cout<< "error!" << endl;
			surfvec[0] = 0.0;
			surfvec[1] = 0.0;
			surfvec[2] = 0.0;
			*nsurfptr = nsurf;
			*Nerrptr = Nerr;
			return(0);
		}
			//end sort
		
			//calculate the normal vector and midpoint of the intersection
		
		double v0 = 0.0;
		double v1 = 0.0;
		double v2 = 0.0;
		
		for(int n=0; n<nsurf; n++)
		{
			v0 = v0 + cuts[n][0];
			v1 = v1 + cuts[n][1];
			v2 = v2 + cuts[n][2];
		}
		
		Vmid[0] = v0/nsurf;
		Vmid[1] = v1/nsurf;
		Vmid[2] = v2/nsurf;
		
		for(int n=0; n<nsurf; n++)
		{
			int m = (n+1)%nsurf;
			for(int i=0;i<3;i++)
			{
				ad[i] = cuts[m][i]-cuts[n][i];
				bd[i] = Vmid[i]-cuts[n][i];
			}
			
				//calculate the covariant components of the surface vector
			double su0l = 0.5*(ad[2]*bd[1]-ad[1]*bd[2]);
			double su1l = 0.5*(ad[0]*bd[2]-ad[2]*bd[0]);
			double su2l = 0.5*(ad[1]*bd[0]-ad[0]*bd[1]);
			
				//choose the positive direction to be towards lower energy
			double vsum = vd0*su0l + vd1*su1l + vd2*su2l;
			int sign = int(vsum/fabs(vsum));
			
			surf[n][0] = sign*su0l;
			surf[n][1] = sign*su1l;
			surf[n][2] = sign*su2l;
		}
		
		for(int j=0;j<3;j++)
			for(int i=0;i<nsurf;i++)
				surfvec[j] += surf[i][j];
		
		/*     for(int i=0;i<nsurf;i++)
		 {
		 for(int j=0;j<3;j++)
         cout<<surf[i][j] << "   ";
		 cout<<endl;
		 }
		 for(int i=0;i<3;i++)
		 {
		 cout<<Vmid[i] <<endl;
		 }
		 cout<< nsurf<<endl;*/
		*nsurfptr = nsurf;
		*Nerrptr = Nerr;
		
		return(0);
		
}

	//statics 
parameterMap* CMesh::pMap;
int CMesh::mOctant;
bool CMesh::mPureBjorken, CMesh::mBjorken, CMesh::mPrintMs;
bool CMesh::mSVTrimInit, CMesh::mSVTrim, CMesh::mSawtooth, CMesh::bGSat, CMesh::bBGK;
double CMesh::mFOTemp, CMesh::mDeadT, CMesh::mT0, CMesh::mDx, CMesh::mDy, CMesh::mDn;
int CMesh::mNSize, CMesh::mXSize, CMesh::mYSize;
int CMesh::mNSizeOrig, CMesh::mXSizeOrig, CMesh::mYSizeOrig;
int CMesh::mPrintN, CMesh::mPrintX, CMesh::mPrintY;
double CMesh::wnRho0, CMesh::wnRAu, CMesh::wnXi, CMesh::wnSigma;
double CMesh::wnA, CMesh::wnB, CMesh::wnKTau, CMesh::gSatN;
double CMesh::mE0, CMesh::mRn, CMesh::mRt, CMesh::mRa;
double CMesh::mRx, CMesh::mRy, CMesh::mRSig, CMesh::mWnBinRatio;
double CMesh::mInitFlow, CMesh::mInitNS, CMesh::mInitNSTxxP;
double CMesh::icNiceA, CMesh::icNiceB, CMesh::collRootS;
double CMesh::fracPPE, CMesh::sigmaAtt, CMesh::wnAttRatio;
int CMesh::activeX, CMesh::activeY, CMesh::activeN;
double CMesh::wnTx, CMesh::wnTy;
string CMesh::mDataRoot;
FILE *CMesh::osuFOS;
Cornelius CMesh::mCornelius;

double ****CMesh::localE4;
double ***CMesh::localE3;
double *CMesh::dxCorn;

#ifdef HYDRO_BOOST_THREADS
CEos* CMesh::eosVector[NTHREADS];
CCellHelper* CMesh::helperVector[NTHREADS];
#endif