/*  Cmesh.cpp
 *  isHydro3
 *  Created by Joshua Vredevoogd on 2/11/09.
 */

#include "octMesh.h"

octMesh::octMesh(parameterMap* pM) : CMesh(pM) {
	
	mCells = new CCell***[mNSize+3];	
	for (int i=0;i<mNSize+3;i++){
		mCells[i] = new CCell**[mXSize+3];
		for (int j=0;j<mXSize+3;j++){
			mCells[i][j] = new CCell*[mYSize+3];
		}
	}
	
	if (parameter::getB(*pMap,"HYDRO_KLN_INPUT",false)) {
		klnInput();
		return;
	}
	
		// make cells at eta = 0
	for (int j=-1;j<=mXSize+1;j++)
		for (int k=-1;k<=mYSize+1;k++) {
			mCells[0][j+1][k+1] = new CCell(mT0, 
											-mDn,
											mDx*(double)(j),
											mDy*(double)(k));
			
			mCells[1][j+1][k+1] = new CCell(mT0, 
											0.,
											mDx*(double)(j),
											mDy*(double)(k));
			initialCondition(mCells[1][j+1][k+1]);
		}
	
		// now do the rest of the cells, using eta=0 IC for eta!=0				
	for (int i=1;i<=mNSize+1;i++)
		for (int j=-1;j<=mXSize+1;j++)
			for (int k=-1;k<=mYSize+1;k++){
				mCells[i+1][j+1][k+1] = new CCell(mT0, 
												  mDn*(double)(i),
												  mDx*(double)(j),
												  mDy*(double)(k));
				initialCondition(mCells[i+1][j+1][k+1],mCells[1][j+1][k+1]);
			}
	
		// Connect cells
	for (int i=0;i<=mNSize;i++)
		for (int j=0;j<=mXSize;j++)
			for (int k=0;k<=mYSize;k++) {
				mCells[i+1][j+1][k+1]->setEtaNeighbors( mCells[i][j+1][k+1], 
													   mCells[i+2][j+1][k+1]);
				mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
													 mCells[i+1][j+2][k+1]);
				mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
													 mCells[i+1][j+1][k+2]);
				
			}
	
		// Create a layer of dead cells around the edges
	for (int i=-1;i<=mNSize+1;i++) 
		for (int j=-1;j<=mXSize+1;j++) {
			mCells[i+1][j+1][0]->setActive(false);
			mCells[i+1][j+1][mYSize+2]->setActive(false);
		}
	
	for (int i=-1;i<=mNSize+1;i++)
		for (int j=-1;j<=mYSize+1;j++) {
			mCells[i+1][0][j+1]->setActive(false);
			mCells[i+1][mXSize+2][j+1]->setActive(false);
		}
	
	for (int i=-1;i<=mXSize+1;i++)
		for (int j=-1;j<=mYSize+1;j++) {
			mCells[0][i+1][j+1]->setActive(false);
			mCells[mNSize+2][i+1][j+1]->setActive(false);
		}
	
	flipEdge();
}

octMesh::octMesh(parameterMap* pM, const char *fName) {
	mCells = new CCell***[mNSize+2];	
	for (int i=0;i<mNSize+2;i++){
		mCells[i] = new CCell**[mXSize+2];
		for (int j=0;j<mXSize+2;j++){
			mCells[i][j] = new CCell*[mYSize+2];
		}
	}
	
	pMap = pM;
	CCell* dummyCell = new CCell(pMap);
	
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
	
	mWnBinRatio = parameter::getD(*pMap,"GLAUBER_WNBIN_RATIO",1.0);
	mInitFlow = parameter::getD(*pMap,"HYDRO_INIT_FLOW",0.0);
	
	mInitNS = parameter::getD(*pMap,"HYDRO_INIT_NS",0.0);
	mInitNSTxxP = parameter::getD(*pMap,"HYDRO_INIT_NS_TXXP",0.0);
	
	if (mInitNS != 0.0 && mInitNSTxxP != 0.0) {
		std::cout << "\nCMesh::CMesh(parameterMap* pM)....\nConflicting instructions for shear tensor initialization...\n";
		std::cout << "Use either HYDRO_INIT_NS or HYDRO_INIT_NS_TXXP but not both!\n  Aborting.\n\n";
		fflush(stdout);
		exit(1);
	}
	
	icNiceA = parameter::getD(*pMap, "IC_NICE_A", 1.4/ROOT2);
	icNiceB = parameter::getD(*pMap, "IC_NICE_B", 2.);
	
	
		//CMesh::CMesh(const string fName) {
	double mX[4];
	double mS[11];
	
	int active;
	
	std::ifstream mFile(fName);
	
		// make cells
	mFile >> mX[0] >> mNSize >> mXSize >> mYSize;
	
	if (!mPureBjorken) {
		for (int i=-1;i<=mNSize;i++)
			for (int j=-1;j<=mXSize;j++) 
				for (int k=-1;k<=mYSize;k++) {
					
					mFile >> active;
					
					for (int l=1;l<4;l++) 
						mFile >> mX[l];
					for (int l=0;l<11;l++) 
						mFile >> mS[l];
					
					mFile.ignore(256,'\n');
					
					mCells[i+1][j+1][k+1] = new CCell(mX[0],mX[3],mX[1],mX[2]);
						//printf("%d %d %d : %g %g %g\n",i,j,k,mX[1],mX[2],mX[3]);
					
					for (int l=0;l<11;l++) 
						mCells[i+1][j+1][k+1]->setS(l,mS[l]);
					
					if (active == 1)
						setActive(i,j,k,true);
					else 
						setActive(i,j,k,false);
				}
	}
	else 
		for (int j=-1;j<=mXSizeOrig;j++)
			for (int k=-1;k<=mYSizeOrig;k++) {
				for (int l=1;l<4;l++) 
					mFile >> mX[l];
				for (int l=0;l<11;l++)
					mFile >> mS[l];
				for (int i=-1;i<=mNSizeOrig;i++) {
					mCells[i+1][j+1][k+1] = new CCell(mX[0],mX[1],mX[2],mX[3]);
					for (int l=1;l<11;l++) mCells[i+1][j+1][k+1]->setS(l,mS[l]);
				}
			}
	
		// connect your cells...
	if (!mPureBjorken)
		for (int i=0;i<mNSize;i++)
			for (int j=0;j<mXSize;j++)
				for (int k=0;k<mYSize;k++) {
					mCells[i+1][j+1][k+1]->setEtaNeighbors( mCells[i][j+1][k+1], 
														   mCells[i+2][j+1][k+1]);
					mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
														 mCells[i+1][j+2][k+1]);
					mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
														 mCells[i+1][j+1][k+2]);
				}
	else 
		for (int i=0;i<mNSizeOrig;i++)
			for (int j=0;j<mXSizeOrig;j++)
				for (int k=0;k<mYSizeOrig;k++) {
					mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
														 mCells[i+1][j+2][k+1]);
					mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
														 mCells[i+1][j+1][k+2]);
				}
	
	delete dummyCell;
}

octMesh::octMesh(CMesh* mesh) {
	
	mCells = new CCell***[mNSizeOrig+3];	
	for (int i=0;i<mNSizeOrig+3;i++){
		mCells[i] = new CCell**[mXSizeOrig+3];
		for (int j=0;j<mXSizeOrig+3;j++){
			mCells[i][j] = new CCell*[mYSizeOrig+3];
		}
	}
	
	for (int i=-1;i<=mNSizeOrig+1;i++)
		for (int j=-1;j<=mXSizeOrig+1;j++)
			for (int k=-1;k<=mYSizeOrig+1;k++)
				mCells[i+1][j+1][k+1] = new CCell(mesh->mCells[i+1][j+1][k+1]);
	
		// connect your cells and (antiquated..) reduce boundary problem...
	for (int i=0;i<=mNSizeOrig;i++)
		for (int j=0;j<=mXSizeOrig;j++)
			for (int k=0;k<=mYSizeOrig;k++) {
				mCells[i+1][j+1][k+1]->setEtaNeighbors( mCells[i][j+1][k+1], 
													   mCells[i+2][j+1][k+1]);
				mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
													 mCells[i+1][j+2][k+1]);
				mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
													 mCells[i+1][j+1][k+2]);
			}
}

	// destructor
octMesh::~octMesh() {
	for (int i=0;i<mNSizeOrig+3;i++){
		for (int j=0;j<mXSizeOrig+3;j++){
			for (int k=0;k<mYSizeOrig+3;k++){
				delete mCells[i][j][k];
			}
			delete [] mCells[i][j];
		}
		delete [] mCells[i];
	}
	delete [] mCells;
}
					
void octMesh::klnInput() {
	
	string fName = parameter::getS(*pMap,"HYDRO_KLN_FILENAME"," ");
	FILE *profile = fopen(fName.c_str(),"r");
	
	char dummy;
	double density, dummy1, dummy2, dummy3, eta0, x0, y0;
	double deta=0.25, dx=0.08311, dy=0.08311;
	int neta=0, nx=60, ny=60;
	
		// loop over the whole lattice and initialize values:
		//    int bytes_read=fscanf(profile,"%s %s %s %s %d %s %d %s %d %s %lf %s %lf %s %lf",
		//					  &dummy,&dummy,&dummy,&dummy,&neta,&dummy,&nx,&dummy,&ny,&dummy,
		//					  &deta,&dummy,&dx,&dummy,&dy);
	int bytes_read=0;
	
	
	for (int i=0;i<=mXSize;i++)
		for (int j=0;j<=mYSize;j++) 
			for (int k=0;k<=mNSize;k++)
				mCells[k+1][i+1][j+1] = NULL;
	
	for (int i=0;i<=nx;i++)
		for (int j=0;j<=ny;j++) 
			for (int k=0;k<=neta;k++){
				
				bytes_read=fscanf(profile,"%lf %lf %lf %lf",&dummy1,&dummy2,&density,&dummy3);
				
				if (dummy1 < -.5*dx || dummy2 < -.5*dy){
					k--;
					continue;
				}
				
				dummy3 = k*deta;
				mCells[k+1][i+1][j+1] = new CCell(mT0, dummy3, dummy2, dummy1);
				
				if (i==1) mDx = mCells[k+1][i+1][j+1]->getX(1);
				if (j==1) mDy = mCells[k+1][i+1][j+1]->getX(2);
				if (k==1) mDn = mCells[k+1][i+1][j+1]->getX(3);
				
					// fill e with entropy density
				if (density != 0.) 
						//mCells[k+1][i+1][j+1]->setE(6.00689*density);
					mCells[k+1][i+1][j+1]->setE(density);
				else 
					mCells[k+1][i+1][j+1]->setE(1E-5);
			}
	
	for (int i=0;i<=mXSize;i++)
		for (int j=0;j<=mYSize;j++) 
			for (int k=0;k<=mNSize;k++)
				if (mCells[k+1][i+1][j+1] == NULL){
					mCells[k+1][i+1][j+1] = new CCell(mT0,k*mDn,i*mDx,j*mDy);
					mCells[k+1][i+1][j+1]->setE(1E-5);
				}
	
	fclose(profile);
	
	for (int i=-1;i<=mNSize+1;i++)
		for (int j=-1;j<=mXSize+1;j++) {
			mCells[i+1][j+1][0] = new CCell(mT0, mDn*i, mDx*j, -mDy);
			mCells[i+1][j+1][mYSize+2] = new CCell(mT0, mDn*i, mDx*j, mDy*(mYSize+1));
		}
	
	for (int i=-1;i<=mNSize+1;i++)
		for (int j=-1;j<=mYSize+1;j++) {
			mCells[i+1][0][j+1] = new CCell(mT0, mDn*i, -mDx, mDx*j);
			mCells[i+1][mXSize+2][j+1] = new CCell(mT0, mDn*i, mDx*(mXSize+1), mDy*j);
		}
	
	for (int i=-1;i<=mXSize+1;i++)
		for (int j=-1;j<=mYSize+1;j++) {
			mCells[0][i+1][j+1] = new CCell(mT0, -mDn, mDx*i, mDy*j);
			mCells[mNSize+2][i+1][j+1] = new CCell(mT0, mDn*(mNSize+1), mDx*i, mDy*j);
		}
	
		// Connect cells
	for (int i=0;i<=mNSize;i++)
		for (int j=0;j<=mXSize;j++)
			for (int k=0;k<=mYSize;k++) {
				mCells[i+1][j+1][k+1]->setEtaNeighbors( mCells[i][j+1][k+1], 
													   mCells[i+2][j+1][k+1]);
				mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
													 mCells[i+1][j+2][k+1]);
				mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
													 mCells[i+1][j+1][k+2]);
				
				if ( mCells[i+1][j+1][k+1]->getE() == 1E-5)
					mCells[i+1][j+1][k+1]->setActive(false);
			}
	
		// Create a layer of dead cells around the edges
	for (int i=-1;i<=mNSize+1;i++) 
		for (int j=-1;j<=mXSize+1;j++) {
			mCells[i+1][j+1][0]->setActive(false);
			mCells[i+1][j+1][mYSize+2]->setActive(false);
		}
	
	for (int i=-1;i<=mNSize+1;i++)
		for (int j=-1;j<=mYSize+1;j++) {
			mCells[i+1][0][j+1]->setActive(false);
			mCells[i+1][mXSize+2][j+1]->setActive(false);
		}
	
	for (int i=-1;i<=mXSize+1;i++)
		for (int j=-1;j<=mYSize+1;j++) {
			mCells[0][i+1][j+1]->setActive(false);
			mCells[mNSize+2][i+1][j+1]->setActive(false);
		}
	
	mCells[0][0][0]->setDx(1,mDx);
	mCells[0][0][0]->setDx(2,mDy);
	mCells[0][0][0]->setDx(3,mDn);
	flipEdge();
}

void octMesh::fillActiveCells() {
	if (mDeadT == 0.) {
		if (!mPureBjorken)
			for (int i=0;i<=mNSize;i++)
				for (int j=0;j<=mXSize;j++)
					for (int k=0;k<=mYSize;k++)
						activeCells.push_back(getCell(i,j,k));
		else 
			for (int j=0;j<=mXSize;j++)
				for (int k=0;k<=mYSize;k++)
					activeCells.push_back(getCell(0,j,k));
	}
	
	if (!mPureBjorken){
		for (int i=mNSize; i>=0; i--) 
			for (int j=mXSize; j>=0; j--) 
				for (int k=mYSize; k>=0; k--) {
#ifdef HYDRO_BOOST_THREADS
					setEos(i,j,k,eosVector[0]);
#endif
					if (!getActive(i,j,k)) 
						continue;
					else if (getT(i,j,k) > mDeadT) 
						activeCells.push_back(getCell(i,j,k));
					else {
						CCell *mc = deactivate(i,j,k);
						if (mc != NULL)
							activeCells.remove(mc);
					}
				}
	}
	else {
		for (int j=mXSize; j>=0; j--) 
			for (int k=mYSize; k>=0; k--) {
#ifdef HYDRO_BOOST_THREADS
				setEos(0,j,k,eosVector[0]);
#endif
				if (!getActive(0,j,k)) 
					continue;
				else if (getT(0,j,k) > mDeadT) 
					activeCells.push_back(getCell(0,j,k));
				else {
					CCell *mc = deactivate(0,j,k);
					if (mc != NULL)
						activeCells.remove(mc);
				}
			}
	}			

	for (int i=mNSize; i>0; i--) 
		if (!getActive(i,0,0)) mNSize = i;
		else break;
	
	for (int i=mXSize; i>0; i--) 
		if (!getActive(0,i,0)) mXSize = i;
		else break;
	
	for (int i=mYSize; i>0; i--) 
		if (!getActive(0,0,i)) mYSize = i;
		else break;
	
#ifdef HYDRO_BOOST_THREADS
	marks[0][0] = activeCells.begin();
	marks[NTHREADS-1][1] = activeCells.end();
	
	std::list<CCell*>::iterator it, end;
	int cells_per_thread = int(activeCells.size())/NTHREADS;
	int i=0;
	for (it=activeCells.begin(), end=activeCells.end();
		 it != end; i++, ++it){
		if (i%cells_per_thread == 0)  {
			int j = i/int(cells_per_thread);
			marks[j][0] = it;
			marks[j-1][1] = it;
			if (j==NTHREADS-1) break;
		}
			
	}
#endif
}

void octMesh::average(CMesh *m1, CMesh *m2) {
		// protection for averaging with yourself and another mesh (since getTau takes the central cell time, which will be changed)
	double mTau = getTau();
	double mTau1 = m1->getTau();
	double mTau2 = m2->getTau();
	
	for (int i=-1;i<=mNSize;i++)
		for (int j=-1;j<=mXSize;j++)
			for (int k=-1;k<=mYSize;k++)
				for (int s=1;s<11;s++)
					mCells[i+1][j+1][k+1]->setS(s, m1->getS(i,j,k,s) + ((mTau - mTau1)/(mTau2 - mTau1))
												*(m2->getS(i,j,k,s) - m1->getS(i,j,k,s)));	
}

void octMesh::copyMesh(CMesh *m) {
	for (int j=0;j<=mXSize;j++)
		for (int k=0;k<=mYSize;k++){
			if (!mPureBjorken){
				for (int i=0;i<=mNSize;i++)
					mCells[i+1][j+1][k+1]->copy(m->getCell(i,j,k));
			}
			else 
				mCells[1][j+1][k+1]->copy(m->getCell(0,j,k));
		}
}

	// converts a mesh in entropy density to one in energy density
void octMesh::fixEntropyMesh() {
	for (int i=-1;i<=mNSize+1;i++)
		for (int j=-1;j<=mXSize+1;j++)
			for (int k=-1;k<=mYSize+1;k++) {
				double mE = mCells[i+1][j+1][k+1]->getE();
				if (mE != 1E-5)
					mCells[i+1][j+1][k+1]->setE(wnKTau/getTau() * getEos()->getEGivenS(mE));
				else
					mCells[i+1][j+1][k+1]->setE(1E-10);
			}
}

void octMesh::update(int eta, int x, int y) {
	activeN=eta+1;
	activeX=x+1;
	activeY=y+1; 
	return getCell(eta,x,y)->update();
}

void octMesh::selfUpdate(int eta, int x, int y) {
	activeN=eta+1;
	activeX=x+1;
	activeY=y+1; 
	return getCell(eta,x,y)->selfUpdate();
}

void octMesh::deaden() {	
	
	if (mDeadT == 0.) return;
	bool flagCheckActive = false;  
	
	double tempS=0., tempE=0., tempEp=0., tempEx=0., tempVr=0.;
	if (!mPureBjorken)
		for (int i=mNSize; i>=0; i--) 
			for (int j=mXSize; j>=0; j--) 
				for (int k=mYSize; k>=0; k--) {
#ifdef HYDRO_BOOST_THREADS
					getCell(i,j,k)->setEos(eosVector[0]);
#endif
					if (!getActive(i,j,k)) 
						continue;
					if (getT(i,j,k) > mDeadT) 
						break;
					else {						
#ifdef HYDRO_BOOST_THREADS						
						for (list<CCell*>::iterator it  = activeCells.begin(), end = activeCells.end();
							 it != end; ++it) {
							if (*it == getCell(i,j,k)){
								removeCell(it);
								break;
							}
						}
#else	
						activeCells.remove(getCell(i,j,k));
#endif
									 
						CCell* mc = deactivate(i,j,k);
						if (mc == getCell(i,j,k))
							flagCheckActive = true;
						else if (mc != NULL) {
#ifdef HYDRO_BOOST_THREADS						
							for (list<CCell*>::iterator it  = activeCells.begin(), end = activeCells.end();
								 it != end; ++it) {
								if (*it == mc){
									removeCell(it);
									break;
								}
							}
#else	
							activeCells.remove(mc);
#endif
						}

						tempS  += getS(i,j,k)*getS(i,j,k,0);					
						tempE  += getTxy(i,j,k,0,0);
						tempEp += getTxy(i,j,k,2,2) - getTxy(i,j,k,1,1);
						tempEx -= (mDx*mDx*j*j - mDy*mDy*k*k)*getE(i,j,k);
						tempVr += getE(i,j,k) * getGammaTrans(i,j,k) * getVr(i,j,k);
					}
				}
	else {
		for (int j=mXSize; j>=0; j--) 
			for (int k=mYSize; k>=0; k--) {
#ifdef HYDRO_BOOST_THREADS
				getCell(0,j,k)->setEos(eosVector[0]);
#endif
				if (!getActive(0,j,k)) 
					continue;
				if (getT(0,j,k) > mDeadT) 
					break;
				else {
#ifdef HYDRO_BOOST_THREADS						
					for (list<CCell*>::iterator it  = activeCells.begin(), end = activeCells.end();
						 it != end; ++it) {
						if (*it == getCell(0,j,k)){
							removeCell(it);
							break;
						}
					}
#else	
					activeCells.remove(getCell(0,j,k));
#endif
					
					CCell* mc = deactivate(0,j,k);
					if (mc == getCell(0,j,k))
						flagCheckActive = true;
					else if (mc != NULL) {
#ifdef HYDRO_BOOST_THREADS						
						for (list<CCell*>::iterator it  = activeCells.begin(), end = activeCells.end();
							 it != end; ++it) {
							if (*it == mc){
								removeCell(it);
								break;
							}
						}
#else	
						activeCells.remove(mc);
#endif
					}
					
					tempS  += getS(0,j,k)*getS(0,j,k,0);					
					tempE  += getTxy(0,j,k,0,0);
					tempEp += getTxy(0,j,k,2,2) - getTxy(0,j,k,1,1);
					tempEx -= (mDx*mDx*j*j - mDy*mDy*k*k)*getE(0,j,k);
					tempVr += getE(0,j,k) * getGammaTrans(0,j,k) * getVr(0,j,k);
				}
			}
	}			
	
	if (flagCheckActive)
		cleanActiveCells();
							
	if (mBjorken){
		eDead  -= mDn*mDx*mDy * getTau() * tempE;
		sDead  -= mDn*mDx*mDy * getTau() * tempS;
		epDead -= mDn*mDx*mDy * getTau() * tempEp;
		exDead -= mDn*mDx*mDy * getTau() * tempEx;
		vrDead -= mDn*mDx*mDy * getTau() * tempVr;
	} 
	else {
		eDead  -= mDn*mDx*mDy * tempE;
		sDead  -= mDn*mDx*mDy * tempS;
		epDead -= mDn*mDx*mDy * tempEp;
		exDead -= mDn*mDx*mDy * tempEx;
		vrDead -= mDn*mDx*mDy * tempVr;
	}
	
	for (int i=mNSize; i>0; i--) 
		if (!getActive(i,0,0)) mNSize = i;
		else break;
	
	for (int i=mXSize; i>0; i--) 
		if (!getActive(0,i,0)) mXSize = i;
		else break;
	
	for (int i=mYSize; i>0; i--) 
		if (!getActive(0,0,i)) mYSize = i;
		else break;
	
}

void octMesh::flipEdge() {
	double mTau = getTau();
	
	if (!mPureBjorken) {
		for (int i=0;i<=mXSize;i++)
			for (int j=0;j<=mYSize;j++) {	
				flipEdge(3,mCells[2][i+1][j+1],
						 mCells[1][i+1][j+1],
						 mCells[0][i+1][j+1]);
				mCells[0][i+1][j+1]->setTau(mTau);
			}
			//inner rad bound
		for (int i=-1;i<=mNSize;i++) {
			for (int j=0;j<=mYSize;j++) {
				flipEdge(1,mCells[i+1][2][j+1],
						 mCells[i+1][1][j+1],
						 mCells[i+1][0][j+1]);
				mCells[i+1][0][j+1]->setTau(mTau);
			}
			for (int j=-1;j<=mXSize;j++) {	
				flipEdge(2,mCells[i+1][j+1][2],
						 mCells[i+1][j+1][1],
						 mCells[i+1][j+1][0]);
				mCells[i+1][j+1][0]->setTau(mTau);
			}
		}
	}
	else {
		for (int j=0;j<=mYSize;j++) {
			flipEdge(1,mCells[1][2][j+1],
					 mCells[1][1][j+1],
					 mCells[1][0][j+1]);
			mCells[1][0][j+1]->setTau(mTau);
		}
		for (int j=-1;j<=mXSize;j++) {	
			flipEdge(2,mCells[1][j+1][2],
					 mCells[1][j+1][1],
					 mCells[1][j+1][0]);
			mCells[1][j+1][0]->setTau(mTau);
		}	
	}
}

void octMesh::flipEdge(int fE, CCell* c1, CCell* c2, CCell* cOut) {
	cOut->copy(c1);
	cOut->flip(fE);
}

	// freezeout surface at eta=0
void octMesh::getFOS(double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], int &foSize) {
	foSize = 0;
	int j = mYSize;  //keep track of where we are in j
	
	for (int i=0; i<=mXSize; i++) {
		double v=0.;
		if (j+3 <= mYSize) j+=3;
		else j=mYSize;
		
		while (j>0) {
			if (getT(0,i,j-1) > mFOTemp) {
				mFOS[foSize][1] = getX(0,i,j-1,2) + (getX(0,i,j,2)-getX(0,i,j-1,2))
				*((mFOTemp-getT(0,i,j-1))/(getT(0,i,j)-getT(0,i,j-1)));
				mFOS[foSize][0] = getX(0,i,j-1,1);
				
				if (foSize==0) {
					mFOSigma[0][0] = 0.;
					mFOSigma[0][1] = 1.;
				}
				else if (foSize>1) {
					v = (mFOS[foSize-2][1] - mFOS[foSize][1])/(mFOS[foSize][0] - mFOS[foSize-2][0]);
					mFOSigma[foSize-1][0] = v/sqrt(1.+v*v);
					mFOSigma[foSize-1][1] = 1./sqrt(1.+v*v);
				}
				foSize++;
				break;
			}
			else j--;
		}
		if (v>1.) break;
		
		if (j==0 && i==0){
			foSize=0;
			return;
		}
	}
	
	j--;
	
	if (j==-2)
		j = floor( (int)(mFOS[foSize-1][1]/mDy));

	while(j>=0) {
		for (int i=mXSize;i>0;i--) {
			if (getT(0,i-1,j) > mFOTemp){
				mFOS[foSize][0] = getX(0,i-1,j,1) + (getX(0,i,j,1)-getX(0,i-1,j,1))*((mFOTemp-getT(0,i-1,j))/(getT(0,i,j)-getT(0,i-1,j)));
				mFOS[foSize][1] = getX(0,i-1,j,2);

				if (j==0) {
					mFOSigma[foSize][0] = 1.;
					mFOSigma[foSize][1] = 0.;
				}
				
				if (foSize > 1) {
					double v = (mFOS[foSize-2][0] - mFOS[foSize][0])/(mFOS[foSize][1] - mFOS[foSize-2][1]);
					mFOSigma[foSize-1][1] = v/sqrt(1.+v*v);
					mFOSigma[foSize-1][0] = 1./sqrt(1.+v*v);
				}
				
				foSize++;
				break;
			}
		}
		j--;
	}
		//	goodFOS(mFOS,foSize);
}

	// freezeout surface at eta=0
bool octMesh::getFOS(double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], 
					 double mFOVelo[XSIZE+YSIZE][2], double mFODiss[XSIZE+YSIZE][6], int &foSize) {

	foSize = 0;
	int j = mYSize;  //keep track of where we are in j
	double v=0.;
	
	for (int i=0; i<=mXSize; i++) {
		if (j+3 <= mYSize) j+=3;
		else j=mYSize;
		
		while (j>0) {
			if (getT(0,i,j-1) > mFOTemp) {

				double mTR = (mFOTemp-getT(0,i,j-1))/(getT(0,i,j)-getT(0,i,j-1));

				mFOS[foSize][1] = getX(0,i,j-1,2) + (getX(0,i,j,2) - getX(0,i,j-1,2))*mTR;
				mFOS[foSize][0] = getX(0,i,j-1,1);

				CCell foCell(getTau(),mFOS[foSize][0],mFOS[foSize][1],0.);

				for (int k=0;k<11;k++)
					foCell.setS(k,getS(0,i,j-1,k) + (getS(0,i,j,k)-getS(0,i,j-1,k))*mTR);

#ifdef HYDRO_BOOST_THREADS
				foCell.setEos(eosVector[0]);
				foCell.setHelper(helperVector[0]);
#endif
				foCell.selfUpdate();

				for (int k=1;k<3;k++) 
					mFOVelo[foSize][k-1] = foCell.getS(k);

				mFODiss[foSize][0] = 0.5*(foCell.getPixyLocal(1,1)-foCell.getPixyLocal(2,2));
				mFODiss[foSize][1] = 0.5*(foCell.getPixyLocal(1,1)+foCell.getPixyLocal(2,2)
										  -2*pow(foCell.getTau(),2)*foCell.getPixyLocal(3,3))/ROOT3;
				mFODiss[foSize][2] = foCell.getPixyLocal(1,2);
				mFODiss[foSize][3] = foCell.getTau()*foCell.getPixyLocal(1,3);
				mFODiss[foSize][4] = foCell.getTau()*foCell.getPixyLocal(2,3);
				mFODiss[foSize][5] = foCell.getB();

				if (foSize==0) {
					v=0.;
					mFOSigma[foSize][0] = 0.;
					mFOSigma[foSize][1] = 1.;
				}
				else if (foSize>1) {
					v = (mFOS[foSize-2][1] - mFOS[foSize][1])/(mFOS[foSize][0] - mFOS[foSize-2][0]);
					mFOSigma[foSize-1][0] = v/sqrt(1.+v*v);
					mFOSigma[foSize-1][1] = 1./sqrt(1.+v*v);
				}
				foSize++;
				break;
			}
			else j--;
		}
		if (v>1.) 
			break;
		
		if (j==0 && i==0){
			foSize=0;
			return false;
		}
	}
	
	j--;
	
	if (j==-2)
		j = (int)(mFOS[foSize-1][1]/mDy);
	
	while(j>=0) {
		for (int i=mXSize;i>0;i--) {
			if (getT(0,i-1,j) > mFOTemp) {
				double mTR = (mFOTemp-getT(0,i-1,j))/(getT(0,i,j)-getT(0,i-1,j));
				mFOS[foSize][0] = getX(0,i-1,j,1) + (getX(0,i,j,1)-getX(0,i-1,j,1))*mTR;
				mFOS[foSize][1] = getX(0,i-1,j,2);
				
				CCell foCell(getTau(),mFOS[foSize][0],mFOS[foSize][1],0.);
				
				for (int k=0;k<11;k++)
					foCell.setS(k, getS(0,i-1,j,k) + (getS(0,i,j,k)-getS(0,i-1,j,k))*mTR);

#ifdef HYDRO_BOOST_THREADS
				foCell.setEos(eosVector[0]);
				foCell.setHelper(helperVector[0]);
#endif
				
				foCell.selfUpdate();
				
				for (int k=1;k<3;k++) 
					mFOVelo[foSize][k-1] = foCell.getS(k);
				
				mFODiss[foSize][0] = 0.5*(foCell.getPixyLocal(1,1)-foCell.getPixyLocal(2,2));
				mFODiss[foSize][1] = 0.5*(foCell.getPixyLocal(1,1)+foCell.getPixyLocal(2,2)
										  -2*pow(foCell.getTau(),2)*foCell.getPixyLocal(3,3))/ROOT3;
				mFODiss[foSize][2] = foCell.getPixyLocal(1,2);
				mFODiss[foSize][3] = foCell.getTau()*foCell.getPixyLocal(1,3);
				mFODiss[foSize][4] = foCell.getTau()*foCell.getPixyLocal(2,3);
				mFODiss[foSize][5] = foCell.getB();
				
				if (j==0) {
					mFOSigma[foSize][0] = 1.;
					mFOSigma[foSize][1] = 0.;
				}

				if (foSize > 1){
					v = (mFOS[foSize-2][0] - mFOS[foSize][0])/(mFOS[foSize][1] - mFOS[foSize-2][1]);
					mFOSigma[foSize-1][1] = v/sqrt(1.+v*v);
					mFOSigma[foSize-1][0] = 1./sqrt(1.+v*v);
				}
				
				foSize++;
				break;
			}
		}
		j--;
	}
	
	if (!goodFOS(mFOS,foSize)) 
		return false;
	
	return true;
}

	// freezeout surface at eta
bool octMesh::getFOS(int eta, double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], 
					 double mFOVelo[XSIZE+YSIZE][2], double mFODiss[XSIZE+YSIZE][6], int &foSize) {
	
	if (getT(eta,0,0) < mFOTemp)
		return false;
	else if (!getActive(eta,0,0)){
		printf("CMesh::getFOS(eta) found inactive before below FOtemp??\n\n");
		return false;
	}
	
	
	foSize = 0;
	int j = mYSize;  //keep track of where we are in j
	double v=0.;
	
	for (int i=0; i<=mXSize; i++) {
		if (j+3 <= mYSize) j+=3;
		else j=mYSize;
		
		while (j>0) {
			if (getT(eta,i,j-1) > mFOTemp) {
				
				double mTR = (mFOTemp-getT(eta,i,j-1))/(getT(eta,i,j)-getT(eta,i,j-1));
				
				mFOS[foSize][1] = getX(eta,i,j-1,2) + (getX(eta,i,j,2) - getX(eta,i,j-1,2))*mTR;
				mFOS[foSize][0] = getX(eta,i,j-1,1);
				
				CCell foCell(getTau(),mFOS[foSize][0],mFOS[foSize][1],0.);
				
				for (int k=0;k<11;k++)
					foCell.setS(k,getS(eta,i,j-1,k) + (getS(eta,i,j,k)-getS(eta,i,j-1,k))*mTR);
				
#ifdef HYDRO_BOOST_THREADS
				foCell.setEos(eosVector[0]);
				foCell.setHelper(helperVector[0]);
#endif
				foCell.selfUpdate();
				
				for (int k=1;k<3;k++) 
					mFOVelo[foSize][k-1] = foCell.getS(k);
				
				mFODiss[foSize][0] = 0.5*(foCell.getPixyLocal(1,1)-foCell.getPixyLocal(2,2));
				mFODiss[foSize][1] = 0.5*(foCell.getPixyLocal(1,1)+foCell.getPixyLocal(2,2)
										  -2*pow(foCell.getTau(),2)*foCell.getPixyLocal(3,3))/ROOT3;
				mFODiss[foSize][2] = foCell.getPixyLocal(1,2);
				mFODiss[foSize][3] = foCell.getTau()*foCell.getPixyLocal(1,3);
				mFODiss[foSize][4] = foCell.getTau()*foCell.getPixyLocal(2,3);
				mFODiss[foSize][5] = foCell.getB();
				
				if (foSize==0) {
					v=0.;
					mFOSigma[foSize][0] = 0.;
					mFOSigma[foSize][1] = 1.;
				}
				else if (foSize>1) {
					v = (mFOS[foSize-2][1] - mFOS[foSize][1])/(mFOS[foSize][0] - mFOS[foSize-2][0]);
					mFOSigma[foSize-1][0] = v/sqrt(1.+v*v);
					mFOSigma[foSize-1][1] = 1./sqrt(1.+v*v);
				}
				foSize++;
				break;
			}
			else j--;
		}
		if (v>1.) 
			break;
		
		if (j==0 && i==0){
			foSize=0;
			return false;
		}
	}
	
	j--;
	
	if (j==-2)
		j = (int)(mFOS[foSize-1][1]/mDy);
	
	while(j>=0) {
		for (int i=mXSize;i>0;i--) {
			if (getT(eta,i-1,j) > mFOTemp) {
				double mTR = (mFOTemp-getT(eta,i-1,j))/(getT(eta,i,j)-getT(eta,i-1,j));
				mFOS[foSize][0] = getX(eta,i-1,j,1) + (getX(eta,i,j,1)-getX(eta,i-1,j,1))*mTR;
				mFOS[foSize][1] = getX(eta,i-1,j,2);
				
				CCell foCell(getTau(),mFOS[foSize][0],mFOS[foSize][1],0.);
				
				for (int k=0;k<11;k++)
					foCell.setS(k, getS(eta,i-1,j,k) + (getS(eta,i,j,k)-getS(eta,i-1,j,k))*mTR);
				
#ifdef HYDRO_BOOST_THREADS
				foCell.setEos(eosVector[0]);
				foCell.setHelper(helperVector[0]);
#endif
				
				foCell.selfUpdate();
				
				for (int k=1;k<3;k++) 
					mFOVelo[foSize][k-1] = foCell.getS(k);
				
				mFODiss[foSize][0] = 0.5*(foCell.getPixyLocal(1,1)-foCell.getPixyLocal(2,2));
				mFODiss[foSize][1] = 0.5*(foCell.getPixyLocal(1,1)+foCell.getPixyLocal(2,2)
										  -2*pow(foCell.getTau(),2)*foCell.getPixyLocal(3,3))/ROOT3;
				mFODiss[foSize][2] = foCell.getPixyLocal(1,2);
				mFODiss[foSize][3] = foCell.getTau()*foCell.getPixyLocal(1,3);
				mFODiss[foSize][4] = foCell.getTau()*foCell.getPixyLocal(2,3);
				mFODiss[foSize][5] = foCell.getB();
				
				if (j==0) {
					mFOSigma[foSize][0] = 1.;
					mFOSigma[foSize][1] = 0.;
				}
				
				if (foSize > 1){
					v = (mFOS[foSize-2][0] - mFOS[foSize][0])/(mFOS[foSize][1] - mFOS[foSize-2][1]);
					mFOSigma[foSize-1][1] = v/sqrt(1.+v*v);
					mFOSigma[foSize-1][0] = 1./sqrt(1.+v*v);
				}
				
				foSize++;
				break;
			}
		}
		j--;
	}
	
	if (!goodFOS(mFOS,foSize)) 
		return false;
	
	return true;
}

bool octMesh::goodFOS(double mFOS[XSIZE+YSIZE][2], int foSize) {
	double lastPhi = atan2(mFOS[1][1],mFOS[1][0]);
	for (int i=2;i<foSize;i++){
		double newPhi = atan2(mFOS[i][1],mFOS[i][0]);
			//		printf("%g ?< %g  (%g,%g)\n",lastPhi,newPhi,mFOS[i][1],mFOS[i][0]);
		if (lastPhi < newPhi)
			return false;
		lastPhi = newPhi;
	}
	return true;
}

int octMesh::genFOS(CMesh *mMesh) {
	
#ifdef HYDRO_XDMF_OUTPUT
	static int genFOSCalls = 0;
#endif
	
	if (!mPureBjorken) {
		double dTau = getTau() - mMesh->getTau();
		double fosE = getEos()->getEGivenT(mFOTemp);
		double fosP = getEos()->getP(fosE);
		static CCell mC(-1.,-1.,-1.,-1.);
		
		dxCorn[0] = dTau;
		for (int i=1;i<4;i++) 
			dxCorn[i] = mC.dx[i];
		
		mCornelius.init(4,fosE,dxCorn);
		
#ifdef HYDRO_BOOST_THREADS
		mC.setEos(eosVector[0]);
		mC.setHelper(helperVector[0]);
#endif
		
		for (int i=0;i<mNSize;i++)
			for (int j=0;j<mXSize;j++)
				for (int k=0;k<mYSize;k++){					
					CCell *cubeCells[2][2][2][2];
					
					for (int m=0;m<2;m++)
						for (int n=0;n<2;n++)
							for (int p=0;p<2;p++)
								cubeCells[0][m][n][p] = mMesh->getCell(i+m,j+n,k+p);
					for (int m=0;m<2;m++)
						for (int n=0;n<2;n++)
							for (int p=0;p<2;p++)
								cubeCells[1][m][n][p] = getCell(i+m,j+n,k+p);
					
					bool mSnap = false;
					for (int m=0;m<2;m++)
						for (int n=0;n<2;n++)
							for (int p=0;p<2;p++)
								for (int q=0;q<2;q++)
									if (!cubeCells[m][n][p][q]->getActive())
										mSnap = true;
#ifdef HYDRO_BOOST_THREADS
									else  {
										cubeCells[m][n][p][q]->setEos(eosVector[0]);
										cubeCells[m][n][p][q]->setHelper(helperVector[0]);
									}
#endif
					if (mSnap) break;
					
					for (int i1=0; i1 < 2; i1++) 
						for (int i2=0; i2 < 2; i2++) 
							for (int i3=0; i3 < 2; i3++) 
								for (int i4=0; i4 < 2; i4++) 
									localE4[i1][i2][i3][i4] = cubeCells[i1][i2][i3][i4]->getE();
					
					if (!containsSurf(fosE,localE4)) 
						continue;
					
					double foPoint[4];
					double surfVec[4];
					mCornelius.find_surface_4d(localE4);
					
					if (mCornelius.get_Nelements() > 1)
						printf("\n ********* MULTI_ELEMENT CELL (%g %g %g %g) ************\n\n",
							   cubeCells[0][0][0][0]->getTau(),cubeCells[0][0][0][0]->getX(1),
							   cubeCells[0][0][0][0]->getX(2) ,cubeCells[0][0][0][0]->getX(3));

					for (int elem=0; elem < mCornelius.get_Nelements(); elem++) {
						for (int dim=0; dim < 4; dim++) {
							surfVec[dim]  = mCornelius.get_normal_elem(elem,dim);
							
							foPoint[dim]  = mCornelius.get_centroid_elem(elem,dim);
							foPoint[dim] += cubeCells[0][0][0][0]->getX(dim);
						}
						
						mC.setX(foPoint);
						
						double foS[11];
						hcubeInterp(foS, foPoint, cubeCells, dTau);
						mC.setS(foS);
						mC.selfUpdate();

						/*
						printf("%g %g %g %g : %g %g %g %g\n",
							   cubeCells[0][0][0][0]->getTau(),cubeCells[0][0][0][0]->getX(1),
							   cubeCells[0][0][0][0]->getX(2), cubeCells[0][0][0][0]->getX(3),
							   dxCorn[0],dxCorn[1],dxCorn[2],dxCorn[3]);
						for (int m=0;m<2;m++) {
							for (int n=0;n<2;n++)
								for (int p=0;p<2;p++)
									for (int q=0;q<2;q++)
										printf("%g ",localE[m][n][p][q]-fosE);
							printf("\n");
						}
						printf("\n%g %g %g %g\n\n",mC.getTau(),mC.getX(1),mC.getX(2),mC.getX(3));
						 */
						 
						fprintf(osuFOS,"%g %g %g %g ",mC.getTau(),mC.getX(1),mC.getX(2),mC.getX(3)); //1-4
						fprintf(osuFOS,"%g %g %g %g ",mC.getTau()*surfVec[0],
								-mC.getTau()*surfVec[1],-mC.getTau()*surfVec[2],-surfVec[3]/mC.getTau());     //5-8
						fprintf(osuFOS,"%g %g %g ",mC.getUx(),mC.getUy(),mC.getUz());                //9-11
						fprintf(osuFOS,"%g 0. %g 0. 0. %g ",fosE,mFOTemp,fosP);                      //12-17
						
							// 18-21, 22-24, 25-26, 27
						for (int i=0;i<4;i++)
							for (int j=i;j<4;j++)
								fprintf(osuFOS,"%g ",mC.getPixy(i,j));
						
#ifdef HYDRO_XDMF_OUTPUT
						fprintf(osuFOS,"%g %d",mC.getB(),genFOSCalls); // 28-29
#else
						fprintf(osuFOS,"%g",mC.getB()); // 28
#endif
						fprintf(osuFOS,"\n");
					
					}
				}
					
#ifdef HYDRO_XDMF_OUTPUT
		genFOSCalls++;
#endif
		
		return 0;
	}
	else {
		double dTau = getTau() - mMesh->getTau();
		double fosE = getEos()->getEGivenT(mFOTemp);
		double fosP = getEos()->getP(fosE);
	
		static CCell mC(0.,0.,0.,0.);
		
#ifdef HYDRO_BOOST_THREADS
		mC.setEos(eosVector[0]);
		mC.setHelper(helperVector[0]);
#endif
		
		int nSurf=0; int *nSurfptr = &nSurf;
		int nErr=0;  int *nErrptr = &nErr;
	
		for (int i=0;i<mXSize;i++)
			for (int j=0;j<mYSize;j++) {
				
				CCell *cube[8];
				cube[0] = mMesh->getCell(0,i,j);
				cube[1] = mMesh->getCell(0,i+1,j);
				cube[2] = mMesh->getCell(0,i+1,j+1);
				cube[3] = mMesh->getCell(0,i,j+1);
				cube[4] = getCell(0,i,j);
				cube[5] = getCell(0,i+1,j);
				cube[6] = getCell(0,i+1,j+1);
				cube[7] = getCell(0,i,j+1);
				
					// only among active cells
				for (int k=0;k<8;k++)
					if (!cube[k]->getActive())
						break;
				
					// fill energy density of cells
				double localE[8];
				for (int k=0;k<8;k++)
					localE[k] = cube[k]->getE();
				
				if (!containsSurf(fosE,localE)) 
					continue;
				
				double surfVec[3], midX[3];
				
				int cubeOut = cubical2p1(fosE, localE, surfVec, nSurfptr,
										 midX, dTau, mDx, mDy, nErrptr);
				
				if (*nErrptr != 0 && false){
					printf("errors? (%d)\n",*nErrptr);
				}
				
				double mS[11];
				cubeInterp(mS, midX, cube, dTau);
				mC.setS(mS);
				mC.setTau(cube[0]->getTau()+midX[0]);
				mC.setX(cube[0]->getX(1)+midX[1]);
				mC.setY(cube[0]->getX(2)+midX[2]);
				mC.selfUpdate();
				
					//printf("\n\nMc print::\n");
					//mC.selfPrint();
					//mC.verify();
				
				fprintf(osuFOS,"%g %g %g %g %g %g ",mC.getTau(),mC.getX(1),mC.getX(2),surfVec[0],-surfVec[1],-surfVec[2]); 
				fprintf(osuFOS,"%g %g %g 0. %g 0. 0. %g ",mC.getVx(),mC.getVy(),fosE,mFOTemp,fosP); // 7-14
				/*
				fprintf(osuFOS,"%g %g %g %g %g %g %g\n",getTau()*getTau()*mC.getPixy(3,3),
						mC.getPixy(0,0),mC.getPixy(0,1),mC.getPixy(0,2),
						mC.getPixy(1,1),mC.getPixy(1,2),mC.getPixy(2,2)); //15-22
				 */
				fprintf(osuFOS,"%g %g %g %g %g %g %g\n",getTau()*getTau()*mC.getPixyBulk(3,3),
						mC.getPixyBulk(0,0),mC.getPixyBulk(0,1),mC.getPixyBulk(0,2),
						mC.getPixyBulk(1,1),mC.getPixyBulk(1,2),mC.getPixyBulk(2,2)); //15-22
			}
		if (nErr > 1)
			printf("\ngenFOS errors (%d)... continuing...\n",nErr);
		return *nErrptr;
	}
}

void octMesh::checkAzimuthalSymmetry() {
	double tolerance = 1E-10;
	
	for (int i=0;i<=mNSize;i++) 
		for (int j=0;j<=mXSize;j++)
			for (int k=j;k<=mYSize;k++) {
				if (!getActive(i,j,k)) break;
				if ( abs(getS(i,j,k,1) - getS(i,k,j,2)) > tolerance*abs(getS(i,j,k,1)))
					printf("velo sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,1),getS(i,k,j,2));
				if ( abs(getS(i,j,k,2) - getS(i,k,j,1)) > tolerance*abs(getS(i,j,k,2)))
					printf("velo sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,2),getS(i,k,j,1));
				if ( abs(getS(i,j,k,4) - getS(i,k,j,4)) > tolerance*abs(getS(i,j,k,4)))
					printf("energy sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,4),getS(i,k,j,4));
				
				if ( abs(getS(i,j,k,5) + getS(i,k,j,5)) > tolerance*abs(getS(i,j,k,5)) && k!=j)
					printf("a1 sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,5),-getS(i,k,j,5));
				
				if ( abs(getS(i,j,k,6) - getS(i,k,j,6)) > tolerance*abs(getS(i,j,k,6)))
					printf("a2 sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,6),getS(i,k,j,6));
				if ( abs(getS(i,j,k,7) - getS(i,k,j,7)) > tolerance*abs(getS(i,j,k,7)))
					printf("a3 sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,7),getS(i,k,j,7));
				if ( abs(getS(i,j,k,8) - getS(i,k,j,9)) > tolerance*abs(getS(i,j,k,8)))
					printf("a4/a5 sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,8),getS(i,k,j,9));
			}
}

void octMesh::writeIEPvEta(){
	string fn = mDataRoot;
	fn.append("/fIEPvEta.txt");
	FILE *fIEPvEta = fopen(fn.c_str(),"w");
	
	double iEP[mNSize], iEPD[mNSize];
	
	for (int i=0;i<mNSize;i++) {
		iEP[i]  = 0.;
		iEPD[i] = 0.;
		for (int j=0;j<mXSize;j++)
			for (int k=0;k<mYSize;k++) {
				if (!getActive(i,j,k)) break;
#ifdef HYDRO_BOOST_THREADS
				getCell(i,j,k)->setEos(eosVector[0]);
				getCell(i,j,k)->setHelper(helperVector[0]);
#endif
				selfUpdate(i,j,k);
				double mTxx = getTxy(1,1);
				double mTyy = getTxy(2,2);
				iEP[i]  += mTxx - mTyy;
				iEPD[i] += mTxx + mTyy;
			}
		fprintf(fIEPvEta,"%g %g\n",getX(i,0,0,3),iEP[i]/iEPD[i]);
	}
	fclose(fIEPvEta);
	
	fn = mDataRoot;
	fn.append("/fIEP_ideal_vEta.txt");
	fIEPvEta = fopen(fn.c_str(),"w");
	
	for (int i=0;i<=mNSize;i++) {
		iEP[i]  = 0.;
		iEPD[i] = 0.;
		for (int j=0;j<=mXSize;j++)
			for (int k=0;k<=mYSize;k++) {
				if (!getActive(i,j,k)) break;
#ifdef HYDRO_BOOST_THREADS
				getCell(i,j,k)->setEos(eosVector[0]);
				getCell(i,j,k)->setHelper(helperVector[0]);
#endif
				selfUpdate(i,j,k);
				double mTxx = getTxy(1,1) - getPixy(1,1);
				double mTyy = getTxy(2,2) - getPixy(2,2);
				iEP[i]  += mTxx - mTyy;
				iEPD[i] += mTxx + mTyy;
			}
		fprintf(fIEPvEta,"%g %g\n",getX(i,0,0,3),iEP[i]/iEPD[i]);
	}
	fclose(fIEPvEta);
}

double octMesh::getdNd3p(double px, double py, double Y, double m) {
	double mTemp = 0.;
	double pz,T,E,P,p0,pu,puu;
		//  p0 = sqrt(m*m + px*py + py*py)*cosh(Y);
		//  pz = p0 * tanh(Y);
	
		//  printf("%g %g %g %g:\n\n",px,py,Y,m);
	
	for (int eta=0;eta<=mNSize;eta++) {
		p0 = sqrt(m*m + px*py + py*py) * cosh(Y - mDn*eta);
		pz = p0 * tanh(Y - mDn*eta);
		
			//	printf("\n%d %g %g:\n\n",eta,p0,pz);
		
		for (int i=0;i<=mXSize;i++)
			for (int j=0;j<=mYSize;j++) {
				pu = p0*getS(eta,i,j,0) - px*getS(eta,i,j,1) - py*getS(eta,i,j,2) - pz*getS(eta,i,j,3);
				
				puu = getPixy(eta,i,j,0,0)*p0*p0;
				puu -= getPixy(eta,i,j,0,1)*p0*px + getPixy(eta,i,j,0,2)*p0*pz + getPixy(eta,i,j,0,3)*p0*pz;
				puu +=  getPixy(eta,i,j,1,1)*px*px +getPixy(eta,i,j,2,2)*py*py +getPixy(eta,i,j,3,3)*pz*pz;
				puu += 2.*(getPixy(eta,i,j,1,2)*px*py +getPixy(eta,i,j,1,3)*px*pz +getPixy(eta,i,j,2,3)*py*pz);
				
				T = getT(eta,i,j);
				E = getE(eta,i,j);
				P = getP(eta,i,j);
					//		printf("%d %d:  %g %g %g\n",i,j,exp(-pu/T),0.5*puu/((E+P)*T*T),cosh(Y - mDn*eta));
				mTemp += p0 * exp( - pu/T) * (1. + 0.5*puu/((E+P)*T*T));// * cosh(Y - mDn*eta);
			}
	}
	
	return mTemp*mDx*mDy*mDn*getTau()*JIFPI*JIFPI*JIFPI*8.;
}

	// for fixed tau freeze-out
void octMesh::calcdNdY(double m, double d, double *output, double dRap, int nRap){
	if (mPureBjorken) 
		return;
	
	double temp = 0.;
	
	mFBose = true;
	mFFerm = false;
	
	double dPt = 0.2;  int nPt = 15;
	double dPhi = (JPI)/(16.);  int nPhi = 8;
	
	double dNd3p[nPt][nPhi+1][nRap];
	double c = d*getTau()*pow(2.*JPI*JHBARC,-3);
	
	for (int iRap=0;iRap<nRap;iRap++) 
		for (int iPt=0; iPt<nPt; iPt++)
			for (int iPhi=0;iPhi<=nPhi;iPhi++)  {
				dNd3p[iPt][iPhi][iRap] = 0.;
				for (int iN= 0;iN<=mNSize;iN++)
					for (int iX= 0;iX<=mXSize;iX++)
						for (int iY= 0;iY<=mYSize;iY++) {
							if (!getActive(iN,iX,iY)) break;
							dNd3p[iPt][iPhi][iRap] += c*getF(iX, iY, iN, dPt*iPt, dPhi*iPhi, dRap*iRap, m);
						}
			}
	
	for (int iRap=0;iRap<nRap;iRap++) {
		temp = 0.;
		for (int iPt=0; iPt<nPt; iPt++){
			for (int iPhi=0;iPhi<=nPhi;iPhi++) 
				temp += 4. * dPt * dPhi * (dPt*iPt) * dNd3p[iPt][iPhi][iRap];
			temp -= 2.* dPt * dPhi * (dPt*iPt) * (dNd3p[iPt][0][iRap] + dNd3p[iPt][nPhi][iRap]);
		}
		output[iRap] = temp;
	}
	
	string fn = "fdNd3p.txt";
	fn = mDataRoot + fn;
	FILE *fdNd3p = fopen( fn.c_str(),"w");
	fprintf(fdNd3p,"%d %g %d %g %d %g\n",nRap,dRap,nPt,dPt,nPhi,dPhi);
	
	for (int iRap=0;iRap<nRap;iRap++) 
		for (int iPt=0; iPt<nPt; iPt++)
			for (int iPhi=0;iPhi<=nPhi;iPhi++) 
				fprintf(fdNd3p,"%g %g %g %g\n",dRap*iRap,dPt*iPt,dPhi*iPhi,dNd3p[iPt][iPhi][iRap]);
	
}

double octMesh::getF( int x, int y, int n, double pt, double phi, double rap, double m) {
	double f = 0.;
	
	double dV = mDx*mDy*mDn;
	
	double T = getT(n,x,y);
	double mT = sqrt(m*m + pt*pt);
	
	double u[4];
	for (int i=1;i<4;i++)  u[i] = getS(n,x,y,i);
	u[0] = sqrt(1.+u[1]*u[1]+u[2]*u[2]+pow(getTau()*u[3],2));
	
	double mX[4];
	for (int i=1;i<4;i++) mX[i] = getX(n,x,y,i);
	mX[0] = getTau();
	
	double p[4];
	p[0] = mT*cosh(rap - mX[3]);
	p[1] = pt*cos(phi);
	p[2] = pt*sin(phi);
	p[3] = mT*sinh(rap - mX[3]);
	
	double diffTotal = 0.;
	
		// always do the positive
	double uDotp = u[0]*p[0] - u[1]*p[1] - u[2]*p[2] - mX[0]*u[3]*p[3];
	double dF = 0.;
	
	f += p[0]/mT * CMesh::getF(uDotp/T,dF); //Q1
	diffTotal += dV;
	
		//	return 8.*mT*f*diffTotal;
	
	if (x>0) {
		mX[1] = -mX[1];
		u[1] = -u[1];
		uDotp = u[0]*p[0] - u[1]*p[1] - u[2]*p[2] - mX[0]*u[3]*p[3];
		f += p[0]/mT * CMesh::getF(uDotp/T,dF); //Q2
		diffTotal += dV;
		
		if (y>0) {
			mX[2] = -mX[2];
			u[2] = -u[2];
			uDotp = u[0]*p[0] - u[1]*p[1] - u[2]*p[2] - mX[0]*u[3]*p[3];
			f += p[0]/mT * CMesh::getF(uDotp/T,dF); //Q3
			diffTotal += dV;
			
			if (mX[3] > 0.) {
				mX[3] = -mX[3];
				u[3] = -u[3];
				p[0] = mT*cosh(rap - mX[3]);
				p[3] = mT*sinh(rap - mX[3]);
				
				uDotp = u[0]*p[0] - u[1]*p[1] - u[2]*p[2] - mX[0]*u[3]*p[3];
				f += p[0]/mT * CMesh::getF(uDotp/T,dF); //Q4
				diffTotal += dV;
				
				mX[3] = -mX[3];
				u[3] = -u[3];
			}
			mX[2] = -mX[2];
			u[2] = -u[2];
		}
		
		if (n>0) {
			mX[3] = -mX[3];
			u[3] = -u[3];
			p[0] = mT*cosh(rap - mX[3]);
			p[3] = mT*sinh(rap - mX[3]);
			uDotp = u[0]*p[0] - u[1]*p[1] - u[2]*p[2] - mX[0]*u[3]*p[3];
			f += p[0]/mT * CMesh::getF(uDotp/T,dF); //Q5
			diffTotal += dV;
			
			mX[3] = -mX[3];
			u[3] = -u[3];
		}
		
		mX[1] = -mX[1];
		u[1] = -u[1];
	}
	
	if (y>0) {
		mX[2] = -mX[2];
		u[2] = -u[2];
		p[0] = mT*cosh(rap - mX[3]);
		p[3] = mT*sinh(rap - mX[3]);
		uDotp = u[0]*p[0] - u[1]*p[1] - u[2]*p[2] - mX[0]*u[3]*p[3];
		f += p[0]/mT * CMesh::getF(uDotp/T,dF); //Q6
		
		diffTotal += dV;
		
		if (n>0) {
			mX[3] = -mX[3];
			u[3] = -u[3];
			p[0] = mT*cosh(rap - mX[3]);
			p[3] = mT*sinh(rap - mX[3]);
			uDotp = u[0]*p[0] - u[1]*p[1] - u[2]*p[2] - mX[0]*u[3]*p[3];
			f += p[0]/mT * CMesh::getF(uDotp/T,dF); //Q7
			diffTotal += dV;
			
			mX[3] = -mX[3];
			u[3] = -u[3];
		}
		mX[2] = -mX[2];
		u[2] = -u[2];
	}
	
	if (n>0) {
		mX[3] = -mX[3];
		u[3] = -u[3];
		p[0] = mT*cosh(rap - mX[3]);
		p[3] = mT*sinh(rap - mX[3]);
		uDotp = u[0]*p[0] - u[1]*p[1] - u[2]*p[2] - mX[0]*u[3]*p[3];
		f += p[0]/mT * CMesh::getF(uDotp/T,dF); //Q8
		diffTotal += dV;
		
		mX[3] = -mX[3];
		u[3] = -u[3];
	}
	
	return mT*f*diffTotal;
}
