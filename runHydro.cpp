#include "CEos.h"
#include "CMesh.h"
#include "CHydro.h"

int main (int argc, char * const argv[]) {

	parameterMap* pMap = new parameterMap();
	for (int i=1;i<argc;i++) 
		parameter::ReadParsFromFile(*pMap, argv[i]);
		
		//parameter::PrintPars(*pMap);
	CHydro* mHydro = new CHydro(pMap);
	int status = mHydro->runHydro();
	delete mHydro;

	return status;
	
	/*
	CMesh* m1 = new octMesh(pMap);
	CMesh* m2 = new octMesh(m1);
	
	CEos* mEos = new CEos(pMap);
    m1->setEos(mEos);
	m2->setTau(m1->getTau()+0.1);
	
	FILE* fos = fopen("output/osuFOS.dat","w");
	m1->setOsuFos(fos);
	m2->genFOS(m1);
	
	return 0;
	 */
}

