# enter make, make runHydro, make install, make clean, make uninstall
# local executables can be made with make runHydro or make multiEos
include ../../makefile_defs.mk
### the variables below are set by the include statement above
CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS} 
CORAL_INCLUDE=${INSTALLDIR}/include
#e.g. /Users/scottepratt/git -- 
CORAL_HOME=${MADAI_HOME}/rhic/CorAL/trunk
#e.g. ${MADAI_HOME}/rhic/CorAL/trunk
INSTALLDIR=${MADAI_INSTALLDIR}
#e.g. /Users/scottepratt/local
GSLPATH=${MADAI_GSLPATH}
#e.g. /opt/local
#########################################################################
#HOPEFULLY, everything below need not be touched

INC=-I ${GSLPATH}/include -Iinclude -I${CORAL_INCLUDE} -I${MADAI_HDF5_HOME}/include -I${MADAI_XDMF_HOME}/include
LIB=-Llib -L${GSLPATH}/lib -L${INSTALLDIR}/lib -L${MADAI_XDMF_HOME}/lib

ISHYDRO3_HFILES = include/CHydro.h include/CEos.h include/CCell.h include/CMesh.h include/hydroDef.h include/fullMesh.h include/octMesh.h include/cornelius.h
ISHYDRO3_OBJFILES = build/CEos.o build/CCell.o build/CMesh.o build/octMesh.o build/fullMesh.o build/CHydro.o build/cornelius.o

ishydro3 : ishydro3_dirs ${ISHYDRO3_HFILES} lib/libishydro3.a

lib/libishydro3.a : ${ISHYDRO3_OBJFILES} ${SURF_OBJFILES}
	rm -f lib/libishydro3.a;\
	ar r lib/libishydro3.a ${ISHYDRO3_OBJFILES} ${SURF_OBJFILES};

####################################
# for development
runHydro : lib/libishydro3.a src/hydroDef.h runHydro.cpp
	${CPP} runHydro.cpp -o runHydro ${OPT} ${INC} -Llib -lishydro3 ${LIB} -lcoralutils #-lhdf5 -lhdf5_cpp -lhdf5_hl -lXdmf -lmpi -lmpi_cxx -lboost_thread-mt

####################################

include/CHydro.h : src/CHydro.h
	cp -f src/CHydro.h include/

include/CEos.h : src/CEos.h
	cp -f src/CEos.h include/

include/CCell.h : src/CCell.h
	cp -f src/CCell.h include/

include/CMesh.h : src/CMesh.h
	cp -f src/CMesh.h include/

include/octMesh.h : src/octMesh.h
	cp -f src/octMesh.h include/
	
include/fullMesh.h : src/fullMesh.h
	cp -f src/fullMesh.h include/

include/hydroDef.h : src/hydroDef.h
	cp -f src/hydroDef.h include/

include/cornelius.h : src/cornelius.h
	cp -f src/cornelius.h include/

######################################

build/CEos.o : ${ISHYDRO3_HFILES} src/CEos.cpp 
	${CPP} ${OPT} ${INC} -c src/CEos.cpp -o build/CEos.o

build/CCell.o : ${ISHYDRO3_HFILES} src/CCell.cpp
	${CPP} -c ${OPT} ${INC} src/CCell.cpp -o build/CCell.o 

build/CMesh.o : ${ISHYDRO3_HFILES} src/CMesh.cpp
	${CPP} -c ${OPT} ${INC} -c src/CMesh.cpp -o build/CMesh.o
	
build/fullMesh.o : ${ISHYDRO3_HFILES} src/fullMesh.cpp
	${CPP} -c ${OPT} ${INC} -c src/fullMesh.cpp -o build/fullMesh.o
	
build/octMesh.o : ${ISHYDRO3_HFILES} src/octMesh.cpp
	${CPP} -c ${OPT} ${INC} -c src/octMesh.cpp -o build/octMesh.o

build/CHydro.o : ${ISHYDRO3_HFILES} src/CHydro.cpp
	${CPP} -c ${OPT} ${INC} src/CHydro.cpp -o build/CHydro.o

build/cornelius.o : ${SURF_HFILES} src/cornelius.cpp
	${CPP} ${INC} -c src/cornelius.cpp -o build/cornelius.o

############################

install : ishydro3 ${ISHYDRO3_HFILES}
	cp -f lib/libishydro3.a ${INSTALLDIR}/lib/;\
	cp -f include/CCell.h ${INSTALLDIR}/include/;\
	cp -f include/CEos.h ${INSTALLDIR}/include/;\
	cp -f include/CHydro.h ${INSTALLDIR}/include/;\
	cp -f include/CMesh.h ${INSTALLDIR}/include/;\
	cp -f include/octMesh.h ${INSTALLDIR}/include;\
	cp -f include/fullMesh.h ${INSTALLDIR}/include;\
	cp -f include/hydroDef.h ${INSTALLDIR}/include/;\
	cp -f include/cornelius.h ${INSTALLDIR}/include/;\
	cp -f -r EosData ${INSTALLDIR}/progdata/

uninstall : 
	rm -f ${INSTALLDIR}/lib/libishydro3.a ${INSTALLDIR}/include/CCell.h ${INSTALLDIR}/include/CEos.h;\
	rm -f ${INSTALLDIR}/include/CHydro.h ${INSTALLDIR}/include/CMesh.h ${INSTALLDIR}/include/hydroDef.h;\
	rm -f ${INSTALLDIR}/include/fullMesh.h ${INSTALLDIR}/include/octMesh.h;\
	rm -f ${INSTALLDIR}/include/cornelius.h;\
	rm -r -f ${INSTALLDIR}/progdata/EosData

clean :
	rm -f include/*.h;\
	rm -f build/*.o;\
	rm -f lib/libishydro3.a

#######################

ishydro3_dirs :
	mkdir -p build;\
	mkdir -p include;\
	mkdir -p lib;\
	mkdir -p ${INSTALLDIR}/include;\
	mkdir -p ${INSTALLDIR}/lib;\
	mkdir -p ${INSTALLDIR}/progdata
