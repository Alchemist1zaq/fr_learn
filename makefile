CXX           = g++
LINK          = g++

BLAS_INC_DIR = /usr/include
BLAS_LIB_DIR = /usr/lib

CXX_BASE    = -pipe -Wunused-parameter -Wuninitialized -std=c++11 -I./include -I$(TIOGA_INC) $(DEFINES)
CXX_BLAS = -I$(BLAS_INC_DIR) -L$(BLAS_LIB_DIR)
LD_BLAS = -L$(BLAS_LIB_DIR) -L/usr/lib/x86_64-linux-gnu -lopenblas

CXX_BASE += $(CXX_BLAS)
LIBS += $(LD_BLAS)

CXXFLAGS += -g -O2
FFLAGS = -Ofast

#CXXFLAGS += -Wno-unknown-pragmas
CXXFLAGS += -I./include/

OBJECTS_DIR   = ./obj
DESTDIR       = ./bin

OBJECTS = obj/global.o \
		obj/funcs.o \
		obj/points.o \
		obj/matrix.o \
		obj/input.o \
		obj/ele.o \
		obj/polynomials.o \
		obj/operators.o \
		obj/geo.o \
		obj/output.o \
		obj/face.o \
		obj/intFace.o \
		obj/boundFace.o \
		obj/mpiFace.o \
		obj/overFace.o \
		obj/flux.o \
		obj/flurry.o \
		obj/solver.o \
		obj/solver_overset.o \
		obj/multigrid.o \
		obj/superMesh.o \
		obj/overComm.o

TARGET        = Flurry

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

$(TARGET):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(DESTDIR)/$(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS) $(DBG)

clean:
	cd obj && rm -f *.o && cd .. && rm -f bin/Flurry

obj/global.o: src/global.cpp include/global.hpp \
		include/error.hpp \
		include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/global.o src/global.cpp

obj/funcs.o: src/funcs.cpp include/funcs.hpp include/global.hpp include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/funcs.o src/funcs.cpp

obj/points.o: src/points.cpp include/points.hpp include/funcs.hpp include/global.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/points.o src/points.cpp

obj/matrix.o: src/matrix.cpp include/matrix.hpp \
		include/error.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/matrix.o src/matrix.cpp

obj/input.o: src/input.cpp include/input.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/input.o src/input.cpp

obj/ele.o: src/ele.cpp include/ele.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/solver.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/polynomials.hpp \
		include/flux.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ele.o src/ele.cpp

obj/polynomials.o: src/polynomials.cpp include/polynomials.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/polynomials.o src/polynomials.cpp

obj/operators.o: src/operators.cpp include/operators.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/solver.hpp \
		include/ele.hpp \
		include/face.hpp \
		include/polynomials.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/operators.o src/operators.cpp

obj/geo.o: src/geo.cpp include/geo.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/input.hpp \
		include/solver.hpp \
		include/ele.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/polynomials.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/geo.o src/geo.cpp

obj/output.o: src/output.cpp include/output.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/solver.hpp \
		include/ele.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/polynomials.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/output.o src/output.cpp

obj/face.o: src/face.cpp include/face.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/ele.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/solver.hpp \
		include/operators.hpp \
		include/polynomials.hpp \
		include/flux.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/face.o src/face.cpp

obj/intFace.o: src/intFace.cpp include/intFace.hpp  include/face.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/intFace.o src/intFace.cpp

obj/boundFace.o: src/boundFace.cpp include/boundFace.hpp include/face.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/boundFace.o src/boundFace.cpp

obj/mpiFace.o: src/mpiFace.cpp include/mpiFace.hpp include/face.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/mpiFace.o src/mpiFace.cpp

obj/overFace.o: src/overFace.cpp include/overFace.hpp include/face.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/overFace.o src/overFace.cpp

obj/flux.o: src/flux.cpp include/flux.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/input.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/flux.o src/flux.cpp

obj/flurry.o: src/flurry.cpp include/flurry.hpp \
		include/global.hpp \
		include/error.hpp \
		include/matrix.hpp \
		include/input.hpp \
		include/geo.hpp \
		include/solver.hpp \
		include/ele.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/polynomials.hpp \
		include/output.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/flurry.o src/flurry.cpp

obj/solver.o: src/solver.cpp include/solver.hpp \
		include/global.hpp \
		include/error.hpp \
    include/funcs.hpp \
		include/matrix.hpp \
		include/ele.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/face.hpp \
		include/operators.hpp \
		include/overComm.hpp \
		include/polynomials.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/solver.o src/solver.cpp

obj/solver_overset.o: src/solver_overset.cpp include/solver.hpp \
  		include/global.hpp \
		include/input.hpp \
		include/geo.hpp \
		include/matrix.hpp \
		include/overComm.hpp \
		include/overFace.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/solver_overset.o src/solver_overset.cpp

obj/multigrid.o: src/multigrid.cpp include/multigrid.hpp \
	include/input.hpp \
	include/solver.hpp \
	include/ele.hpp \
	include/geo.hpp \
	include/operators.hpp \
	include/matrix.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/multigrid.o src/multigrid.cpp

obj/superMesh.o: src/superMesh.cpp include/superMesh.hpp \
	include/global.hpp \
	include/matrix.hpp \
	include/geo.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/superMesh.o src/superMesh.cpp
	
obj/overComm.o: src/overComm.cpp include/overComm.hpp \
		include/global.hpp \
		include/funcs.hpp \
		include/geo.hpp \
		include/input.hpp \
		include/matrix.hpp \
		include/operators.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/overComm.o src/overComm.cpp