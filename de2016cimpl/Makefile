CPP	= g++
CPPFLAGS	= -O3 -std=c++11 -Wall -fexceptions -msse4.2 -fopenmp
LFLAGS	= -lgomp -pthread
INCPATH	= -I../boost_1_61_0 -I../boost.simd/include

DEPS = setup.h molecule.h solver.h neighbourhood.h stats.h stratify.h structural_similarity.h descriptor.h power_spectrum.h local_similarity.h run.h
OBJ = main.o setup.o molecule.o solver.o neighbourhood.o stats.o stratify.o structural_similarity.o descriptor.o power_spectrum.o local_similarity.o run.o

COBJ = soap_c_wrap.o
CDEPS = soap_c_wrap.h

default: main.exe

%.o: %.cpp $(DEPS)
	$(CPP) $(CPPFLAGS) $(INCPATH) -c -o $@ $<

main.exe: $(OBJ)
	$(CPP) -o $@ $^ $(LFLAGS)

mainc.exe: $(OBJ) mainc.c soap_c_wrap.cpp soap_c_wrap.h
	$(CPP) $(CPPFLAGS) $(INCPATH) -c soap_c_wrap.cpp -o soap_c_wrap.o
	gcc -c mainc.c -o mainc.o
	g++ mainc.o soap_c_wrap.o $(OBJ) -o mainc.exe

clean:
	rm -rf *.o *.exe
