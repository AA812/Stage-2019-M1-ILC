COMPIL=g++ -std=c++11
INCLUDE=`pythia8-config --cxxflags` -I`root-config --incdir`
LIB=`pythia8-config --libs` `root-config --libs`


all : main.exe 

clean:
	rm *.o *.exe

main.exe: PlotEnergy.o main.o 
	$(COMPIL) main.o $(INCLUDE) $(LIB) -o main.exe

main.o: main.cc
	$(COMPIL) $(INCLUDE) main.cc -c -o main.o

PlotEnergy.o: PlotEnergy.cc 
	$(COMPIL) $(INCLUDE) PlotEnergy.cc -c -o PlotEnergy.o



