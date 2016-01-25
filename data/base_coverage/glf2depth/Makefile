CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3
INCS = -L/opt/local/lib/
LIBS = -lz

OBJECTS = glf2depth.o

glf2depth: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCS) -o $@ $(OBJECTS) $(LIBS)
	
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<
	
clean:
	rm -f *.o 