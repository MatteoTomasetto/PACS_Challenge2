CXXFLAGS = -fPIC
LDFLAGS = -L. -Wl,-rpath=${PWD}
LIBS = -lZeroFun

.PHONY: all clean distclean

all: main

main: main.o libZeroFun.so
	$(CXX) $(LDFLAGS) main.o -o main $(LIBS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

libZeroFun.so: ZeroFun.o
	$(CXX) $(LDFLAGS) -shared -Wl,-soname,libZeroFun.so \
	ZeroFun.o -o libZeroFun.so

ZeroFun.o: ZeroFun.cpp
	$(CXX) $(CXXFLAGS) -c ZeroFun.cpp

clean:
	$(RM) *.o 

distclean: clean
	$(RM) libZeroFun.so main
