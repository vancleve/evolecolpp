CXXFLAGS = -fPIC -Wall -Wextra -I. `python-config --includes` -std=c++11
LDFLAGS = `python-config --ldflags` -lboost_python -lboost_system -lgsl -lgslcblas 

ifdef FWDPP
	CXXFLAGS += -I$(FWDPP)
endif

nodebug: cooperation_snowdrift cooperation_snowdrift_metapop
nodebug: CXXFLAGS += -O3 -DNDEBUG

debug: cooperation_snowdrift cooperation_snowdrift_metapop
debug: CXXFLAGS += -g -O0

cooperation_snowdrift: cooperation_snowdrift.o
	$(CXX) -shared -o cooperation_snowdrift.so cooperation_snowdrift.o $(LDFLAGS)

cooperation_snowdrift_metapop: cooperation_snowdrift_metapop.o
	$(CXX) -shared -o cooperation_snowdrift_metapop.so cooperation_snowdrift_metapop.o $(LDFLAGS)

clean:
	rm -f *.o *.so
