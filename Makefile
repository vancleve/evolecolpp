CXXFLAGS = -fPIC -Wall -Wextra -I. `python-config --includes` -std=c++11
LDFLAGS = `python-config --ldflags` -lboost_python -lboost_system -lgsl -lgslcblas 

ifdef FWDPP
	CXXFLAGS += -I$(FWDPP)
endif

nodebug: cooperation_snowdrift
nodebug: CXXFLAGS += -O3 -DNDEBUG

debug: cooperation_snowdrift
debug: CXXFLAGS += -g -O0

cooperation_snowdrift: cooperation_snowdrift.o
	$(CXX) -shared -o cooperation_snowdrift.so cooperation_snowdrift.o $(LDFLAGS)

clean:
	rm -f cooperation_snowdrift.o cooperation_snowdrift.so
