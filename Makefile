CXXFLAGS=-fPIC -Wall -W -O3 -I. -I/Users/vancleve/science/code/fwdpp `python-config --includes` -std=c++11 
LDFLAGS=`python-config --ldflags` -lboost_python -lboost_system -lgsl -lgslcblas 

all: social_evol.o
	$(CXX) -shared -o social_evol.so social_evol.o $(LDFLAGS)

