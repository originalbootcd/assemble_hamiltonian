CXX			= icpc
CXXFLAGS	= -std=c++11 -g
INCLUDE		= -I$(HOME)/anton/boost_1_49_0 -I$(HOME)/include
LIBS		= -L$(HOME)/lib -lconfig++

all:
	@$(CXX) main.cpp -o main $(INCLUDE) $(CXXFLAGS) $(LIBS)
