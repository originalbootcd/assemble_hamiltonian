CXX			= g++
CXXFLAGS	= -std=c++11 -g
INCLUDE		= ../../../boost_1_49_0


all:
	@$(CXX) main.cpp -o main -I$(INCLUDE) $(CXXFLAGS)