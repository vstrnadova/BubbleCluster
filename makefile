INC = $(BOOST_DIR)
FLAGS = -std=c++0x #-std=c++11

bc: BubbleClusterV4.cpp MarkersCluster.h MarkersCluster.cpp MCHelpers.h
	g++ -O3 $(FLAGS) -o bc_v4 BubbleClusterV4.cpp
