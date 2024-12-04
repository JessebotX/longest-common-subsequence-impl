#compiler setup
CXX = g++
MPICXX = mpic++
CXXFLAGS = -std=c++14 -pthread -O3

SERIAL = lcs_serial
PARALLEL = lcs_parallel_test
DISTRIBUTED = lcs_distributed
ALL= $(SERIAL) $(PARALLEL) $(DISTRIBUTED)
COMMON= common/utils.h common/cxxopts.h common/get_time.h

all : $(ALL)

$(SERIAL): %: %.cpp $(COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(PARALLEL): %: %.cpp $(COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(DISTRIBUTED): %: %.cpp $(COMMON)
	$(MPICXX) $(CXXFLAGS) -o $@ $<

$(COMMON):

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)
