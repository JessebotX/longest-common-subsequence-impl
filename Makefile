#compiler setup
CXX = g++
MPICXX = mpic++
CXXFLAGS = -std=c++14 -O3

SERIAL = lcs_serial
PARALLEL = lcs_parallel
DISTRIBUTED = lcs_distributed
ALL= $(SERIAL) $(PARALLEL) $(DISTRIBUTED)

all : $(ALL)

$(SERIAL): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(PARALLEL): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(DISTRIBUTED): %: %.cpp
	$(MPICXX) $(CXXFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)
