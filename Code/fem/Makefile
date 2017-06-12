include Makefile.inc

DIRS := errorestimator fem tests tools
SRCS := triangle.cpp basis.cpp edge.cpp element.cpp \
				elementfinder.cpp elementtree.cpp matchable.cpp math.cpp matrix.cpp \
				partition.cpp refinable.cpp tritype.cpp triclass.cpp vertex.cpp \
				solvable.cpp system.cpp reader.cpp approximator.cpp nearbest.cpp \
				dofs.cpp piecewisepolynomial.cpp poly.cpp dofhandler.cpp \
				elementmatrices.cpp errors.cpp
LIBS := 
BINS := 

include $(patsubst %, %/module.mk, $(DIRS))

OBJS := $(patsubst %.cpp, %.o, $(SRCS))

all: $(BINS)

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@ $(INC)

$(BINS): %: %.o $(OBJS)
	$(CPP) $(CPPFLAGS) -o $@ $^ $(INC) $(LIBS)

clean:
	rm -rf *.o **/*.o $(BINS) *.dSYM callgrind.out.*

.PHONY: all clean callgrind valgrind
