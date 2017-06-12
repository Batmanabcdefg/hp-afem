CPP=g++
CPPFLAGS=-std=c++11 -Wall -ggdb3
INC= -lginac -lcln -I /usr/local/include/eigen3 -I /usr/include/eigen3
SRCS=
OBJS=$(subst .cpp,.o,$(SRCS))

.PHONY: clean

default: createbases

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<

createbezierbases: createbezierbases.o $(OBJS)
	$(CPP) $(INC) $(CPPFLAGS) -o createbezierbases createbezierbases.o $(OBJS)

createszabobases: createszabobases.o $(OBJS)
	$(CPP) $(INC) $(CPPFLAGS) -o createszabobases createszabobases.o $(OBJS)

createbases: createbases.o $(OBJS)
	$(CPP) $(INC) $(CPPFLAGS) -o createbases createbases.o $(OBJS)

clean:
	rm -rf *.o createbases createszabobases createbezierbases *.dSYM callgrind.out.*
