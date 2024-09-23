TARGET = cmsa cplex greedy

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = /home/cblum/ILOG/CPLEX_Studio221/cplex
CONCERTDIR    = /home/cblum/ILOG/CPLEX_Studio221/concert
GCC = gcc
CCC = g++
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++11 -fpermissive -w
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread -lpthread -ldl
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread -lpthread -ldl
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

all: ${TARGET}

cmsa: cmsa.o $(OBJS)
	$(CCC) $(CCFLAGS) cmsa.o -o cmsa $(CCLNFLAGS)

cmsa.o: cmsa.cpp
	$(CCC) -c $(CCFLAGS) cmsa.cpp -o cmsa.o 

cplex: cplex.o $(OBJS)
	$(CCC) $(CCFLAGS) cplex.o -o cplex $(CCLNFLAGS)

cplex.o: cplex.cpp
	$(CCC) -c $(CCFLAGS) cplex.cpp -o cplex.o 

mds_cplex: mds_cplex.o $(OBJS)
	$(CCC) $(CCFLAGS) mds_cplex.o -o mds_cplex $(CCLNFLAGS)

mds_cplex.o: mds_cplex.cpp
	$(CCC) -c $(CCFLAGS) mds_cplex.cpp -o mds_cplex.o 

greedy: greedy.o $(OBJS)
	$(CCC) $(CCFLAGS) greedy.o -o greedy $(CCLNFLAGS)

greedy.o: greedy.cpp
	$(CCC) -c $(CCFLAGS) greedy.cpp -o greedy.o 
	
clean:
	@rm -f *~ *.o ${TARGET} core

