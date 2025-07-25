LIBS = -lExternalMinimizer `root-config --glibs` -lSpectrum -lTreePlayer -lMathMore
CPP = g++
CFLAGS = -std=c++11 -g -fPIC `root-config --cflags` `gsl-config --cflags`
INCDIR = ../../include
LIBDIR = ../../lib
SRCDIR = ../../bin
CURDIR = $(shell pwd)
BINDIR = $(CURDIR)

SRCS := $(wildcard *.cxx)
EXES := $(patsubst %.cxx,$(BINDIR)/%,$(SRCS))
#OBJS := $(patsubst %,%.o,$(EXES))

all : $(EXES)

$(BINDIR)/%: $(CURDIR)/%.cxx
	@mkdir -p $(BINDIR)
	@printf "Creating executable $@ from $<\n"
	$(CPP) $< -o $@ $(CFLAGS) -I$(INCDIR) -L$(SRCDIR) -L$(LIBDIR) $(LIBS)

clean:
	@printf "Cleaning up\n"
	@rm $(EXES)
