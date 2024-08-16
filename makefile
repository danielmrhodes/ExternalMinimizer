CURDIR = $(shell pwd)
SRCDIR = $(CURDIR)/src
INCDIR = $(CURDIR)/include
BINDIR = $(CURDIR)/bin

ROOT_LIBS = `root-config --glibs` -lSpectrum -lTreePlayer -lMathMore

LIBRS = $(ROOT_LIBS)
INCLUDE = $(INCDIR)

CPP        = g++
CFLAGS = -std=c++11 -g -fPIC `root-config --cflags` `gsl-config --cflags`

HEAD = $(wildcard include/*.h)
OBJECTS = $(patsubst include/%.h,lib/%.so,$(HEAD))

$(info   OBJECTS = $(OBJECTS))

TARGET = bin/libExternalMinimizer.so

main: $(TARGET)
	@printf "Make complete\n"

$(TARGET): $(OBJECTS) bin/DictOutput.cxx lib
	@printf "Now compiling shared library $@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -I. $(LIBRS) -o $@ -shared bin/DictOutput.cxx $(OBJECTS)

bin/DictOutput.cxx: $(HEAD)
	@printf "Linking libraries\n"
	@rootcint -f $@ $(HEAD) lib/linkdef.h

bin:
	@mkdir bin

obj:
	@mkdir obj

lib/%.so: src/%.cc include/%.h 
	@printf "Now compiling library $@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -L$(LIBRS) -o $@ -shared -c $< 

clean:  
	@printf "Tidying up...\n"
	@rm $(OBJECTS)
