######################################################################
# Makefile for disorder
#
#  - "make" creates the vbfh executable
#  - "make clean" removes all files created by the "make" command
######################################################################

# include the settings from Makefile.inc, generated by configure script
include Makefile.inc

# main program and modules to be compiled
MAIN = disorder
MODULES = integration io_utils lcl_dec mod_parameters mod_phase_space mod_matrix_element mod_dsigma disent-lib mod_analysis
ANALYSIS = fastjetfortran pwhg_bookhist-multi 
#ANALYSIS += sigmaR
ANALYSIS += cut_Ecur 
#ANALYSIS += jets_lab_frame

FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STD)
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)

#c++ compiler
CXX = g++
CXXFLAGS = -O3	
CXXFLAGS += $(shell $(FASTJET_CONFIG) --cxxflags)
#CXXLIBS= -lstdc++ -L$(PWD) -Wl,-rpath=$(PWD) -lpanscales
CXXLIBS= -lstdc++

# librairies and flags needed for compilation
#FF = ifort
AUX=$(PWD)/aux
SRC=$(PWD)/src
OBJ=$(PWD)/obj
ANA=$(PWD)/analysis
VPATH=./:$(SRC):/$(OBJ):/$(ANA)

FFLAGS= -O3 -ffixed-line-length-132 #-Wunused
FFLAGS+= $(shell $(HPEXEC) --fflags)
INCLUDE= -I$(ANA) -I$(SRC) $(wildcard *.h)
FFLAGS+= $(INCLUDE) -J$(OBJ)
LDFLAGS= $(shell $(HPEXEC) --libs) $(shell $(LHEXEC) --libs) 

all: disorder mergegrids mergedata 

# main program
$(MAIN): %: %.o $(addsuffix .o,$(MODULES))  $(addsuffix .o,$(ANALYSIS)) Makefile
	$(FF) $(OBJ)/$@.o $(patsubst %,$(OBJ)/%,$(addsuffix .o,$(MODULES))) \
	$(patsubst %,$(OBJ)/%,$(addsuffix .o,$(ANALYSIS))) $(CXXLIBS) $(FFLAGS) $(LDFLAGS) $(LIBSFASTJET) -o $@

mergegrids:
	$(FF) $(FFLAGS) -mcmodel=large -o $(AUX)/$@ $(AUX)/mergegrids.f

mergedata:
	$(FF) $(FFLAGS) -mcmodel=large -o $(AUX)/$@ $(AUX)/mergedata.f

# object files
%.o: %.f Makefile 
	$(FF) $(FFLAGS) -o $(OBJ)/$@ $< -c

%.o: %.f90 Makefile
	$(FF) $(FFLAGS) -o $(OBJ)/$@ $< -c

%.o: %.cc Makefile
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $(OBJ)/$@ $<


# f90 module dependencies
mod_dsigma.o: mod_parameters.o mod_phase_space.o mod_matrix_element.o mod_analysis.o
mod_phase_space.o: mod_parameters.o 
disorder.o: mod_matrix_element.o mod_parameters.o mod_phase_space.o integration.o analysis.o mod_dsigma.o
mod_parameters.o: io_utils.o lcl_dec.o integration.o
mod_matrix_element.o: mod_parameters.o 
#histo.o: types.o consts.o
analysis.o: mod_parameters.o
mod_analysis.o: mod_parameters.o
pwhg_bookhist-multi.o: mod_parameters.o
$(addsuffix .o,$(ANALYSIS)): mod_parameters.o

# make clean
clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod *~ *.log fort* $(MAIN)

