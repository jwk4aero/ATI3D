##########################################
#######	        Makefile GRID       #######
##########################################

SD = .

PROGRAM = SOLVER
#PROGRAM = GRID

#include $(SD)/MAKE_default
include $(SD)/MAKE_bechlars



BUILDDIR = $(SD)/$(PROGRAM)
BUILDOBJ = $(SD)/objects_$(PROGRAM)
SRCGRID = $(SD)/grid
SRCSOL = $(SD)/solver_files
SRCCOMMON = $(SD)/common
SRCLIB = $(SD)/lib

    
#	A target must be supplied.

no-target:
	@echo "Please supply a target (e.g. make jet)"

folders:
	mkdir -p $(BUILDDIR)
	mkdir -p $(BUILDOBJ)
	cp Makefile $(BUILDDIR)/.
	sed -i 's/SD = ./SD = ../g' $(BUILDDIR)/Makefile

rm-folders:
	rm -rf $(BUILDDIR)
	rm -rf $(BUILDOBJ)

#
#	Southampton 'iridis2' Pentium Beowulf cluster using MPICH.
#
#=========================================================

#	hector - common 

#=========================================================

hector-cr-omp:
	make    $(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -s real64 -O3 -Oaggress -Omsgs " \
		"LFLAGS= -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
		"LIBS=/usr/local/packages/nag/tecio/XT/12.0/tecio.a -lstdc++ -lfftw3 -lm" \
		"LIBS2="

hector-cr-omp-debug:
	make    $(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -Rbcps -ecCD -s real64 -K trap=fp -g -O0 -Omsgs" \
		"LFLAGS= -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
		"LIBS= -lstdc++ -lfftw3 -lm" \
		"LIBS2="


hector-cr:
	make    $(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -O3 -h noomp  " \
		"LFLAGS= -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
		"LIBS= -lstdc++ -lfftw3 -lm" \
		"LIBS2="

hector-cr-debug:
	make    $(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -Rbcps -ecCD -s real64 -K trap=fp -g -O0 -Omsgs -h noomp " \
		"LFLAGS= -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
		"LIBS= -lstdc++ -lfftw3 -lm" \
		"LIBS2="
#=========================================================

#	iridis - common (intel compiler - intelmpi)

#=========================================================

iridis-intel-omp:
	make 	$(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS) -DIRIDIS" \
		"OPTIONS2= -DCRAY -D_HPCX_" \
		"CPP=/lib/cpp" \
		"CPPFLAGS=-P -traditional " \
		"FC=mpiifort"\
		"FFLAGS= -r8 -mt_mpi -O3 -heap-arrays" \
		"LFLAGS=-O3" \
		"LIBS=-L/lib -lfftw3" \
		"LIBS2=-L$(HOME)/lib -lstdc++ -ltecio -lm "


iridis-intel:
	make 	$(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS) -DIRIDIS" \
		"OPTIONS2= -DCRAY -D_HPCX_" \
		"CPP=/lib/cpp" \
		"CPPFLAGS=-P -traditional " \
		"FC=mpiifort"\
		"FFLAGS= -r8 -O3 -heap-arrays" \
		"LFLAGS=-O3" \
		"LIBS=-L/lib -mkl -I$(FFTW_INC_DIR) " \
		"LIBS2=-L$(HOME)/lib -lstdc++ -lm "

iridis-intel-deb:
	make 	$(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS) -DIRIDIS" \
		"OPTIONS2= -DCRAY -D_HPCX_" \
		"CPP=/lib/cpp" \
		"CPPFLAGS=-P -traditional " \
		"FC=mpiifort"\
		"FFLAGS= -g -C -fpe0 -traceback -r8 -check all -heap-arrays" \
		"LFLAGS=-O0" \
		"LIBS=-L/lib -lfftw3 -mkl" \
		"LIBS2=-L$(HOME)/lib -lstdc++ -ltecio -lm "

#=========================================================

#	local - common (gnu compiler - openmpi)

#=========================================================
		
test:
	make	 $(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS) " \
		"OPTIONS2= -DCRAY -D_HPCX_" \
		"CPP=/lib/cpp" \
		"CPPFLAGS=-P -traditional -I/usr/include" \
		"FC=mpif90 -fcray-pointer"\
		"FFLAGS= -O3" \
		"LFLAGS=-O3 " \
		"LIBS=-L/lib -llapack" \
		"LIBS2=-L$(HOME)/lib -lstdc++ -lfftw3 -lm "

test-debug:
	make    $(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS) -Dlocal" \
		"OPTIONS2= -DCRAY -D_HPCX_" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=mpif90 -fcray-pointer"\
		"FFLAGS= -g -fbacktrace -O0 -fbounds-check -Wconversion" \
		"LFLAGS=-O0" \
		"LIBS=-L/lib -llapack" \
		"LIBS2=-L$(HOME)/lib -lstdc++ -lfftw3 -lm "
#=========================================================

#	user defined, individual make targets

#=========================================================

hector-pgi:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -O3 -fastsse -r8  -Minfo=mp" \
		"LFLAGS= -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
		"LIBS=/usr/local/packages/nag/tecio/XT/12.0/tecio.a -lstdc++ -lfftw3 -lm" \
		"LIBS2="

iridis:
	make 	$(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS) -DIRIDIS" \
		"OPTIONS2= -DCRAY -D_HPCX_" \
		"CPP=/lib/cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=mpif90"\
		"FFLAGS= -r8 -O3" \
		"LFLAGS=-O3" \
		"LIBS=-L/lib -lfftw3" \
		"LIBS2=-L$(HOME)/lib -lstdc++ -ltecio -lm "

iridis-gnu:
	make 	$(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS) -DIRIDIS" \
		"OPTIONS2= -DCRAY -D_HPCX_" \
		"CPP=/lib/cpp" \
		"CPPFLAGS=-P -traditional -I/usr/include" \
		"FC=mpif90 -fcray-pointer"\
		"FFLAGS= -O3" \
		"LFLAGS=-O3" \
		"LIBS=-L/lib -lfftw3" \
		"LIBS2=-L$(HOME)/lib -lstdc++ -ltecio -lm "

iridis-gnu-deb:
	make 	$(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS) -DIRIDIS" \
		"OPTIONS2= -DCRAY -D_HPCX_" \
		"CPP=/lib/cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=mpif90 -fcray-pointer"\
		"FFLAGS= -g -fbacktrace -O0 -fbounds-check -Wconversion" \
		"LFLAGS=-O0" 
		"LIBS=" \
		"LIBS2=-L$(HOME)/lib -lstdc++ -ltecio -lm "

iridis-debug:
	make	$(PROGRAM) "EXE=$(EXE).x" \
		"OPTIONS=$(OPTIONS) -DIRIDIS" \
		"OPTIONS2= -DCRAY -D_HPCX_" \
		"CPP=/lib/cpp" \
		"CPPFLAGS=-P -traditional -I/local/mpich-pgi/include" \
		"FC=mpif90"\
		"FFLAGS= -g -check all -traceback -check bounds" \
		"LFLAGS= " \
		"LIBS=" \
		"LIBS2=-L$(HOME)/lib -lstdc++ -ltecio -lm "



hector-gnu:
	make    $(PROGRAM) "EXE=$(EXE)" \
      	"OPTIONS=$(OPTIONS)" \
      	"OPTIONS2=-DCRAY -D_HPCX_ -I/opt/cray/mpt/5.3.4/xt/gemini/mpich2-gnu/46/include/" \
      	"CPP=cpp" \
      	"CPPFLAGS=-P -traditional" \
      	"FC=ftn" \
      	"FFLAGS=-O3 -fcray-pointer" \
      	"LFLAGS=" \
		"LIBS=/usr/local/packages/nag/tecio/XT/12.0/tecio.a -lstdc++ -lfftw3 -lm" \
		"LIBS2="

hector-gnu-omp:
	make    $(PROGRAM) "EXE=$(EXE)" \
      	"OPTIONS=$(OPTIONS)" \
      	"OPTIONS2=-DCRAY -D_HPCX_ -I/opt/cray/mpt/5.3.4/xt/gemini/mpich2-gnu/46/include/" \
      	"CPP=cpp" \
      	"CPPFLAGS=-P -traditional" \
      	"FC=ftn" \
      	"FFLAGS=-O3 -fcray-pointer" \
      	"LFLAGS=-fopenmp" \
		"LIBS=/usr/local/packages/nag/tecio/XT/12.0/tecio.a -lstdc++ -lfftw3 -lm" \
		"LIBS2="


hector-pgi-debug:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -g -Mbounds -O1 -fastsse -r8  -Minfo=mp" \
		"LFLAGS= -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
		"LIBS=/usr/local/packages/nag/tecio/XT/12.0/tecio.a -lstdc++ -lfftw3 -lm" \
		"LIBS2="
hector-hybrid-pgi:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -O3 -fastsse -r8 -mp=nonuma -Minfo=mp" \
		"LFLAGS=-mp=nonuma -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
		"LIBS=/usr/local/packages/nag/tecio/XT/12.0/tecio.a -lstdc++ -lfftw3 -lm" \
		"LIBS2="

hector-pgi-deb:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -O3 -fastsse -r8  -Minfo=mp -g " \
		"LFLAGS= -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
		"LIBS=/usr/local/packages/nag/tecio/XT/12.0/tecio.a -lstdc++ -lfftw3 -lm" \
		"LIBS2="

test-hector-debug:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS=-Rbcps -ecCD -s real64 -K trap=fp" \
                "LFLAGS=-ecCD -s real64 -K trap=fp -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
                "LIBS=-L$(HOME)/lib -ltecio -lstdc++ -lfftw3 -lm" \
                "LIBS2="

test-hector:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS=-s real64 -O3" \
                "LFLAGS=-s real64 -O3 -I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
                "LIBS=-L$(HOME)/lib -ltecio -lstdc++ -lfftw3 -lm" \
                "LIBS2="

#	"FFLAGS= -g -C -fcray-pointer -ffpe-trap=invalid -fbacktrace" \

test-hector-gfortran:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -fcray-pointer -O2" \
                "LFLAGS=-I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
                "LIBS=-L$(HOME)/lib -ltecio -lstdc++ -lfftw3 -lm" \
                "LIBS2="

#	"FFLAGS= -g -fpe0 -traceback -r8 " \

test-iridis:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_ " \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=mpiifort"\
		"FFLAGS=" \
		"LFLAGS=" \
		"LIBS=-L$(HOME)/lib -L/local/software/rh53/fftw/3.2.2/intel/double/lib -ltecio -lstdc++ -lfftw3 -lm " \
		"LIBS2="

test-iridis-debug:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_ " \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=mpiifort"\
		"FFLAGS= -g -C -fpe0 -traceback -r8 " \
		"LFLAGS=" \
		"LIBS=-L$(HOME)/lib -L/local/software/rh53/fftw/3.2.2/intel/double/lib -ltecio -lstdc++ -lfftw3 -lm" \
		"LIBS2="

test-hector-debug-pgi:
	make    $(PROGRAM) "EXE=$(EXE)" \
		"OPTIONS=$(OPTIONS)" \
		"OPTIONS2=-DCRAY -D_HPCX_  -I/opt/xt-mpt/default/mpich2-64/P/include/" \
		"CPP=cpp" \
		"CPPFLAGS=-P -traditional" \
		"FC=ftn"\
		"FFLAGS= -g -Mbounds -Mchkstk -traceback -r8 " \
		"LFLAGS=-I/usr/local/packages/nag/tecio/XT/12.0/tecsrc" \
		"LIBS=-L$(HOME)/lib -ltecio -lstdc++ -lfftw3 -lm" \
		"LIBS2="
		



# Load the correct libraries depending on the pre-compiler flags
# cpb
ifneq (,$(findstring -Dlapack ,$(OPTIONS)))
LIBS3 :=  -llapack 
endif	



OBJECTS_SOLVER = $(BUILDOBJ)/stdtypes.o \
	$(BUILDOBJ)/Main3D.o \
	$(BUILDOBJ)/MainVar3D.o	\
	$(BUILDOBJ)/Subroutineso.o	\
	$(BUILDOBJ)/Subroutines3D.o	\
	$(BUILDOBJ)/AerofoilTurbulence3D.o	\
	$(BUILDOBJ)/GridAerofoil.o	\
	$(BUILDOBJ)/Subspace.o	\
	$(BUILDOBJ)/SubspaceVar.o  \
	$(BUILDOBJ)/stopnow.o	

.SUFFIXES:	.F90 .f90 .o

$(BUILDDIR)/%.f90 :  $(SRCLIB)/%.F90
	$(CPP) $(CPPFLAGS) $(OPTIONS) $(OPTIONS2) -Wno-endif-labels $(SRCLIB)/$*.F90 $(BUILDDIR)/$*.f90
$(BUILDDIR)/%.f90 :  $(SRCCOMMON)/%.F90
	$(CPP) $(CPPFLAGS) $(OPTIONS) $(OPTIONS2) -Wno-endif-labels $(SRCCOMMON)/$*.F90 $(BUILDDIR)/$*.f90 
$(BUILDDIR)/%.f90 :  $(SRCGRID)/%.F90
	$(CPP) $(CPPFLAGS) $(OPTIONS) $(OPTIONS2) -Wno-endif-labels $(SRCGRID)/$*.F90 $(BUILDDIR)/$*.f90 
$(BUILDDIR)/%.f90 :  $(SRCSOL)/%.F90
	$(CPP) $(CPPFLAGS) $(OPTIONS) $(OPTIONS2) -Wno-endif-labels $(SRCSOL)/$*.F90 $(BUILDDIR)/$*.f90 

$(BUILDOBJ)/%.o :  $(BUILDDIR)/%.f90
	$(FC) $(FFLAGS) -c -o $(BUILDOBJ)/$*.o $(BUILDDIR)/$*.f90 -J$(BUILDOBJ)

GRID:	$(OBJECTS_GRID)
	$(FC) $(LFLAGS)  -o $(EXE) $(OBJECTS_GRID) $(LIBS) $(LIBS2) $(LIBS3) 

SOLVER:	$(OBJECTS_SOLVER)
	$(FC) $(LFLAGS)  -o $(EXE) $(OBJECTS_SOLVER) $(LIBS) $(LIBS2) $(LIBS3) 

clean:
	rm -f $(BUILDDIR)/*.f90 $(BUILDOBJ)/*.o $(BUILDOBJ)/*.mod

clean-complete:
	rm -rf $(BUILDDIR)/ $(BUILDOBJ)/

clean-objects:
	rm -f $(BUILDOBJ)/*.o

$(BUILDDIR)/stdtypes.f90:		$(SRCLIB)/stdtypes.F90
$(BUILDDIR)/Main3D.f90:		$(SRCSOL)/Main3D.F90
$(BUILDDIR)/MainVar3D.f90:		$(SRCSOL)/MainVar3D.F90
$(BUILDDIR)/Subroutineso.f90:		$(SRCSOL)/Subroutineso.F90
$(BUILDDIR)/Subroutines3D.f90:		$(SRCSOL)/Subroutines3D.F90
$(BUILDDIR)/AerofoilTurbulence3D.f90:		$(SRCSOL)/AerofoilTurbulence3D.F90
$(BUILDDIR)/GridAerofoil.f90:		$(SRCSOL)/GridAerofoil.F90
$(BUILDDIR)/Subspace.f90:		$(SRCSOL)/Subspace.F90
$(BUILDDIR)/SubspaceVar.f90:		$(SRCSOL)/SubspaceVar.F90
$(BUILDDIR)/stopnow.f90:		$(SRCLIB)/stopnow.F90

$(BUILDOBJ)/stdtypes.o: $(BUILDDIR)/stdtypes.f90	
$(BUILDOBJ)/Main3D.o: $(BUILDDIR)/Main3D.f90	$(BUILDOBJ)/Subroutineso.o $(BUILDOBJ)/Subroutines3D.o $(BUILDOBJ)/AerofoilTurbulence3D.o $(BUILDOBJ)/SubspaceVar.o
$(BUILDOBJ)/MainVar3D.o: $(BUILDDIR)/MainVar3D.f90	
$(BUILDOBJ)/Subroutineso.o: $(BUILDDIR)/Subroutineso.f90	$(BUILDOBJ)/MainVar3D.o
$(BUILDOBJ)/Subroutines3D.o: $(BUILDDIR)/Subroutines3D.f90	$(BUILDOBJ)/stdtypes.o $(BUILDOBJ)/Subroutineso.o $(BUILDOBJ)/MainVar3D.o
$(BUILDOBJ)/AerofoilTurbulence3D.o: $(BUILDDIR)/Subroutines3D.f90	$(BUILDOBJ)/Subroutineso.o $(BUILDOBJ)/GridAerofoil.o
$(BUILDOBJ)/GridAerofoil.o: $(BUILDDIR)/GridAerofoil.f90	$(BUILDOBJ)/Subroutineso.o
$(BUILDOBJ)/SubspaceVar.o:		$(BUILDDIR)/SubspaceVar.f90 $(BUILDOBJ)/SubspaceVar.o
$(BUILDOBJ)/Subspace.o:		$(BUILDDIR)/Subspace.f90 $(BUILDOBJ)/Subspace.o $(BUILDOBJ)/SubspaceVar.o
$(BUILDOBJ)/stopnow.o: $(BUILDDIR)/stopnow.f90	

#
# Directive file for the Cray linker to preset all memory to NaN
#
debug.cld:
	echo "preset=nan" > debug.cld
