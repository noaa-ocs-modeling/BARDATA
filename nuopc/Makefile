###########################################################################
## GNU Makefile for the BARDATA NUOPC Cap
##
##!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
##!>
##!> @date 03/07/2022
###########################################################################

SHELL := /bin/bash

clist := clean distclean
ifeq (,$(filter $(clist), $(MAKECMDGOALS)))
  ###########################################################################
  # Mandatory environment variables NETCDFHOME, ESMFMKFILE and possibly PAHM_BINDIR
  ifneq ($(origin NETCDFHOME), environment)
    $(error Environment variable NETCDFHOME was not set)
  endif
  #-----
  ifneq ($(origin ESMFMKFILE), environment)
    $(error Environment variable ESMFMKFILE was not set)
  endif
  #-----

  include $(ESMFMKFILE)

  INCDIRS += -I${NETCDFHOME}/include
  LIBDIRS += -L${NETCDFHOME}/lib

  ifneq ($(origin HDF5HOME), environment)
    INCDIRS += -I${HDF5HOME}/include   
    LIBDIRS += -L${HDF5HOME}/lib
  endif

  INCDIRS += -I.
  LIBDIRS += -L.

  # For reference only
  LIBS   :=  -lnetcdf  -lnetcdff -lhdf5 -lhdf5_fortran
  ################################################################################
endif

.SUFFIXES: .f90 .F90 .c .C

%.o : %.F90
	@echo ""
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) \
	$(ESMF_F90COMPILECPPFLAGS) $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<
	@echo ""

%.o : %.F
	@echo ""
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) \
	$(ESMF_F90COMPILECPPFLAGS) $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<
	@echo ""

.PRECIOUS: %.o

# -----------------------------------------------------------------------------
vars:
	@echo " >>>>>>>>>>>>>>>> vars >>>>>>>>>>>>>>"
	@echo "INCDIRS $(INCDIRS)"
	@echo "LIBDIRS $(LIBDIRS)"
	@echo "ESMF_F90COMPILER $(ESMF_F90COMPILER)"
	@echo "ESMF_F90LINKER $(ESMF_F90LINKER)"
	@echo "ESMF_F90LINKOPTS $(ESMF_F90LINKOPTS)"
	@echo "ESMF_F90LINKPATHS  $(ESMF_F90LINKPATHS)"
	@echo "ESMF_F90LINKRPATHS $(ESMF_F90LINKRPATHS)"
	@echo "ESMF_F90ESMFLINKLIBS $(ESMF_F90ESMFLINKLIBS)"	

# -----------------------------------------------------------------------------

PWDDIR := $(shell pwd)
INSTALLDIR := $(INSTDIR)    # coming from NEMS/src/incmake/component_BARDATA.mk

MAKEFILE = bardata.mk

LIBRARY  = libbardata_cap.a

MODULES  = bardata.o bardata_mod.o 

nuopc: $(LIBRARY) 

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $?

install: $(LIBRARY)
	@echo ""
	mkdir -p $(strip $(INSTALLDIR))
	@echo "# ESMF self-describing build dependency makefile fragment" > \
	      $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo                                                            >> \
	      $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_FRONT     = BARDATA    "                         >> \
	      $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_INCPATH   = $(PWDDIR) $(strip $(INSTALLDIR)) "   >> \
	      $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_CMPL_OBJS = "                                    >> \
	      $(strip $(INSTALLDIR))/$(MAKEFILE)
	@echo "ESMF_DEP_LINK_OBJS =  -L$(strip $(INSTALLDIR)) $(strip $(INSTALLDIR))/$(LIBRARY)  " >> \
	      $(strip $(INSTALLDIR))/$(MAKEFILE)
	cp -f *.a *.mod $(strip $(INSTALLDIR)) 
	@echo ""

# module dependencies:
bardata.o:   bardata_mod.o 

# -----------------------------------------------------------------------------
.PHONY: clean distclean

clean:
	@echo "Cleaning Model build..."
	@rm -v -f $(MODULES) *.mod $(LIBRARY)

distclean: clean
	@rm -v -f $(LIBRARY)
