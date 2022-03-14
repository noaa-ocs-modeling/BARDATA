###########################################################################
## GNU Makefile for the BARDATA NUOPC Cap
##
##!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
##!>
##!> @date 03/07/2022
###########################################################################

SHELL := /bin/bash

PWDDIR := $(shell pwd)

clist := clean distclean
ifeq (,$(filter $(clist), $(MAKECMDGOALS)))
  include compiler_defs.mk
endif

# Coming from NEMS/src/incmake/component_BARDATA.mk (option for this Makefile)
INSTALLDIR := $(INSTDIR)

MAKEFILE   := bardata.mk

LIBRARY    := libbardata_cap.a

MODULES    := bardata_cap.o bardata_mod.o 


########################################
### BEG:: Compilation Rules
########################################

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

########################################
### END:: Compilation Rules
########################################


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

### Module Dependencies:
bardata_cap.o:   bardata_mod.o 


########################################
### BEG:: Miscellaneous Rules
########################################

.PHONY: clean distclean

clean:
	@echo "Cleaning Model build..."
	@rm -v -f $(MODULES) *.mod $(LIBRARY)

distclean: clean
	@rm -v -f $(LIBRARY)

########################################
### END:: Miscellaneous Rules
########################################
