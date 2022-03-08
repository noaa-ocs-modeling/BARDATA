###########################################################################
## GNU Makefile for OC-GTSM
##
##!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
##!>
##!> @date 03/07/2022
###########################################################################

SHELL := /bin/sh

clist := clean distclean
ifeq (,$(filter $(clist), $(MAKECMDGOALS)))
  include compiler_defs.mk
endif


SRCS :=  OGCM_DL.f90
	
OBJS :=   $(patsubst %.f90, %.o, $(SRCS) )

PROGRAM := ogcmdl


########################################
### BEG:: Compilation Rules
########################################

.SUFFIXES: .f90 .F90 .o

.f90.o:
	$(F90) -cpp -c $(FCFLAGS) $(INCDIRS) -o $@ $<

.F90.o:
	$(F90) -c $(FCFLAGS) $(INCDIRS) -o $@ $<

########################################
### END:: Compilation Rules
########################################


all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F90) $(OBJS) $(LDFLAGS) -o $@

### Module Dependencies:
OGCM_DL.o :  OGCM_DL.f90


########################################
### BEG:: Miscellaneous Rules
########################################
.PHONY: clean distclean

clean:
	@rm -f *.o *.mod

distclean: clean
	@rm -f $(PROGRAM)

########################################
### END:: Miscellaneous Rules
########################################

