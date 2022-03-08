###########################################################################
## GNU Makefile for OC-GTSM
##
##!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
##!>
##!> @date 03/07/2022
###########################################################################


###########################################################################
# Mandatory environment variables NETCDFHOME, HDF5HOME, ...

# For NetCDF
ifneq ($(origin NETCDFHOME), environment)
  $(warning Environment variable NETCDFHOME was not set, searching for NETCDF_HOME)
  ifeq ($(origin NETCDF_HOME), environment)
    NETCDFHOME := ${NETCDF_HOME}
  else
    $(error No suitable NETCDF environment variable found. Please set NETCDFHOME before running this Makefile)
  endif
endif

INCDIRS += -I${NETCDFHOME}/include
LIBDIRS += -L${NETCDFHOME}/lib
LIBS    += -lnetcdf -lnetcdff

# For HDF5
ifeq ($(origin HDF5HOME), environment)
  INCDIRS += -I${HDF5HOME}/include
  LIBDIRS += -L${HDF5HOME}/lib
  LIBS    += -lhdf5 -lhdf5_fortran
endif

# For datetime-fortran
ifneq ($(origin DATETIMEHOME), environment)
  DATETIMEHOME := ../datetime-fortran/build
endif
ifeq ("$(wildcard $(DATETIMEHOME))", "")
  $(error The essential datetime-fortran directory ${DATETIMEHOME} was not found.)
endif
INCDIRS += -I$(DATETIMEHOME)/include
LIBDIRS += -L$(DATETIMEHOME)/lib
LIBS    += -ldatetime

# For For GSW-Fortran
ifneq ($(origin DATETIMEHOME), environment)
  GSWHOME := ../GSW-Fortran/build
endif
ifeq ("$(wildcard $(GSWHOME))", "")
  $(error The essential GSW-Fortran directory ${GSWHOME} was not found.)
endif
INCDIRS += -I$(GSWHOME)/gsw
LIBDIRS += -L$(GSWHOME)
LIBS    += -lgsw
################################################################################

# We need a parallel compiler here
F90 ?=  mpif90

FCFLAGS ?= -g -O2 -ffree-line-length-none -fbacktrace #-g -check bounds 

LDFLAGS ?= -O2 $(LIBDIRS) $(LIBS)
