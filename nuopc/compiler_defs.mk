###########################################################################
## GNU Makefile for NUOPC
##
##!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
##!>
##!> @date 03/07/2022
###########################################################################


###########################################################################
# Mandatory environment variables NETCDFHOME, HDF5HOME, ...

############################
###### For NetCDF      #####
############################
ifeq ($(origin NETCDFHOME), environment)
  NETCDF_INCDIR := $(NETCDFHOME)/include
  NETCDF_LIBDIR := $(NETCDFHOME)/lib
else
  $(warning Environment variable NETCDFHOME was not set, searching for NETCDF, NETCDF_HOME or NETCDF_ROOT)
  ifeq ($(origin NETCDF), environment)
    NETCDFHOME := $(NETCDF)
    NETCDF_INCDIR := $(NETCDF)/include
    NETCDF_LIBDIR := $(NETCDF)/lib
  else
    ifeq ($(origin NETCDF_HOME), environment)
      NETCDFHOME := $(NETCDF_HOME)
      NETCDF_INCDIR := $(NETCDF_HOME)/include
      NETCDF_LIBDIR := $(NETCDF_HOME)/lib
    else
      ifeq ($(origin NETCDF_ROOT), environment)
        NETCDFHOME := $(NETCDF_ROOT)
        NETCDF_INCDIR := $(NETCDF_ROOT)/include
        NETCDF_LIBDIR := $(NETCDF_ROOT)/lib
      else
        NETCDFHOME:=$(shell nc-config --prefix 2>/dev/null)
        NETCDF_INCDIR:=$(shell nc-config --includedir 2>/dev/null)
        NETCDF_LIBDIR:=$(shell nc-config --libdir 2>/dev/null)
      endif
    endif
  endif
endif

ifeq ($(wildcard $(NETCDFHOME)),)
  $(error No suitable NETCDF environment variable found. Please set NETCDFHOME before running this Makefile)
else
  INCDIRS += -I$(NETCDFHOME)/include
  LIBDIRS += -L$(NETCDFHOME)/lib
  LIBS    += -lnetcdf -lnetcdff
endif


############################
###### For HDF5        #####
############################
ifeq ($(origin HDF5HOME), environment)
  HDF5_INCDIR := $(HDF5HOME)/include
  HDF5_LIBDIR := $(HDF5HOME)/lib
else
  ifeq ($(origin HDF5), environment)
    HDF5HOME := $(HDF5)
    HDF5_INCDIR := $(HDF5)/include
    HDF5_LIBDIR := $(HDF5)/lib
  else
    ifeq ($(origin HDF5_HOME), environment)
      HDF5HOME := $(HDF5_HOME)
      HDF5_INCDIR := $(HDF5_HOME)/include
      HDF5_LIBDIR := $(HDF5_HOME)/lib
    else
      ifeq ($(origin HDF5_ROOT), environment)
        HDF5HOME := $(HDF5_ROOT)
        HDF5_INCDIR := $(HDF5_ROOT)/include
        HDF5_LIBDIR := $(HDF5_ROOT)/lib
      endif
    endif
  endif
endif

ifeq ($(wildcard $(HDF5HOME)),)
  $(warning No suitable HDF5 environment variable found. If needed, please set HDF5HOME before running this Makefile)
else
  INCDIRS += -I$(HDF5HOME)/include
  LIBDIRS += -L$(HDF5HOME)/lib
  LIBS    += -lhdf5 -lhdf5_fortran
endif


############################
###### For ESMF        #####
############################
ifneq ($(origin ESMFMKFILE), environment)
  $(error No ESMFMKFILE environment variable found. Please set ESMFMKFILE before running this Makefile)
else
  include $(ESMFMKFILE)
endif

#INCDIRS += -I.
#LIBDIRS += -L.

################################################################################

