#Makefile include paths.mk

export

# RAMS root directory.
BUILD_DIR =$(shell pwd)
BIN_DIR   =$(BUILD_DIR)
RAMS_ROOT =$(BIN_DIR)/..

# Versions.
BRAMS_VERSION  = 5.2
JULES_VERSION  = 3.0
RAMS_VERSION   = 5.04
REVU_VERSION   = 2.3.1
UTILS_VERSION  = 2.0
HYPACT_VERSION = 1.2.0
POST_VERSION   = 2.0

# Source directories.#
# definition of aerosol source code:
ifeq ($(AER), MATRIX)
AER_DIR=matrix
AERLEVEL=$(AER)$(MATRIXLEV)
else
AER_DIR=ccatt
AERLEVEL=SIMPLE
endif

#__ (EXTERNAL LIBRARIES) _______

# JULES #
JULESID    =-jules-$(JULES_VERSION)
JULES_DIR  =$(RAMS_ROOT)/src/jules

JULES_LIB  =$(JULES_DIR)/LIB/libjules_fast.a \
            $(JULES_DIR)/LIB/UTILS/drhook_dummy/yomhook.o \
            $(JULES_DIR)/LIB/UTILS/netcdf_dummy/jules_netcdf_dummy.o
JULES_INC  =-I$(JULES_DIR)/LIB/MODS \
            -I$(JULES_DIR)/LIB/UTILS/drhook_dummy/ \
	    -I$(JULES_DIR)/LIB/UTILS/netcdf_dummy/
	    


# UTILS #
UTILS_INCS =$(RAMS_ROOT)/src/utils/include
UTILS_MODS =$(RAMS_ROOT)/src/utils/lib/modules
LIBUTILS   =$(BIN_DIR)/libutils-$(UTILS_VERSION).a

# POST #
POST_SRC      = $(RAMS_ROOT)/src/post
POST_INCS     = $(UTILS_INCS)
#
#
#IFS routines
IFS_DIR           =/Users/saulofreitas/work/ecmwf-closure/convect_ifs_40r1_tiedtke

UTILS_DUMP    =$(RAMS_ROOT)/src/utils/dump
UTILS_MDTOOLS    =$(RAMS_ROOT)/src/utils/model-devel-tools/src
#__ (BRAMS RELATED) ____________
# New paths like rams 6.x

BC            =$(RAMS_ROOT)/src/brams/bc
MODEL         =$(RAMS_ROOT)/src/brams/model
ADVC          =$(RAMS_ROOT)/src/brams/advect
CCATT         =$(RAMS_ROOT)/src/brams/ccatt
AEROSOL       =$(RAMS_ROOT)/src/brams/$(AER_DIR)
MODEL_CHEM    =$(CCATT)/$(CHEM)
ISAN_CHEM     =$(RAMS_ROOT)/src/brams/isan_chem
TUV           =$(CCATT)/TUV
CUPARM        =$(RAMS_ROOT)/src/brams/cuparm
FDDA          =$(RAMS_ROOT)/src/brams/fdda
INIT          =$(RAMS_ROOT)/src/brams/init
IO            =$(RAMS_ROOT)/src/brams/io
ISAN          =$(RAMS_ROOT)/src/brams/isan
MEMORY        =$(RAMS_ROOT)/src/brams/memory
MICRO         =$(RAMS_ROOT)/src/brams/micro
MKSFC         =$(RAMS_ROOT)/src/brams/mksfc
MPI           =$(RAMS_ROOT)/src/brams/mpi
NESTING       =$(RAMS_ROOT)/src/brams/nesting
RADIATE       =$(RAMS_ROOT)/src/brams/radiate
RRTMG_SW_SRC  =$(RAMS_ROOT)/src/brams/radiate/rrtmg_sw_v3.8/gcm_model/src
RRTMG_SW_MOD  =$(RAMS_ROOT)/src/brams/radiate/rrtmg_sw_v3.8/gcm_model/modules
RRTMG_LW_SRC  =$(RAMS_ROOT)/src/brams/radiate/rrtmg_lw_v4.85/gcm_model/src
RRTMG_LW_MOD  =$(RAMS_ROOT)/src/brams/radiate/rrtmg_lw_v4.85/gcm_model/modules
SIB           =$(RAMS_ROOT)/src/brams/sib
SOIL_MOISTURE =$(RAMS_ROOT)/src/brams/soil_moisture
SURFACE       =$(RAMS_ROOT)/src/brams/surface
TEB_SPM       =$(RAMS_ROOT)/src/brams/teb_spm
TURB          =$(RAMS_ROOT)/src/brams/turb
STILT         =$(RAMS_ROOT)/src/brams/stilt
MATRIX        =$(RAMS_ROOT)/src/brams/matrix
UTILS_LIB     =$(RAMS_ROOT)/src/utils/lib
ISAN          =$(RAMS_ROOT)/src/brams/isan
ISAN_MODS     =$(RAMS_ROOT)/src/brams/isan

## RMF - DEPRECIATED ROUTINES ###
#EFF           =$(RAMS_ROOT)/src/utils/eff
#NCARGD        =$(RAMS_ROOT)/src/utils/ncargd
#UTILS_DUMP    =$(RAMS_ROOT)/src/utils/plot
#################################

# executable & libraries names
#EXE     =$(BIN_DIR)/brams-$(BRAMS_VERSION)-CHEM_$(CHEM)$(JULESID)-AER_$(AER)
#LIBMODEL=$(BIN_DIR)/brams-$(BRAMS_VERSION)-CHEM_$(CHEM)$(JULESID)-AER_$(AER).a

EXE     =$(BIN_DIR)/brams-$(BRAMS_VERSION)
LIBMODEL=$(BIN_DIR)/brams-$(BRAMS_VERSION).a

## TERMINAL COLOURS PARAMETERS ##
NO_COLOUR   = "\033[0m"
BOLD        = "\033[1m"
RED         = "\033[0;31m"
BOLD_RED    = "\033[1;31m"
GREEN       = "\033[0;32m"
BOLD_GREEN  = "\033[1;32m"
YELLOW      = "\033[0;33m" 
BOLD_YELLOW = "\033[1;33m" 
BLUE        = "\033[0;34m" 
BOLD_BLUE   = "\033[1;34m" 
PURPLE      = "\033[0;35m" 
BOLD_PURPLE = "\033[1;35m" 
CYAN        = "\033[0;36m" 
BOLD_CYAN   = "\033[1;36m" 

