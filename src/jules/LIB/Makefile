###############################################################################
##
## This is the top-level makefile for building a JULES executable.
##
## This makefile is maintained to be compliant with the GNU make utility ONLY, 
## although it will most probably work with other make utilities.  The GNU make
## utility is available free-of-charge from http://www.gnu.org/software/make/, 
## but most Unix/Linux systems have it lurking somewhere anyway.  (Type
## 'gmake -v' or 'which gmake' to check.)
##
## To build, you will need to set the following three variables to point to 
## the directory in which you are building JULES (most probably the one 
## containing this makefile) and the location of the pre-compiled netCDF 
## Fortran interface.  The netCDF interface library is available free-of-charge
## from Unidata at http://www.unidata.ucar.edu/software/netcdf/.  Note that
## this makefile does NOT build the netCDF library, that must be done 
## beforehand as a separate process.
##
## CDF_LIB_PATH should point to the directory containing the netCDF library
## archive (e.g. libnetcdf.a).
## CDF_MOD_PATH should point to the directory containing the netCDF Fortran 90
## module file (e.g. netcdf.mod).
##
## There are two options for the build that can be passed to the makefile from
## the command line.  The variable COMPILER allows the user to tell the make 
## which compiler they are using and the makefile will pick up compiler-specific
## variables from the separate text files, e.g. Makefile.comp.sun, that come 
## with this distribution.  Currently, there are four allowed values for 
## COMPILER: sun, intel, g95, nag. 
##
## The variable BUILD allows the user to maintain several JULES libraries and 
## executables built with different compiler flags.  This is particularly useful
## when making any modifications to the JULES source code, as it allows the user
## to switch to a debugging version of the model without having to recompile
## all of the source code or edit any makefiles.  This makefile is set up to 
## build a "normal" version, a debugging version and an optimised/fast version.
## Currently, there are 3 allowed values for BUILD: run, debug, fast.
## 
## These two options can be specified from the command line, e.g.,
## 
##    make BUILD=run COMPILER=sun
##
## KNOWN BUGS.
## * Using Intel compiler on Linux, GNU make sometimes returns root "/" from 
##   $(PWD) command.  Get around this by specifying JULESDIR=/home/jules on 
##   the command line.
## * Using Intel compiler on Linux, GNU make does not move the .mod files to 
##   MOD_DIR after each sub-make completes.  Get around this by adding 
##   -module $(JULESDIR)/MODS to FF_* variables in Makefile.comp.intel.
##
###############################################################################
JULESDIR=$(PWD)


.SUFFIXES: .C .f90 .F90 .c .o


#package        = ../src/ramsmain.f90
version        = 5.0.1
PREFIX         = /usr/local
FC             = @FC@

CHEM	       = RELACS_TUV

CCOMP	       = /opt/pgi/osx86-64/12.6/bin/pgcc
CPCOMP	       = /opt/hdf5-pgi12.6/bin/h5pcc
FCOMP	       = /opt/pgi/osx86-64/12.6/bin/pgf90
FPCOMP	       = /opt/hdf5-pgi12.6/bin/h5pfc


## The user is advised not to modify code below this line. ##


######################################
## Make some decisions about the    ##
## netCDF library.                  ##
######################################
CDFDUMMY=true
OBJ_CDF_DUMMY=
ifeq ($(CDFDUMMY),true)
DIR_CDF_DUMMY=$(JULESDIR)/UTILS/netcdf_dummy
CDF_MOD_PATH=$(DIR_CDF_DUMMY)
CDF_LIB_PATH=$(DIR_CDF_DUMMY)
OBJ_CDF_DUMMY=$(DIR_CDF_DUMMY)/jules_netcdf_dummy.o
CLEAN_DEPENDENCIES += clean_cdfdummy
endif

ifndef CDF_LIB_PATH
$(error You need to set the variable CDF_LIB_PATH either in the JULES Makefile or from the command line)
endif

ifndef CDF_MOD_PATH
$(error You need to set the variable CDF_MOD_PATH either in the JULES Makefile or from the command line)
endif

#######################################
## Set up variables for drhook dummy ##
#######################################
DIR_DRHOOK_DUMMY=$(JULESDIR)/UTILS/drhook_dummy
DRHOOK_OBJECTS=$(DIR_DRHOOK_DUMMY)/parkind1.o \
               $(DIR_DRHOOK_DUMMY)/yomhook.o
CLEAN_DEPENDENCIES += clean_drhook_dummy

######################################
## Rules for making library archive ##
## members from source files.       ##
######################################
include $(JULESDIR)/Makefile.common.mk

######################################
## Choose compiler.                 ##
######################################
## Set default values for these variables else you can get 
## some nasty "rm -f *" commands issued during a clean.


ifndef FPCOMP
$(error You need to set the variable FPCOMP either in configure command line)
endif

ifeq ($(notdir $(FPCOMP)),mpif90)
F_COMP=$(FPCOMP) -f90=$(FCOMP)
C_COMP=$(CPCOMP) -cc=$(CCOMP)	
endif

ifeq ($(notdir $(FPCOMP)),ftn)
F_COMP=$(FPCOMP) 
C_COMP=$(CPCOMP) 
endif

ifeq ($(notdir $(FPCOMP)),h5pfc)
F_COMP=$(FPCOMP) 
C_COMP=$(CPCOMP) 
endif


ifeq ($(notdir $(FCOMP)),pgf90)
L_OPTS=-fast -O0 -fastsse -Munroll=n:4 -Mipa=fast,inline,safe  -Msmart
endif 

ifeq ($(notdir $(FCOMP)),ifort)
L_OPTS=-O0
endif 

ifeq ($(notdir $(FCOMP)),gfortran)
L_OPTS=-O0
endif 


F_OPTS=$(L_OPTS) -I../src/jules/LIB/MODS -I../src/jules/LIB/UTILS/netcdf_dummy/ -Wl,--no-as-needed -ldl

C_OPTS=$(F_OPTS)

LOADER_OPTS=$(F_OPTS)



LIB_INC=-L
LIB_PRE=-l
LIB_FPRE=lib
LIB_FSUF=.a
MOD_INC=-M 
MOD_FSUF=.mod

COMPILER=$(notdir $(FCOMP))

ifndef COMPILER
$(error You need to set the variable COMPILER either in the JULES Makefile or from the command line)
endif

ifeq ($(COMPILER),sun)
include $(JULESDIR)/Makefile.comp.sun
endif

ifeq ($(COMPILER),g95)
include $(JULESDIR)/Makefile.comp.g95
endif

ifeq ($(COMPILER),gfortran)
include $(JULESDIR)/Makefile.comp.gfortran.x86
endif

ifeq ($(COMPILER),nag)
include $(JULESDIR)/Makefile.comp.nag
endif

ifeq ($(COMPILER),ifort)
include $(JULESDIR)/Makefile.comp.intel
endif

ifeq ($(COMPILER),pgf90)
include $(JULESDIR)/Makefile.comp.pgf
endif

ifeq ($(COMPILER),misc)
include $(JULESDIR)/Makefile.comp.misc
endif

ifeq ($(COMPILER)$(FC),misc@@)
$(error You can only use COMPILER misc if you have set the variables in Makefile.comp.misc)
endif


FC=$(F_COMP)
CC=$(C_COMP)
LINKER=$(F_COMP)

######################################
## Construct an number of useful    ##
## variables based on user/compiler ##
## options.                         ##
######################################
EXENAME=jules$(LIB_SPECIAL).exe
JULESLIBNOM=$(LIB_FPRE)jules$(LIB_SPECIAL)$(LIB_FSUF)
JULESLIB=$(JULESDIR)/$(JULESLIBNOM)
ARC=$(JULESLIB)

MAINOBJNOM=jules.o
MAINOBJ=$(JULESDIR)/$(MAINOBJNOM)
HEADMAKE=$(JULESDIR)/Makefile $(JULESDIR)/Makefile.common.mk
MOD_DIR=$(JULESDIR)/MODS
PAR_DIR=$(JULESDIR)/SOURCE/IMOGEN/MODULES/PARAMS
PAR_DIR2=$(JULESDIR)/SOURCE/MODULES/PARAMS
PAR_VAR=$(JULESDIR)//SOURCE/IMOGEN/MODULES/VARS/
PAR_CNT=$(JULESDIR)/SOURCE/MODULES/CONTROL/
PAR_SNT=$(JULESDIR)/SOURCE/SUBROUTINES/CONTROL/
PAR_MAP=$(JULESDIR)/SOURCE/IMOGEN/MODULES/MAP/
INCLUDE_DIR=$(JULESDIR)/INCLUDES

LIB_PATH  = $(LIB_INC)$(JULESDIR) 
LIB_PATH += $(LIB_INC)$(CDF_LIB_PATH) 

LIB_LIBS  = $(LIB_PRE)jules$(LIB_SPECIAL)
ifeq ($(CDFDUMMY),false)
LIB_LIBS += $(LIB_PRE)netcdf 
endif

MOD_PATH  = $(MOD_PUT)$(MOD_DIR)
MOD_PATH += $(MOD_INC)$(MOD_DIR)
MOD_PATH += $(MOD_INC)$(CDF_MOD_PATH)
MOD_PATH += $(MOD_INC)$(DIR_DRHOOK_DUMMY)
MOD_PATH += $(MOD_INC)$(PAR_DIR)
MOD_PATH += $(MOD_INC)$(PAR_DIR2)
MOD_PATH += $(MOD_INC)$(PAR_VAR)
MOD_PATH += $(MOD_INC)$(PAR_CNT)
MOD_PATH += $(MOD_INC)$(PAR_SNT)
MOD_PATH += $(MOD_INC)$(PAR_MAP)

CFLAGS=
FFLAGS=$(FF_RUN)
LDFLAGS=$(LIB_PATH) $(LIB_LIBS)

######################################
## Start of compilation/library     ##
## options.                         ##
######################################
# Options are run, debug, fast, condor.
BUILD=fast
LIB_SPECIAL=

ifeq ($(BUILD),debug)
  LIB_SPECIAL=_debug
  FFLAGS=$(FF_DBG)
  LINKER=$(FC)
endif

ifeq ($(BUILD),fast)
  LIB_SPECIAL=_fast
  FFLAGS=$(FF_FAST)
  LINKER=$(FC)
endif

ifeq ($(BUILD),condor)
  LIB_SPECIAL=_condor
  FFLAGS=$(FF_CON)
  LINKER=condor_compile $(FC)
endif

################################################
## Build pre-processor variable definitions   ##
################################################
FPP_VARS=L08_1A L19_1A L19_2A A71_1A SCMA BL_DIAG_HACK
FPP_DEFS=$(foreach FPP_VAR,$(FPP_VARS),$(FPP_FDEF)$(FPP_VAR))

################################################
## Set up the path for pre-processor includes ##
################################################
FPP_INC_PATH=$(FPP_INC)$(INCLUDE_DIR)

#######################################
## Make these variables available to ##
## the sub-makes.                    ##
#######################################
export FC
export FPP
export FPP_DEFS
export FPP_INC_PATH
export MOD_FSUF
export MOD_DIR
export CFLAGS
export FFLAGS
export MOD_PATH
export LIB_PATH
export LIB_LIBS
export JULESDIR
export ARC

######################################
## Variables for executing the      ##
## recursive build and clean.       ##
######################################
SRC_DIRS = SOURCE/MODULES/PARAMS \
           SOURCE/IMOGEN/MODULES/PARAMS \
           SOURCE/IMOGEN/MODULES/VARS \
           SOURCE/MODULES/CONTROL \
           SOURCE/IMOGEN/MODULES/MAP \
           SOURCE/SUBROUTINES/VEGETATION \
           SOURCE/SUBROUTINES/CONTROL \
           SOURCE/SUBROUTINES/INITIALISATION \
           SOURCE/SUBROUTINES/RADIATION \
           SOURCE/SUBROUTINES/SNOW \
           SOURCE/SUBROUTINES/SOIL \
           SOURCE/SUBROUTINES/SURFACE \
           SOURCE/IMOGEN/SUBROUTINES

CLEAN_SRC_DIRS=$(SRC_DIRS:%=clean_%)
MAKE_SRC_DIRS=$(SRC_DIRS:%=make_%)

CLEAN_DEPENDENCIES += $(CLEAN_SRC_DIRS)

######################################
## Dependencies for building library##
## and executable.                  ##
######################################

all : $(EXENAME)
.PHONY : all

$(EXENAME) : $(JULESLIB) $(MAINOBJ) $(OBJ_CDF_DUMMY) drhook $(HEADMAKE)
	$(LINKER) $(FFLAGS) -o $(EXENAME) $(MAINOBJ) \
	    $(OBJ_CDF_DUMMY) \
            $(DRHOOK_OBJECTS) \
	    $(LIB_PATH) \
	    $(MOD_PATH) \
	    $(LIB_LIBS)

$(JULESLIB) : $(OBJ_CDF_DUMMY) drhook $(MAKE_SRC_DIRS)
	touch $(JULESLIB)

$(MAINOBJ) : $(JULESLIB)
	$(EXTRACT) $(JULESLIB) $(MAINOBJNOM)

$(OBJ_CDF_DUMMY) : $(HEADMAKE)
	@$(MAKE) -C $(DIR_CDF_DUMMY)

drhook : $(HEADMAKE)
	@$(MAKE) -C $(DIR_DRHOOK_DUMMY)

$(MAKE_SRC_DIRS) : make_% : $(HEADMAKE)
	@$(MAKE) -C $*

######################################
## Dependencies for cleaning.       ##
######################################
.PHONY : clean clean_cdfdummy clean_drhook_dummy
clean: $(CLEAN_DEPENDENCIES)
	$(RM) $(JULESLIB) $(EXENAME) $(wildcard $(MOD_DIR)/*$(MOD_FSUF)) $(MAINOBJ)
$(CLEAN_SRC_DIRS) : clean_% :
	@$(MAKE) clean -C $*
clean_cdfdummy :
	@$(MAKE) clean -C $(DIR_CDF_DUMMY)
clean_drhook_dummy :
	@$(MAKE) clean -C $(DIR_DRHOOK_DUMMY)

## End of file.
