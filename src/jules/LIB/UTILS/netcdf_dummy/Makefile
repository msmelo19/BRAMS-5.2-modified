include ../../Makefile.common.mk
THISFILE = ./Makefile

OBJ_CDF_DUMMY=jules_netcdf_dummy.o

##############################
## Dependencies             ##
##############################
.PHONY : all clean
all : $(OBJ_CDF_DUMMY) $(THISFILE)

$(OBJ_CDF_DUMMY) : jules_netcdf_dummy.f90 $(THISFILE)
	$(FC) -c $(FFLAGS) jules_netcdf_dummy.f90
clean : 
	@$(RM) $(OBJ_CDF_DUMMY) $(wildcard ./*$(MOD_FSUF))
