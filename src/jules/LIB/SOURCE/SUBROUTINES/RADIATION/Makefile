include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(albpft.o)        $(ARC)(albsnow.o)     $(ARC)(ftsa.o)            \
      $(ARC)(tile_albedo.o)   $(ARC)(canyonalb.o)   $(ARC)(matinv.o)          \
      $(ARC)(urbanemis.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(albpft.o): $(ARC)(nstypes.o) $(ARC)(pftparm.o) $(THISFILE)
$(ARC)(albsnow.o): $(ARC)(nstypes.o) $(ARC)(rad_param.o) $(THISFILE)
$(ARC)(ftsa.o): $(ARC)(c_0_dg_c.o) $(THISFILE)
$(ARC)(tile_albedo.o):  $(ARC)(nstypes.o)   $(ARC)(nvegparm.o)                \
                        $(ARC)(pftparm.o)   $(ARC)(c_0_dg_c.o)                \
                        $(ARC)(rad_param.o) $(ARC)(snow_param.o)              \
                        $(ARC)(switches.o)  $(ARC)(switches_urban.o)          \
                        $(ARC)(p_s_parms.o) $(ARC)(jules_mod.o)               \
                        $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

## End of file.
